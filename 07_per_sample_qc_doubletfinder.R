# Step 1: per-sample QC filter + DoubletFinder for Hou and Absinta.
#
# DoubletFinder must run on each sample individually, before any merging
# (per-sample: load -> QC filter -> normalize/PCA/cluster -> DoubletFinder
# -> remove doublets -> save). Order verified against DoubletFinder's own
# README: quality filtering -> normalization -> PCA -> DoubletFinder, and
# "do not apply DoubletFinder to aggregated data representing multiple
# distinct samples."
#
# Run this AFTER pulling the STARsolo filtered matrices down locally (see
# the rclone commands in the accompanying instructions) -- reading many
# small matrix files over a live Box mount is much slower than local disk.

# ---- One-time package install (skip if already installed) ----
# install.packages(c("Seurat", "Matrix", "dplyr"))
# install.packages("devtools")
# devtools::install_github("chris-mcginnis-ucsf/DoubletFinder")

library(Seurat)
library(DoubletFinder)
library(dplyr)

set.seed(1)

# ---- Reusable per-dataset processing function ----
#
# mito_pattern: gene-symbol regex for mitochondrial genes -- mouse uses
# lowercase "mt-", human uses uppercase "MT-". Case matters, don't reuse
# the same pattern across species.
#
# max_mito_pct defaults to 5, appropriate for single-nucleus data (low
# cytoplasmic contamination expected) -- reconsider if either dataset
# turns out to be whole-cell rather than single-nucleus.
process_dataset <- function(dataset_name,
                             input_root,
                             output_root,
                             mito_pattern,
                             box_remote,
                             min_features = 200,
                             max_features = 6000,
                             min_counts = 500,
                             max_mito_pct = 5,
                             expected_doublet_rate_per_1000 = 0.008,
                             n_pcs = 15,
                             cluster_resolution = 0.6) {

  # Expand "~" to an absolute path now -- shQuote() wraps paths in single
  # quotes for the shell-based rclone push below, and POSIX shells do NOT
  # expand "~" inside single quotes, so an unexpanded "~/..." path gets
  # passed through literally and mis-resolved (confirmed by an actual
  # failure: rclone saw "/home/exouser/~/scratch/..." as the path).
  output_root <- path.expand(output_root)
  input_root <- path.expand(input_root)

  dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
  sample_dirs <- list.dirs(input_root, recursive = FALSE, full.names = TRUE)

  for (sample_dir in sample_dirs) {
    sample_id <- basename(sample_dir)

    out_file <- file.path(output_root, paste0(sample_id, "_clean.rds"))
    if (file.exists(out_file)) {
      message(sample_id, ": already processed (", out_file, " exists) -- skipping")
      next
    }

    message("=== ", dataset_name, ": ", sample_id, " ===")

    # Wrap this sample's processing so any single failure -- a corrupt
    # file, a dropped Box mount, DoubletFinder erroring on a weird cluster
    # structure -- doesn't kill the whole run the way an uncaught error
    # from source() would. next/break aren't affected by tryCatch, so the
    # "next" calls below still correctly move on to the next sample.
    tryCatch({

    matrix_dir <- file.path(sample_dir, "Solo.out", "Gene", "filtered")
    if (!dir.exists(matrix_dir)) {
      warning("No filtered matrix found for ", sample_id, " at ", matrix_dir, " -- skipping")
      next
    }

    # STARsolo writes these uncompressed by default, but Read10X() requires
    # them gzipped -- compress in place if this sample hasn't been already.
    if (file.exists(file.path(matrix_dir, "matrix.mtx")) &&
        !file.exists(file.path(matrix_dir, "matrix.mtx.gz"))) {
      system2("gzip", shQuote(Sys.glob(file.path(matrix_dir, "*.tsv"))))
      system2("gzip", shQuote(file.path(matrix_dir, "matrix.mtx")))
    }

    # ---- Load ----
    counts <- Read10X(data.dir = matrix_dir)
    seu <- CreateSeuratObject(counts = counts, project = sample_id, min.cells = 3)
    seu$sample_id <- sample_id
    seu$dataset <- dataset_name

    # ---- QC filter ----
    # Count passing cells BEFORE calling subset() -- Seurat's subset()
    # throws a hard error ("No cells found") when zero cells match the
    # filter, which would otherwise crash the whole script instead of
    # letting us skip this sample gracefully like we do for the
    # few-but-not-zero case.
    seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = mito_pattern)
    qc_pass <- seu$nFeature_RNA > min_features &
      seu$nFeature_RNA < max_features &
      seu$nCount_RNA > min_counts &
      seu$percent.mt < max_mito_pct

    if (sum(qc_pass) < 50) {
      warning(sample_id, ": only ", sum(qc_pass), " cells survived QC (need >= 50) -- skipping")
      next
    }

    seu <- subset(seu, cells = colnames(seu)[qc_pass])

    # ---- Normalize / PCA / cluster (DoubletFinder needs both) ----
    seu <- NormalizeData(seu)
    seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)
    seu <- ScaleData(seu)
    seu <- RunPCA(seu, npcs = n_pcs)
    seu <- FindNeighbors(seu, dims = 1:n_pcs)
    seu <- FindClusters(seu, resolution = cluster_resolution)

    # ---- DoubletFinder: sweep for the optimal pK ----
    sweep.res <- paramSweep(seu, PCs = 1:n_pcs, sct = FALSE)
    sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
    bcmvn <- find.pK(sweep.stats)
    best_pK <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))

    # ---- Expected doublet count, scaled to this sample's recovered
    # cells (10x's published multiplet-rate table: ~0.8% per 1,000 cells) ----
    n_cells <- ncol(seu)
    doublet_rate <- (n_cells / 1000) * expected_doublet_rate_per_1000
    nExp_poi <- round(doublet_rate * n_cells)
    homotypic.prop <- modelHomotypic(seu$seurat_clusters)
    nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))

    # ---- Run DoubletFinder (homotypic-adjusted) ----
    seu <- doubletFinder(
      seu, PCs = 1:n_pcs, pN = 0.25, pK = best_pK,
      nExp = nExp_poi.adj, reuse.pANN = NULL, sct = FALSE
    )

    # Classification column name is dynamic (embeds pN/pK/nExp) -- find
    # whichever one doubletFinder() just added.
    df_col <- grep("^DF.classifications", colnames(seu@meta.data), value = TRUE)
    df_col <- df_col[length(df_col)]

    n_doublets <- sum(seu@meta.data[[df_col]] == "Doublet")
    message(sample_id, ": ", n_doublets, " / ", n_cells, " cells flagged as doublets")

    # ---- Remove doublets ----
    seu <- subset(seu, cells = rownames(seu@meta.data)[seu@meta.data[[df_col]] == "Singlet"])

    # ---- Save cleaned per-sample object for the later merge step ----
    saveRDS(seu, out_file)
    message(sample_id, ": saved ", ncol(seu), " cells after doublet removal")

    # ---- Push to Box -- nothing should live only on Jetstream's disk ----
    # NOTE: system2() joins command+args into one string and runs it through
    # a shell, so any argument containing spaces (our Box path has "Lab
    # members" in it) MUST be shQuote()'d or the shell word-splits it into
    # multiple arguments (confirmed by an actual failure: rclone received
    # "HUANG_LAB_CENTRAL/Lab" and "members/..." as two separate args).
    push_status <- system2("rclone", c("copy", shQuote(out_file), shQuote(box_remote)))
    if (push_status != 0) {
      warning(sample_id, ": rclone push to Box FAILED (exit status ", push_status, ")")
    } else {
      message(sample_id, ": pushed to Box")
    }

    }, error = function(e) {
      warning(sample_id, ": FAILED with error -- ", conditionMessage(e), " -- skipping and continuing to next sample")
    })
  }
}

BOX_ROOT <- "box:HUANG_LAB_CENTRAL/Lab members/MeghanB./AbsintaHouMelchor_reanalysis/2026_data/Jetstream_moving/01_qc_doubletfinder"

# ---- Hou (mouse, 10x 5' v1) ----
process_dataset(
  dataset_name = "Hou",
  input_root   = "~/scratch/seurat_input/hou",
  output_root  = "~/scratch/seurat_clean/hou",
  mito_pattern = "^mt-",
  box_remote   = paste0(BOX_ROOT, "/Hou")
)

# ---- Absinta (human, 10x 3' v3) ----
process_dataset(
  dataset_name = "Absinta",
  input_root   = "~/scratch/seurat_input/absinta",
  output_root  = "~/scratch/seurat_clean/absinta",
  box_remote   = paste0(BOX_ROOT, "/Absinta"),
  mito_pattern = "^MT-"
)

# ---- aboelnourCPZ (mouse, 10x 3' v3) ----
process_dataset(
  dataset_name = "aboelnourCPZ",
  input_root   = "~/scratch/seurat_input/aboelnourCPZ",
  output_root  = "~/scratch/seurat_clean/aboelnourCPZ",
  box_remote   = paste0(BOX_ROOT, "/aboelnourCPZ"),
  mito_pattern = "^mt-"
)

# ---- aboelnourLPC (mouse, 10x 3' v3.1) ----
process_dataset(
  dataset_name = "aboelnourLPC",
  input_root   = "~/scratch/seurat_input/aboelnourLPC",
  output_root  = "~/scratch/seurat_clean/aboelnourLPC",
  box_remote   = paste0(BOX_ROOT, "/aboelnourLPC"),
  mito_pattern = "^mt-"
)

# ---- Shen (mouse, 10x 3' v2) ----
process_dataset(
  dataset_name = "Shen",
  input_root   = "~/scratch/seurat_input/shen",
  output_root  = "~/scratch/seurat_clean/shen",
  box_remote   = paste0(BOX_ROOT, "/Shen"),
  mito_pattern = "^mt-"
)

# ---- Melchor (mouse, 10x 3' v3.1) ----
process_dataset(
  dataset_name = "Melchor",
  input_root   = "~/scratch/seurat_input/melchor",
  output_root  = "~/scratch/seurat_clean/melchor",
  box_remote   = paste0(BOX_ROOT, "/Melchor"),
  mito_pattern = "^mt-"
)

# ---- Schirmer (human, 10x 3' v2/v3 auto-detected per sample) ----
process_dataset(
  dataset_name = "Schirmer",
  input_root   = "~/scratch/seurat_input/schirmer",
  output_root  = "~/scratch/seurat_clean/schirmer",
  box_remote   = paste0(BOX_ROOT, "/Schirmer"),
  mito_pattern = "^MT-"
)

# ---- Wheeler mouse EAE (mouse, Drop-seq) ----
# QC matches Wheeler et al.'s own reported Methods exactly, rather than
# our generic cross-dataset default: their Drop-seq software was run with
# "min_num_genes_per_cell = 500" for the B6 EAE studies specifically, and
# no UMI-count or mitochondrial cutoff is reported at all. Confirmed this
# matters empirically -- our generic 500-UMI/5%-mito filter was found to
# discard hundreds of EmptyDrops_CR-called real cells per sample in the
# human Drop-seq cohort that Wheeler's own (lenient) threshold would keep
# (see conversation/investigation on the human Drop-seq QC mismatch).
process_dataset(
  dataset_name = "Wheeler_mouse",
  input_root   = "~/scratch/seurat_input/wheeler_mouse",
  output_root  = "~/scratch/seurat_clean/wheeler_mouse",
  box_remote   = paste0(BOX_ROOT, "/Wheeler_mouse"),
  mito_pattern = "^mt-",
  min_features = 500,
  max_features = Inf,
  min_counts   = 0,
  max_mito_pct = 100
)

# ---- Wheeler human 10x (human, 10x 3' v2 -- Control1-3) ----
process_dataset(
  dataset_name = "Wheeler_human_10x",
  input_root   = "~/scratch/seurat_input/wheeler_human_10x",
  output_root  = "~/scratch/seurat_clean/wheeler_human_10x",
  box_remote   = paste0(BOX_ROOT, "/Wheeler_human_10x"),
  mito_pattern = "^MT-"
)

# ---- Wheeler human Drop-seq (human, Drop-seq -- MS patients + Control4/5) ----
# QC matches Wheeler et al.'s own reported Methods: "min_num_genes_per_cell
# = 200" for the human samples, no UMI-count or mitochondrial cutoff
# reported. Confirmed empirically that our generic 500-UMI/5%-mito filter
# was discarding hundreds of EmptyDrops_CR-called real cells per sample
# that this paper-matched threshold recovers -- see check_qc_threshold_
# breakdown.R results. 7 of 18 previously-failed samples cross the
# 50-cell minimum under this threshold; the rest (including MSpatient3-1
# and 3-3, which have 0 cells passing even this lenient filter) are
# genuine low-yield/failed libraries, not a QC artifact.
process_dataset(
  dataset_name = "Wheeler_human_dropseq",
  input_root   = "~/scratch/seurat_input/wheeler_human_dropseq",
  output_root  = "~/scratch/seurat_clean/wheeler_human_dropseq",
  box_remote   = paste0(BOX_ROOT, "/Wheeler_human_dropseq"),
  mito_pattern = "^MT-",
  min_features = 200,
  max_features = Inf,
  min_counts   = 0,
  max_mito_pct = 100
)
