# Wheeler human Drop-seq MS patients: pool each patient's separate
# tissue-region libraries into one object BEFORE running DoubletFinder,
# matching how Wheeler et al.'s own analysis unit was the patient (n=4
# MS patients), not the individual region (n=29 regional libraries).
#
# Regions are still aligned/cell-called separately (already done via
# STARsolo -- can't merge raw FASTQs across regions since Drop-seq
# barcodes are only unique within one physical run) and QC-filtered
# individually using Wheeler's own reported threshold (>=200 genes/cell,
# no UMI/mito cutoff). Only AFTER that per-region filtering do we merge
# a patient's regions together, using add.cell.ids so every cell stays
# traceable to its region and no barcode collisions are possible across
# genuinely different physical runs.
#
# Does NOT touch Control4/Control5 (already single-library-per-patient,
# no pooling needed) or the mouse EAE dataset (replicates are separate
# mice, not regions of one mouse).

library(Seurat)
library(DoubletFinder)
library(dplyr)

set.seed(1)

INPUT_ROOT <- path.expand("~/scratch/seurat_input/wheeler_human_dropseq")
OUTPUT_ROOT <- path.expand("~/scratch/seurat_clean/wheeler_human_dropseq_patient_pooled")
BOX_REMOTE <- "box:HUANG_LAB_CENTRAL/Lab members/MeghanB./AbsintaHouMelchor_reanalysis/2026_data/Jetstream_moving/01_qc_doubletfinder/Wheeler_human_dropseq_patient_pooled"
MITO_PATTERN <- "^MT-"
MIN_FEATURES <- 200  # Wheeler's own reported threshold for human samples

dir.create(OUTPUT_ROOT, recursive = TRUE, showWarnings = FALSE)

PATIENT_GROUPS <- list(
  MSpatient1 = c("MSpatient1-1", "MSpatient1-2", "MSpatient1-3", "MSpatient1-4", "MSpatient1-5", "MSpatient1-6", "MSpatient1-7"),
  MSpatient2 = c("MSpatient2-1", "MSpatient2-2", "MSpatient2-3", "MSpatient2-4", "MSpatient2-5", "MSpatient2-6", "MSpatient2-7", "MSpatient2-8"),
  MSpatient3 = c("MSpatient3-1", "MSpatient3-2", "MSpatient3-3", "MSpatient3-4", "MSpatient3-5"),
  MSpatient4 = c("MSpatient4-1", "MSpatient4-2", "MSpatient4-3", "MSpatient4-4", "MSpatient4-5", "MSpatient4-6", "MSpatient4-7", "MSpatient4-8", "MSpatient4-9")
)

for (patient_id in names(PATIENT_GROUPS)) {
  out_file <- file.path(OUTPUT_ROOT, paste0(patient_id, "_clean.rds"))
  if (file.exists(out_file)) {
    message(patient_id, ": already processed (", out_file, " exists) -- skipping")
    next
  }

  message("=== ", patient_id, " ===")

  tryCatch({

  regions <- PATIENT_GROUPS[[patient_id]]
  region_objs <- list()

  for (region in regions) {
    matrix_dir <- file.path(INPUT_ROOT, region, "Solo.out", "Gene", "filtered")
    if (!dir.exists(matrix_dir)) {
      message(region, ": no filtered matrix found -- excluding from pool")
      next
    }

    if (file.exists(file.path(matrix_dir, "matrix.mtx")) &&
        !file.exists(file.path(matrix_dir, "matrix.mtx.gz"))) {
      system2("gzip", shQuote(Sys.glob(file.path(matrix_dir, "*.tsv"))))
      system2("gzip", shQuote(file.path(matrix_dir, "matrix.mtx")))
    }

    counts <- Read10X(data.dir = matrix_dir)
    seu_region <- CreateSeuratObject(counts = counts, project = region, min.cells = 3)
    seu_region[["percent.mt"]] <- PercentageFeatureSet(seu_region, pattern = MITO_PATTERN)

    qc_pass <- seu_region$nFeature_RNA > MIN_FEATURES
    message(region, ": ", sum(qc_pass), " / ", ncol(seu_region), " barcodes pass QC")

    if (sum(qc_pass) == 0) {
      next
    }
    seu_region <- subset(seu_region, cells = colnames(seu_region)[qc_pass])
    seu_region$region_id <- region
    region_objs[[region]] <- seu_region
  }

  if (length(region_objs) == 0) {
    warning(patient_id, ": no regions had any usable cells -- skipping entire patient")
    next
  }

  seu <- if (length(region_objs) == 1) {
    region_objs[[1]]
  } else {
    merge(region_objs[[1]], y = region_objs[-1], add.cell.ids = names(region_objs))
  }

  # Seurat v5's merge() keeps each original object as a SEPARATE layer
  # (e.g. "counts.MSpatient4-1", "counts.MSpatient4-2", ...) instead of
  # combining them into one real counts matrix. Anything downstream that
  # expects a single unified layer -- DoubletFinder included -- silently
  # only sees the FIRST region's cells otherwise, which either produces
  # quietly wrong (single-region) results or an opaque "subscript out of
  # bounds" error once cell metadata (all regions) and the counts matrix
  # (first region only) stop lining up. JoinLayers() must run before
  # anything else touches this object.
  if (length(region_objs) > 1) {
    seu <- JoinLayers(seu)
  }

  seu$sample_id <- patient_id
  seu$dataset <- "Wheeler_human_dropseq_patient_pooled"

  message(patient_id, ": pooled ", ncol(seu), " cells across ", length(region_objs), " region(s)")

  if (ncol(seu) < 50) {
    warning(patient_id, ": still fewer than 50 cells even after pooling all regions -- skipping")
    next
  }

  # ---- Normalize / PCA / cluster (DoubletFinder needs both) ----
  seu <- NormalizeData(seu)
  seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)
  seu <- ScaleData(seu)
  seu <- RunPCA(seu, npcs = 15)
  seu <- FindNeighbors(seu, dims = 1:15)
  seu <- FindClusters(seu, resolution = 0.6)

  # ---- DoubletFinder: sweep for the optimal pK ----
  sweep.res <- paramSweep(seu, PCs = 1:15, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  best_pK <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))

  n_cells <- ncol(seu)
  doublet_rate <- (n_cells / 1000) * 0.008
  nExp_poi <- round(doublet_rate * n_cells)
  homotypic.prop <- modelHomotypic(seu$seurat_clusters)
  nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))

  seu <- doubletFinder(
    seu, PCs = 1:15, pN = 0.25, pK = best_pK,
    nExp = nExp_poi.adj, reuse.pANN = NULL, sct = FALSE
  )

  df_col <- grep("^DF.classifications", colnames(seu@meta.data), value = TRUE)
  df_col <- df_col[length(df_col)]

  n_doublets <- sum(seu@meta.data[[df_col]] == "Doublet")
  message(patient_id, ": ", n_doublets, " / ", n_cells, " cells flagged as doublets")

  seu <- subset(seu, cells = rownames(seu@meta.data)[seu@meta.data[[df_col]] == "Singlet"])

  saveRDS(seu, out_file)
  message(patient_id, ": saved ", ncol(seu), " cells after doublet removal")

  push_status <- system2("rclone", c("copy", shQuote(out_file), shQuote(BOX_REMOTE)))
  if (push_status != 0) {
    warning(patient_id, ": rclone push to Box FAILED (exit status ", push_status, ")")
  } else {
    message(patient_id, ": pushed to Box")
  }

  }, error = function(e) {
    warning(patient_id, ": FAILED with error -- ", conditionMessage(e), " -- skipping and continuing to next patient")
  })
}

message("Done. Results are in Box under 01_qc_doubletfinder/Wheeler_human_dropseq_patient_pooled/, one file per MS patient.")
