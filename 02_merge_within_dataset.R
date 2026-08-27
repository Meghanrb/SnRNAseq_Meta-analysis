# Step 2: merge the per-sample cleaned Seurat objects within each dataset
# into one combined object per dataset. Absinta gets a qc_flag column
# marking 4 samples with confirmed low cell yield (see investigation:
# GSM5903109, GSM5470486, GSM5470490, GSM5903107 -- all MS lesion tissue,
# retained <35% of cells after standard QC, likely due to genuine
# nuclei-prep difficulty in necrotic/gliotic lesion tissue) -- kept in
# rather than excluded, but flagged for downstream filtering if needed.

library(Seurat)
library(dplyr)

merge_dataset <- function(dataset_name, input_root, output_file, low_yield_gsms = character()) {
  input_root <- path.expand(input_root)
  output_file <- path.expand(output_file)

  if (file.exists(output_file)) {
    message(dataset_name, ": already merged (", output_file, " exists) -- skipping")
    return(invisible(NULL))
  }

  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)

  rds_files <- list.files(input_root, pattern = "_clean\\.rds$", full.names = TRUE)
  message(dataset_name, ": merging ", length(rds_files), " samples")

  seurat_list <- lapply(rds_files, readRDS)

  merged <- if (length(seurat_list) == 1) {
    seurat_list[[1]]
  } else {
    merge(seurat_list[[1]], y = seurat_list[-1])
  }

  # Seurat v5's merge() keeps each input object as a SEPARATE layer
  # instead of one real combined counts matrix -- anything that reads
  # expression data (not just cell-level metadata like ncol()/sample_id,
  # which merge() handles fine on its own) will silently only see the
  # FIRST sample's layer otherwise. Same bug found and fixed in the
  # Wheeler patient-pooling script; must run here too before this object
  # is used for anything beyond counting cells.
  if (length(seurat_list) > 1) {
    merged <- JoinLayers(merged)
  }

  merged$qc_flag <- ifelse(merged$sample_id %in% low_yield_gsms, "low_yield", "normal")

  message(dataset_name, ": merged object has ", ncol(merged), " total cells")
  if (length(low_yield_gsms) > 0) {
    message(dataset_name, ": ", sum(merged$qc_flag == "low_yield"), " cells flagged low_yield across ", length(low_yield_gsms), " samples")
  }

  saveRDS(merged, path.expand(output_file))
  message(dataset_name, ": saved to ", output_file)

  push_status <- system2("rclone", c("copy", shQuote(path.expand(output_file)), shQuote(BOX_MERGE_REMOTE)))
  if (push_status != 0) {
    warning(dataset_name, ": rclone push to Box FAILED (exit status ", push_status, ")")
  } else {
    message(dataset_name, ": pushed to Box")
  }

  merged
}

BOX_MERGE_REMOTE <- "box:HUANG_LAB_CENTRAL/Lab members/MeghanB./AbsintaHouMelchor_reanalysis/2026_data/Jetstream_moving/02_merged"

hou_merged <- merge_dataset(
  dataset_name = "Hou",
  input_root   = "~/scratch/seurat_clean/hou",
  output_file  = "~/scratch/seurat_merged/hou_merged.rds"
)

absinta_merged <- merge_dataset(
  dataset_name = "Absinta",
  input_root   = "~/scratch/seurat_clean/absinta",
  output_file  = "~/scratch/seurat_merged/absinta_merged.rds",
  low_yield_gsms = c("GSM5903109", "GSM5470486", "GSM5470490", "GSM5903107")
)

aboelnourCPZ_merged <- merge_dataset(
  dataset_name = "aboelnourCPZ",
  input_root   = "~/scratch/seurat_clean/aboelnourCPZ",
  output_file  = "~/scratch/seurat_merged/aboelnourCPZ_merged.rds"
)

aboelnourLPC_merged <- merge_dataset(
  dataset_name = "aboelnourLPC",
  input_root   = "~/scratch/seurat_clean/aboelnourLPC",
  output_file  = "~/scratch/seurat_merged/aboelnourLPC_merged.rds"
)

shen_merged <- merge_dataset(
  dataset_name = "Shen",
  input_root   = "~/scratch/seurat_clean/shen",
  output_file  = "~/scratch/seurat_merged/shen_merged.rds"
)

melchor_merged <- merge_dataset(
  dataset_name = "Melchor",
  input_root   = "~/scratch/seurat_clean/melchor",
  output_file  = "~/scratch/seurat_merged/melchor_merged.rds"
)

schirmer_merged <- merge_dataset(
  dataset_name = "Schirmer",
  input_root   = "~/scratch/seurat_clean/schirmer",
  output_file  = "~/scratch/seurat_merged/schirmer_merged.rds"
)

wheeler_mouse_merged <- merge_dataset(
  dataset_name = "Wheeler_mouse",
  input_root   = "~/scratch/seurat_clean/wheeler_mouse",
  output_file  = "~/scratch/seurat_merged/wheeler_mouse_merged.rds"
)

wheeler_human_10x_merged <- merge_dataset(
  dataset_name = "Wheeler_human_10x",
  input_root   = "~/scratch/seurat_clean/wheeler_human_10x",
  output_file  = "~/scratch/seurat_merged/wheeler_human_10x_merged.rds"
)

# Note: input_root here must contain the CONSOLIDATED set of files
# (Control4/5 + the 4 patient-pooled MSpatient objects), not the original
# per-region files -- see the cleanup step that swapped those in before
# this script was updated.
wheeler_human_dropseq_merged <- merge_dataset(
  dataset_name = "Wheeler_human_dropseq",
  input_root   = "~/scratch/seurat_clean/wheeler_human_dropseq",
  output_file  = "~/scratch/seurat_merged/wheeler_human_dropseq_merged.rds"
)
