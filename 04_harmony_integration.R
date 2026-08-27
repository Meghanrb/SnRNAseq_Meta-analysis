# Step 4: Harmony batch integration on the final merged object.
#
# Normalization/scaling/PCA are redone FRESH here on the fully merged,
# joined-layer object -- not reused from any per-sample or per-dataset
# processing done earlier in the pipeline. Batch-aware integration
# should start from a common baseline computed across all cells
# together, not stitch together normalizations computed independently
# per sample/dataset.
#
# Harmony is run with group.by.vars = "dataset" -- correcting for each
# of the 10 source studies as its own batch. This is the starting
# choice; if residual batch structure shows up in the resulting UMAP
# (datasets still separating instead of mixing within shared cell
# types), the next thing to try is adding chemistry as a second
# correction variable rather than assuming dataset alone is sufficient.

library(Seurat)
library(harmony)

set.seed(1)

IN_FILE  <- path.expand("~/scratch/final_merged/all_datasets_merged.rds")
OUT_DIR  <- path.expand("~/scratch/final_merged")
OUT_FILE <- file.path(OUT_DIR, "all_datasets_harmony.rds")
BOX_FINAL <- "box:HUANG_LAB_CENTRAL/Lab members/MeghanB./AbsintaHouMelchor_reanalysis/2026_data/Jetstream_moving/04_final_merged"

message("Loading merged object...")
merged_all <- readRDS(IN_FILE)
message(ncol(merged_all), " cells, ", length(unique(merged_all$dataset)), " datasets")

message("Normalizing / finding variable features / scaling (fresh, on the full merged object)...")
merged_all <- NormalizeData(merged_all)
merged_all <- FindVariableFeatures(merged_all, selection.method = "vst", nfeatures = 2000)
merged_all <- ScaleData(merged_all)

message("Running PCA...")
merged_all <- RunPCA(merged_all, npcs = 30)

message("Running Harmony, correcting for 'dataset'...")
merged_all <- RunHarmony(merged_all, group.by.vars = "dataset")

message("Running UMAP + clustering on the Harmony-corrected embedding (to actually see whether it worked)...")
merged_all <- RunUMAP(merged_all, reduction = "harmony", dims = 1:30)
merged_all <- FindNeighbors(merged_all, reduction = "harmony", dims = 1:30)
merged_all <- FindClusters(merged_all, resolution = 0.6)

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
saveRDS(merged_all, OUT_FILE)
message("Saved to ", OUT_FILE)

push_status <- system2("rclone", c("copy", shQuote(OUT_FILE), shQuote(BOX_FINAL)))
if (push_status != 0) {
  warning("rclone push to Box FAILED (exit status ", push_status, ")")
} else {
  message("Pushed to Box under 04_final_merged/")
}

message("Done. Check DimPlot(merged_all, reduction = 'umap', group.by = 'dataset') to see whether datasets are mixing (good) or still separating (needs more correction).")
