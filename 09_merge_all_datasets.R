# Step 3: final cross-dataset merge. Combines all 10 datasets into one
# object, all sharing human gene symbols:
#   - 6 mouse datasets: pulled from OrthologAL_converted output, using
#     the "RNA_ortho" assay (human-ortholog symbols) -- the original
#     "RNA" assay (mouse symbols) is discarded, not needed downstream.
#   - 4 human datasets: pulled from the within-dataset merge output,
#     using "RNA" as-is (already human symbols, no conversion needed).
#
# Every source object gets its relevant assay explicitly JoinLayers()'d
# before extraction, regardless of whether it's already joined -- this
# is a lesson learned the hard way earlier in this pipeline (Seurat v5's
# merge() silently leaves per-sample layers unjoined, and trusting that
# an upstream step already fixed it caused real, hard-to-detect bugs).

library(Seurat)
library(dplyr)

set.seed(1)

BOX_ORTHO  <- "box:HUANG_LAB_CENTRAL/Lab members/MeghanB./AbsintaHouMelchor_reanalysis/2026_data/Jetstream_moving/03_orthologAL_converted"
BOX_MERGED <- "box:HUANG_LAB_CENTRAL/Lab members/MeghanB./AbsintaHouMelchor_reanalysis/2026_data/Jetstream_moving/02_merged"
BOX_FINAL  <- "box:HUANG_LAB_CENTRAL/Lab members/MeghanB./AbsintaHouMelchor_reanalysis/2026_data/Jetstream_moving/04_final_merged"

LOCAL_ROOT <- path.expand("~/scratch/final_merge_input")
dir.create(LOCAL_ROOT, recursive = TRUE, showWarnings = FALSE)

mouse_datasets <- list(
  Hou           = "OrthologAL_ hou_merged.rds",
  aboelnourCPZ  = "OrthologAL_ aboelnourCPZ_merged.rds",
  aboelnourLPC  = "OrthologAL_ aboelnourLPC_merged.rds",
  Shen          = "OrthologAL_ shen_merged.rds",
  Melchor       = "OrthologAL_ melchor_merged.rds",
  Wheeler_mouse = "OrthologAL_ wheeler_mouse_merged.rds"
)

human_datasets <- list(
  Absinta               = "absinta_merged.rds",
  Schirmer              = "schirmer_merged.rds",
  Wheeler_human_10x     = "wheeler_human_10x_merged.rds",
  Wheeler_human_dropseq = "wheeler_human_dropseq_merged.rds"
)

pull_if_missing <- function(box_remote, box_filename) {
  local_file <- file.path(LOCAL_ROOT, box_filename)
  if (!file.exists(local_file)) {
    message("Pulling ", box_filename, " from Box...")
    system2("rclone", c("copy", shQuote(paste0(box_remote, "/", box_filename)), shQuote(LOCAL_ROOT)))
  }
  local_file
}

prep_mouse_object <- function(dataset_name, box_filename) {
  local_file <- pull_if_missing(BOX_ORTHO, box_filename)
  seu <- readRDS(local_file)

  if (length(Layers(seu[["RNA_ortho"]])) > 1) {
    seu[["RNA_ortho"]] <- JoinLayers(seu[["RNA_ortho"]])
  }

  ortho_counts <- GetAssayData(seu, assay = "RNA_ortho", layer = "counts")
  meta <- seu@meta.data
  meta$original_species <- "mouse"

  new_seu <- CreateSeuratObject(counts = ortho_counts, meta.data = meta)
  new_seu$dataset <- dataset_name
  message(dataset_name, ": ", ncol(new_seu), " cells prepared (ortholog-converted)")
  new_seu
}

prep_human_object <- function(dataset_name, box_filename) {
  local_file <- pull_if_missing(BOX_MERGED, box_filename)
  seu <- readRDS(local_file)

  if (length(Layers(seu[["RNA"]])) > 1) {
    seu[["RNA"]] <- JoinLayers(seu[["RNA"]])
  }

  counts <- GetAssayData(seu, assay = "RNA", layer = "counts")
  meta <- seu@meta.data
  meta$original_species <- "human"

  new_seu <- CreateSeuratObject(counts = counts, meta.data = meta)
  new_seu$dataset <- dataset_name
  message(dataset_name, ": ", ncol(new_seu), " cells prepared (native human)")
  new_seu
}

message("=== Preparing 6 mouse (ortholog-converted) datasets ===")
mouse_objs <- Map(prep_mouse_object, names(mouse_datasets), mouse_datasets)

message("=== Preparing 4 human datasets ===")
human_objs <- Map(prep_human_object, names(human_datasets), human_datasets)

all_objs <- c(mouse_objs, human_objs)

message("=== Merging all ", length(all_objs), " datasets ===")
merged_all <- merge(all_objs[[1]], y = all_objs[-1], add.cell.ids = names(all_objs))

message("Joining layers on the final merged object...")
merged_all <- JoinLayers(merged_all)

message("Final merged object: ", ncol(merged_all), " total cells across ", length(all_objs), " datasets")
print(table(merged_all$dataset))

OUT_DIR <- path.expand("~/scratch/final_merged")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
out_file <- file.path(OUT_DIR, "all_datasets_merged.rds")
saveRDS(merged_all, out_file)
message("Saved to ", out_file)

push_status <- system2("rclone", c("copy", shQuote(out_file), shQuote(BOX_FINAL)))
if (push_status != 0) {
  warning("rclone push to Box FAILED (exit status ", push_status, ")")
} else {
  message("Pushed to Box under 04_final_merged/")
}
