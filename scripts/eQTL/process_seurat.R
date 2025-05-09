# process_seurat.R
library(Seurat)
library(dplyr)
library(readr)

# --- INPUT FILES ---
seurat_path <- "/rds/general/project/puklandmarkproject/live/Users/Parisa/parisa_eqtl/temporary_seurat_object.rds"
assignments_path <- "/rds/general/user/pr422/projects/puklandmarkproject/ephemeral/Parisa/all_assignments.tsv"
output_path <- "/rds/general/project/puklandmarkproject/live/Users/Parisa/parisa_eqtl/updated_seurat_with_individual_id.rds"

cat("Loading Seurat object...\n")
seurat_obj <- readRDS(seurat_path)

cat("Reading demultiplexing assignments...\n")
assignments <- read_tsv(assignments_path, show_col_types = FALSE)

# ---- CLEAN BARCODES FOR MATCHING ----
cat("Cleaning barcodes for matching...\n")

# Clean Seurat cell barcodes (remove everything before last "_")
seurat_obj$clean_barcode <- sub(".*_", "", colnames(seurat_obj))

# Clean assignment barcodes (remove prefix like 'S10A_')
assignments$Clean_Barcode <- sub("^[^_]+_", "", assignments$Barcode)

# ---- MATCHING ----
cat("Finding overlapping barcodes...\n")
intersecting_barcodes <- intersect(seurat_obj$clean_barcode, assignments$Clean_Barcode)
cat("Matched", length(intersecting_barcodes), "barcodes\n")

# ---- FILTER TO COMMON CELLS ----
cat("Filtering Seurat object...\n")
cells_to_keep <- seurat_obj$clean_barcode %in% intersecting_barcodes
seurat_obj <- seurat_obj[, cells_to_keep]

# ---- ADD Individual_ID COLUMN ----
cat("Adding Individual_ID to metadata...\n")
barcode_to_individual <- setNames(assignments$Demuxalot_Individual_Assignment,
                                  assignments$Clean_Barcode)

individual_ids <- barcode_to_individual[seurat_obj$clean_barcode]
names(individual_ids) <- colnames(seurat_obj)
seurat_obj$Individual_ID <- individual_ids

# ---- CLEANUP ----
seurat_obj$clean_barcode <- NULL

# ---- SAVE OUTPUT ----
cat("Saving updated Seurat object to:\n", output_path, "\n")
saveRDS(seurat_obj, file = output_path)
cat("Done! Ready for eQTL pipeline.\n")