# process_seurat.R
library(Seurat)
library(dplyr)
library(readr)

# --- INPUT FILES ---
seurat_path <- "/data/files/merged_seuratobj.rds"
assignments_path <- "/data/files/demuxlet_assignments.txt"
output_path <- "/data/files/matched_seurat_demuxlet.rds"

seurat_obj <- readRDS(seurat_path)
assignments <- read.table(assignments_path, header = TRUE, sep = "", stringsAsFactors = FALSE, fill = TRUE, quote = "\"", check.names = FALSE)

colnames(assignments) <- gsub("\"", "", colnames(assignments))
head(colnames(seurat_obj))
head(assignments)

# Clean Seurat barcodes: remove prefix before last '_'
seurat_obj$clean_barcode <- sub(".*_", "", colnames(seurat_obj))

# Clean assignment barcodes: remove prefix before '_'
assignments$Clean_Barcode <- sub("^[^_]+_", "", assignments$Barcode)

# Match barcodes
intersecting_barcodes <- intersect(seurat_obj$clean_barcode, assignments$Clean_Barcode)
cat("Matched", length(intersecting_barcodes), "barcodes\n")

# Filter Seurat object to matched cells only
cat("Filtering Seurat object to common cells...\n")
cells_to_keep <- seurat_obj$clean_barcode %in% intersecting_barcodes
seurat_obj <- seurat_obj[, cells_to_keep]

# Create named vector for metadata
cat("Preparing Individual_ID vector...\n")
barcode_to_individual <- setNames(assignments$Individual_Assignment, assignments$Clean_Barcode)
individual_ids <- barcode_to_individual[seurat_obj$clean_barcode]

# Name it with actual Seurat barcodes
names(individual_ids) <- colnames(seurat_obj)

# Add Individual_ID metadata (fast and safe)
cat("Adding metadata...\n")
seurat_obj <- AddMetaData(seurat_obj, metadata = individual_ids, col.name = "Individual_ID")

# Sanity check: make sure there are no NAs
na_count <- sum(is.na(seurat_obj$Individual_ID))
cat("NA entries in Individual_ID:", na_count, "\n")
if (na_count > 0) stop("Some cells are missing individual assignments.")

# Drop the helper column
seurat_obj$clean_barcode <- NULL

##Verifying final Seurat object before saving

# 1. Preview metadata
cat("First few rows of metadata:\n")
print(head(seurat_obj@meta.data[, c("Individual_ID"), drop = FALSE]))

# 2. Check that number of Individual_IDs matches number of cells
stopifnot(length(seurat_obj$Individual_ID) == ncol(seurat_obj))

# 3. Check for any NA values (should be 0)
na_count <- sum(is.na(seurat_obj$Individual_ID))
cat("Number of NA values in Individual_ID:", na_count, "\n")
if (na_count > 0) stop("Error: NA values found in Individual_ID. Something went wrong.")

# 4. Show number of unique individuals
cat("Number of unique individuals:\n")
print(length(unique(seurat_obj$Individual_ID)))

# 5. Show a frequency table of first few individuals
cat("Counts per individual (top 10):\n")
print(head(table(seurat_obj$Individual_ID), 10))

colnames(seurat_obj@meta.data)

saveRDS(seurat_obj, output_path)
cat("Final Seurat object saved successfully and is ready for eQTL analysis.\n")