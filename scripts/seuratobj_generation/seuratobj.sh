library(SingleCellExperiment)
library(Seurat)
library(dplyr)

# 1. Read in your RDS
rds_path <- "/data/files/merged.rds"
sce <- readRDS(rds_path)

# 2. Inspect the SCE
cat("Class:", class(sce), "\n")
cat("Dimensions (cells × genes):", dim(sce), "\n")
cat("Assays available:\n"); print(assayNames(sce))
cat("Column metadata fields:\n"); print(colnames(colData(sce)))

# 3. Convert to Seurat
seu <- as.Seurat(
  sce,
  counts = "counts",
  data   = "counts"
)

# 4. Quick Seurat check
cat("\nSeurat assays:\n"); print(Assays(seu))
cat("Seurat dims (cells × features):", dim(seu), "\n")
cat("Seurat metadata columns:\n"); print(colnames(seu@meta.data))
class(seu)

# 5. Rename the two columns
seu@meta.data <- seu@meta.data %>%
  rename(
    Sample_ID = sample,
    cell_type = celltypist.Human_AdultAged_Hippocampus
  )

# 6. Verify
cat("\nAfter renaming, metadata columns:\n")
print(colnames(seu@meta.data))

# save as an RDS
saveRDS(seu, file = "/data/files/merged_seurat.rds")