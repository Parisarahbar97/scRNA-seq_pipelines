# Load necessary packages
library(Seurat)
library(SingleCellExperiment)

# Load your data
tmp <- readRDS("merged.rds")

# Check your assay names to confirm the main assay:
assayNames(tmp)  # You should see something like "X" or "counts"

# Using "X" assay for counts and data
seurat_obj <- CreateSeuratObject(counts = assay(tmp, "X"))

# Add metadata clearly
seurat_obj <- AddMetaData(seurat_obj, metadata = as.data.frame(colData(tmp)))

# Specifically add CellType clearly
seurat_obj$CellType <- colData(tmp)$celltypist.Human_AdultAged_Hippocampus

# Attach UMAP coordinates manually
umap_coords <- reducedDim(tmp, "X_scvi-global_umap")

# Set UMAP embedding inside Seurat object
seurat_obj[["umap"]] <- CreateDimReducObject(embeddings = umap_coords, key = "UMAP_", assay = DefaultAssay(seurat_obj))

# Normalize the data first (standard Seurat step)
seurat_obj <- NormalizeData(seurat_obj)


library(ggplot2)
library(ggrepel)
library(dplyr)

# Expanded marker gene list (yours + previous)
marker_genes <- c(
  # Microglia
  "CX3CR1", "IBA1", "P2RY12", "TMEM119", "ITGAM", "CSF1R", "HEXB", "SIGLEC11", "SALL1", "TLR2", "TLR4", "CD68", "CD45", "TREM2", "CCL4", "CCL3", "CTSS", "TYROBP", "CD83",
  # Astrocytes
  "GFAP", "ALDH1L1", "S100B", "AQP4", "SOX9", "CPE", "CLU", "ALDOC", "SLC1A3", "SLC1A2",
  # Oligodendrocytes
  "MOG", "MBP", "CNP", "PLP1", "MAG", "MAL", "CRYAB",
  # Neurons (general)
  "MAP2", "NEUN", "SYN1", "GAP43", "TUBB3",
  # Endothelial Cells
  "PECAM1", "VWF", "CD34", "CDH5", "CLDN5", "ITM2A", "BSG", "RSG5", "APOLD1",
  # OPCs
  "PDGFRA", "NG2", "SOX10", "OLIG2", "NKX2-2", "CNTN1", "TNR",
  # GABAergic Neurons
  "GAD1", "GAD2", "VGAT",
  # Glutamatergic Neurons
  "VGLUT1", "SLC17A7", "CAMK2A", "VGLUT2",
  # Cajal-Retzius Cells
  "RELN", "CALB2", "TP73",
  # Choroid Plexus
  "TTR", "AQP1", "OTX2"
)

# Calculate cluster centers (assuming plot_df has UMAP1, UMAP2, and CellType)
centers <- plot_df %>%
  group_by(CellType) %>%
  summarise(
    UMAP1 = median(UMAP1),
    UMAP2 = median(UMAP2),
    .groups = "drop"
  )

# Plot each marker gene and display (no saving)
for (gene in marker_genes) {
  
  if (gene %in% rownames(seurat_obj)) {
    
    # Create the FeaturePlot
    p <- FeaturePlot(seurat_obj, features = gene, reduction = "umap",
                     pt.size = 0.8,
                     min.cutoff = "q10",
                     cols = c("lightgrey", "red")) +
         ggtitle(paste("Expression of", gene)) +
         theme(plot.title = element_text(hjust = 0.5))
    
    # Add cluster titles
    p <- p + 
      geom_text_repel(
        data = centers,
        aes(x = UMAP1, y = UMAP2, label = CellType),
        inherit.aes = FALSE,
        color = "black",
        size = 3,
        fontface = "bold",
        max.overlaps = 50,
        show.legend = FALSE
      )
    
    # Just display it
    print(p)
    
  } else {
    print(paste("Gene", gene, "not found in the dataset."))
  }
}



####whith subsamplenig

library(Seurat)
library(ggplot2)
library(ggrepel)
library(dplyr)

# Expanded marker gene list
marker_genes <- c(
  # Microglia
  "CX3CR1", "IBA1", "P2RY12", "TMEM119", "ITGAM", "CSF1R", "HEXB", "SIGLEC11", "SALL1", "TLR2", "TLR4", "CD68", "CD45", "TREM2", "CCL4", "CCL3", "CTSS", "TYROBP", "CD83",
  # Astrocytes
  "GFAP", "ALDH1L1", "S100B", "AQP4", "SOX9", "CPE", "CLU", "ALDOC", "SLC1A3", "SLC1A2",
  # Oligodendrocytes
  "MOG", "MBP", "CNP", "PLP1", "MAG", "MAL", "CRYAB",
  # Neurons (general)
  "MAP2", "NEUN", "SYN1", "GAP43", "TUBB3",
  # Endothelial Cells
  "PECAM1", "VWF", "CD34", "CDH5", "CLDN5", "ITM2A", "BSG", "RSG5", "APOLD1",
  # OPCs
  "PDGFRA", "NG2", "SOX10", "OLIG2", "NKX2-2", "CNTN1", "TNR",
  # GABAergic Neurons
  "GAD1", "GAD2", "VGAT",
  # Glutamatergic Neurons
  "VGLUT1", "SLC17A7", "CAMK2A", "VGLUT2",
  # Cajal-Retzius Cells
  "RELN", "CALB2", "TP73",
  # Choroid Plexus
  "TTR", "AQP1", "OTX2"
)

# Step 1: Subsample ~50,000 cells
set.seed(123)
sampled_cells <- sample(Cells(seurat_obj), size = 50000)
seurat_obj_sub <- subset(seurat_obj, cells = sampled_cells)
plot_df_sub <- plot_df[rownames(plot_df) %in% sampled_cells, ]

# Step 2: Calculate cluster centers
centers <- plot_df_sub %>%
  group_by(CellType) %>%
  summarise(
    UMAP1 = median(UMAP1),
    UMAP2 = median(UMAP2),
    .groups = "drop"
  )

# Step 3: Plot each marker
for (gene in marker_genes) {
  
  if (gene %in% rownames(seurat_obj_sub)) {
    
    # Create FeaturePlot
    p <- FeaturePlot(seurat_obj_sub, features = gene, reduction = "umap",
                     pt.size = 0.5,
                     min.cutoff = "q10",
                     cols = c("lightgrey", "red")) +
         ggtitle(paste("Expression of", gene)) +
         theme(plot.title = element_text(hjust = 0.5))
    
    # Add cluster titles
    p <- p + 
      geom_text_repel(
        data = centers,
        aes(x = UMAP1, y = UMAP2, label = CellType),
        inherit.aes = FALSE,
        color = "black",
        size = 3,
        fontface = "bold",
        max.overlaps = 50,
        show.legend = FALSE
      )
    
    # Display
    print(p)
    
  } else {
    print(paste("Gene", gene, "not found in the dataset."))
  }
}


###dotplot

library(Seurat)
library(ggplot2)

# Define marker genes grouped by cell types
marker_genes_grouped <- list(
  Microglia = c("CX3CR1", "IBA1", "P2RY12", "TMEM119", "ITGAM", "CSF1R", "HEXB", "SIGLEC11", "SALL1", "TLR2", "TLR4", "CD68", "CD45", "TREM2", "CCL4", "CCL3", "CTSS", "TYROBP", "CD83"),
  Astrocytes = c("GFAP", "ALDH1L1", "S100B", "AQP4", "SOX9", "CPE", "CLU", "ALDOC", "SLC1A3", "SLC1A2"),
  Oligodendrocytes = c("MOG", "MBP", "CNP", "PLP1", "MAG", "MAL", "CRYAB"),
  Neurons = c("MAP2", "NEUN", "SYN1", "GAP43", "TUBB3"),
  Endothelial = c("PECAM1", "VWF", "CD34", "CDH5", "CLDN5", "ITM2A", "BSG", "RSG5", "APOLD1"),
  OPCs = c("PDGFRA", "NG2", "SOX10", "OLIG2", "NKX2-2", "CNTN1", "TNR"),
  GABAergic_Neurons = c("GAD1", "GAD2", "VGAT"),
  Glutamatergic_Neurons = c("VGLUT1", "SLC17A7", "CAMK2A", "VGLUT2"),
  Cajal_Retzius = c("RELN", "CALB2", "TP73"),
  ChoroidPlexus = c("TTR", "AQP1", "OTX2")
)

# Set CellType as active identity
Idents(seurat_obj) <- "CellType"

# Loop through each group and plot separately
for (group in names(marker_genes_grouped)) {
  
  genes <- marker_genes_grouped[[group]]
  genes <- genes[genes %in% rownames(seurat_obj)]  # Filter genes that exist
  
  if (length(genes) > 0) {
    p <- DotPlot(seurat_obj, features = genes, dot.scale = 5, cols = c("lightgrey", "red")) +
      RotatedAxis() +
      ggtitle(paste(group, "Marker Gene Expression"))
    
    print(p)
  }
}


###heatmap

library(Seurat)
library(pheatmap)
library(dplyr)

# Define your marker genes (same expanded list you gave)
marker_genes <- c(
  "CX3CR1", "IBA1", "P2RY12", "TMEM119", "ITGAM", "CSF1R", "HEXB", "SIGLEC11", "SALL1", "TLR2", "TLR4", "CD68", "TREM2", "CCL4", "CCL3", "CTSS", "TYROBP", "CD83",
  "GFAP", "ALDH1L1", "S100B", "AQP4", "SOX9", "CPE", "CLU", "ALDOC", "SLC1A3", "SLC1A2",
  "MOG", "MBP", "CNP", "PLP1", "MAG", "MAL", "CRYAB",
  "MAP2", "SYN1", "GAP43", "TUBB3",
  "PECAM1", "VWF", "CD34", "CDH5", "CLDN5", "ITM2A", "BSG", "APOLD1",
  "PDGFRA", "SOX10", "OLIG2", "NKX2-2", "CNTN1", "TNR",
  "GAD1", "GAD2", "VGAT",
  "SLC17A7", "CAMK2A",
  "RELN", "CALB2", "TP73",
  "TTR", "AQP1", "OTX2"
)

# Make sure genes exist
marker_genes <- marker_genes[marker_genes %in% rownames(seurat_obj)]

# Set identities based on CellType annotation
Idents(seurat_obj) <- "CellType"

# Compute average expression per CellType
avg_expr <- AverageExpression(seurat_obj, features = marker_genes, group.by = "CellType", assays = "RNA", slot = "data")$RNA

# Scale across genes (optional but recommended for heatmap visual clarity)
avg_expr_scaled <- t(scale(t(avg_expr)))

# Plot heatmap
pheatmap(avg_expr_scaled,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         fontsize_row =5,
         fontsize_col = 10,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
         main = "Marker Gene Expression Heatmap")