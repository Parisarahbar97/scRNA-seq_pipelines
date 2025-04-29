####BEFORE PLOTTING UMAPS:

library(SingleCellExperiment)
library(tibble)   
library(dplyr)   
library(ggplot2)
library(Seurat)

tmp=readRDS("merged.rds")
umap_coords <- reducedDim(tmp, "X_scvi-global_umap")
plot_df <- data.frame(
  UMAP1 = umap_coords[, 1],
  UMAP2 = umap_coords[, 2],
  CellType = tmp$celltypist.Human_AdultAged_Hippocampus,
  Sample = tmp$sample
)

####PLOTTING:
####LEGEND ON THE RIGHT AND TITLED CLUSTERS
library(ggplot2)
library(ggsci)
library(dplyr)
library(ggrepel)

celltype_colors <- c(
  "Microglia" = "#E41A1C",
  "Astrocytes" = "#377EB8",
  "OPCs" = "#4DAF4A",
  "mOli" = "#984EA3",
  "mGC" = "#FF7F00",
  "CA1_neurons" = "#A65628",
  "CA2-4_neurons" = "#F781BF",
  "CA_neurons" = "#999999",
  "GABA_neurons" = "#66C2A5",
  "Endothelial" = "#E6AB02",
  "Ependymal" = "#A6CEE3",
  "Cajal-Retzius" = "#B2DF8A",
  "Subcul_EntorCtx_neurons" = "#FB9A99",
  "imGC" = "#FDBF6F",
  "ChoroidPlexus" = "#CAB2D6"
)

# Build and assign the full plot to p
p <- ggplot(plot_df, aes(x = UMAP1, y = UMAP2, color = CellType)) +
  geom_point(size = 0.01, alpha = 0.01) +
  geom_text_repel(
    data = centers,
    aes(label = CellType),
    color = "black",
    size = 3,
    fontface = "bold",
    max.overlaps = 50,
    show.legend = FALSE
  ) +
  theme_classic() +
  scale_color_manual(values = celltype_colors) + 
  labs(title = "UMAP of Brain CellTypes", color = "Cell Type") +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  ) +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 3)))

# Explicitly print it in Jupyter (so it’s the exact same object)
print(p)

# Save the same exact object
ggsave(
  filename = "/data/epilep/figures/UMAP_Brain_CellTypist_AdultAged_Hippocampus_2.png",
  plot = p
)

####ONLY LEGEND ON THE RIGHT

# Plot with legend only — no text labels
p <- ggplot(plot_df, aes(x = UMAP1, y = UMAP2, color = CellType)) +
  geom_point(size = 0.01, alpha = 0.01) +
  theme_classic() +
  scale_color_manual(values = celltype_colors) +
  labs(title = "UMAP of Brain CellTypes", color = "Cell Type") +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  ) +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 3)))

# Display the plot in notebook
print(p)

# Save to file
ggsave(
  filename = "/data/epilep/figures/UMAP_Brain_CellTypist_AdultAged_Hippocampus_untitled.png",
  plot = p
)

####ONLY TITLED CLUSTERS

# 2. Build plot with cluster titles, no legend
p <- ggplot(plot_df, aes(x = UMAP1, y = UMAP2, color = CellType)) +
  geom_point(size = 0.01, alpha = 0.01) +
  geom_text_repel(
    data = centers,
    aes(label = CellType),
    color = "black",
    size = 3,
    fontface = "bold",
    max.overlaps = 50,
    show.legend = FALSE
  ) +
  theme_classic() +
  scale_color_manual(values = celltype_colors) +
  labs(title = "UMAP of Brain CellTypes") +
  theme(
    legend.position = "none",                # ❌ No legend
    plot.title = element_text(hjust = 0.5)   # ✅ Centered title
  )

# 3. Display in notebook
print(p)

# 4. Save as PNG
ggsave(
  filename = "/data/epilep/figures/UMAP_Brain_CellTypist_AdultAged_Hippocampus_titled.png",
  plot = p
)

####SUBSAMPLING AND RE-PLOTTING (TITLED CLUSTERS)

library(dplyr)
library(ggplot2)
library(ggrepel)

# 1. Subsample 50,000 cells randomly
set.seed(123)  # for reproducibility
plot_df_sub <- plot_df %>%
  sample_n(50000)

# 2. Recalculate cluster centers (optional: better if based on full plot_df or keep same)
centers_sub <- plot_df_sub %>%
  group_by(CellType) %>%
  summarise(
    UMAP1 = median(UMAP1),
    UMAP2 = median(UMAP2),
    .groups = "drop"
  )

# 3. Plot the subsampled data
p_sub <- ggplot(plot_df_sub, aes(x = UMAP1, y = UMAP2, color = CellType)) +
  geom_point(size = 0.01, alpha = 0.5) +   # Bigger points and less transparency now!
  geom_text_repel(
    data = centers_sub,
    aes(label = CellType),
    color = "black",
    size = 4,
    fontface = "bold",
    max.overlaps = 50,
    show.legend = FALSE
  ) +
  theme_classic() +
  scale_color_manual(values = celltype_colors) +
  labs(title = "Subsampled UMAP (50,000 cells)") +
  theme(
    legend.position = "none",   # No legend
    plot.title = element_text(hjust = 0.5)
  )

# 4. Show in notebook
print(p_sub)

ggsave(
  filename = "/data/epilep/figures/UMAP_CellTypes_TitlesOnly_Subsample50K_2.png",
  plot = p_sub
)