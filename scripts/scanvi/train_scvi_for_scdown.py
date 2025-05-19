#!/usr/bin/env python

import scanpy as sc
import scvi
import numpy as np

# ---- 1. Load the reference AnnData ----
print("Loading AnnData")
adata = sc.read_h5ad("/rds/general/user/pr422/home/PhD/ref_data_humanbrain/Combined.h5ad")

# ---- 2. Check AnnData’s info ----
print("AnnData info:")
print("  shape:", adata.shape)
print("  obs columns:", list(adata.obs.columns))
print("  var columns:", list(adata.var.columns))
print("  Example var head:\n", adata.var.head())

# ---- 3. Set gene symbols from 'Gene' column as feature names ---
if "Gene" in adata.var.columns:
    print("Setting .var_names to 'Gene'")
    adata.var_names = adata.var["Gene"]
    adata.var_names = adata.var_names.astype(str)
    adata.var_names_make_unique()
else:
    raise ValueError("'Gene' column not found in .var — cannot set gene symbols.")

# ---- 4. Filter out 'Cas9' and 'EGFP' ----
remove_names = {"Cas9", "EGFP", "", None}
mask = ~adata.var_names.isin(remove_names)
before = adata.n_vars
adata = adata[:, mask].copy()
after = adata.n_vars
print(f"Filtered out {before - after} non-endogenous or unwanted genes. Remaining genes: {after}")

# ---- 5. Remove duplicate gene symbols ----
print("Checking for duplicate gene symbols")
dups = adata.var_names.duplicated()
if np.any(dups):
    print(f"Found {dups.sum()} duplicates. Removing them.")
    adata = adata[:, ~dups].copy()
else:
    print("No duplicates found.")

# ---- 6. Subsample if necessary ----
print("Subsampling to 10,000 cells if needed")
if adata.n_obs > 10000:
    np.random.seed(0)
    idx = np.random.choice(adata.n_obs, 10000, replace=False)
    adata = adata[idx, :].copy()
    print("Subsampled shape:", adata.shape)
else:
    print(f"Only {adata.n_obs} cells — skipping subsampling.")

# ---- 7. Check cell type column ----
print("Verifying 'cell_type' column...")
assert "cell_type" in adata.obs.columns, " 'cell_type' column not found in adata.obs!"

# ---- 8. Setup for scVI ----
print("Setting up AnnData for SCVI")
scvi.model.SCVI.setup_anndata(adata, batch_key=None)

# ---- 9. Train SCVI ----
print("Training SCVI")
model = scvi.model.SCVI(adata)
model.train()

# ---- 10. Transfer to SCANVI ----
print("Transferring to SCANVI")
scanvi_model = scvi.model.SCANVI.from_scvi_model(
    model,
    unlabeled_category="unlabeled_dummy",
    labels_key="cell_type"
)
scanvi_model.train()

# ---- 11. Save model ----
import os
output_dir = "/rds/general/user/pr422/projects/puklandmarkproject/live/Users/Parisa/parisa_scdownstream/scanvi_output_for_scdown/model"

# Check if the directory already exists
if os.path.exists(output_dir):
    raise FileExistsError(f"Output directory already exists: {output_dir}. Please choose a new one.")

print(f"Saving model and AnnData to {output_dir}")
scanvi_model.save(output_dir, save_anndata=True)

print("Model training and saving complete.")