#!/bin/bash
#PBS -N hvg_selection    
#PBS -l select=1:ncpus=40:mem=480gb
#PBS -l walltime=02:00:00
#PBS -j oe 

# Load or activate your conda environment (adjust to your system)
eval "$(~/anaconda3/bin/conda shell.bash hook)" 
source activate scanvi_update_env

# Make sure we are in the directory from which you submitted the job
cd /rds/general/user/pr422/home/PhD/scanvi_tutorial/human_epilep_scanvi_test/scripts

python <<EOF
import scanpy as sc
import numpy as np
import os

# Path to your concatenated file
input_path = "/rds/general/user/pr422/projects/puklandmarkproject/ephemeral/Parisa/scanvi_output/Combined_symbol_query_concat.h5ad"
# Output file
output_path = "/rds/general/user/pr422/projects/puklandmarkproject/ephemeral/Parisa/scanvi_output/Combined_symbol_query_concat_2000hvg.h5ad"

print("Loading AnnData from:", input_path)
adata = sc.read_h5ad(input_path)
print(f"Original shape: {adata.shape}")
print("Original number of genes:", adata.n_vars)

# -------------------------------------------------------------------
# 1. Highly Variable Gene (HVG) Selection
# -------------------------------------------------------------------
# We select ~2,000 HVGs using Seurat v3 method, accounting for batch
sc.pp.highly_variable_genes(
    adata,
    flavor="seurat_v3",
    n_top_genes=2000,
    batch_key="dataset",
    subset=True
)

print("After HVG selection, shape:", adata.shape)
print("Number of HVGs:", adata.n_vars)

# -------------------------------------------------------------------
# 2. Save the HVG-subset data
# -------------------------------------------------------------------
print(f"Saving HVG-subset data to: {output_path}")
adata.write(output_path)
print("Done.")
EOF