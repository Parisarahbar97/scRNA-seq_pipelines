#!/bin/bash
#PBS -N subset_and_concat    
#PBS -l select=1:ncpus=40:mem=480gb
#PBS -l walltime=02:00:00
#PBS -j oe 


eval "$(~/anaconda3/bin/conda shell.bash hook)" 
source activate scanvi_update_env  

cd /rds/general/user/pr422/home/PhD/scanvi_tutorial/human_epilep_scanvi_test/scripts

python <<EOF
import scanpy as sc
import anndata
import numpy as np
import os

# File paths
ref_path = "/rds/general/user/pr422/projects/puklandmarkproject/ephemeral/Parisa/scanvi_output/Combined_symbol.h5ad"
query_path = "/rds/general/user/pr422/home/PhD/nfcore/scdownstream_epilep_human/output/run3/combine/merge/merged_inner.h5ad"
outdir = "/rds/general/user/pr422/projects/puklandmarkproject/ephemeral/Parisa/scanvi_output"

# 1. Load Updated Reference and Query
print("Loading reference dataset:", ref_path)
adata_ref = sc.read_h5ad(ref_path)
print(f"Reference shape: {adata_ref.shape}")

print("Loading query dataset:", query_path)
adata_query = sc.read_h5ad(query_path)
print(f"Query shape: {adata_query.shape}")

# 2. Find Common Genes
common_genes = np.intersect1d(adata_ref.var_names, adata_query.var_names)
print(f"Number of common genes: {len(common_genes)}")
if len(common_genes) > 0:
    print("Sample of intersecting genes:", common_genes[:20])

# 3. Subset Both Datasets to Common Genes
adata_ref = adata_ref[:, common_genes].copy()
adata_query = adata_query[:, common_genes].copy()

# 4. Save the Subset Data
ref_sub_path = os.path.join(outdir, "Combined_symbol_subset.h5ad")
query_sub_path = os.path.join(outdir, "merged_inner_subset.h5ad")
adata_ref.write(ref_sub_path)
print("Saved subset reference to:", ref_sub_path)

adata_query.write(query_sub_path)
print("Saved subset query to:", query_sub_path)

# 5. Concatenate Datasets
adata_ref.obs['dataset'] = 'reference'
adata_query.obs['dataset'] = 'query'
adata_combined = anndata.concat(
    [adata_ref, adata_query],
    join='inner',
    label='dataset',
    keys=['reference','query'],
    index_unique=None
)

print("Combined shape:", adata_combined.shape)
print(adata_combined)

# 6. Save the Concatenated Result
out_combined = os.path.join(outdir, "Combined_symbol_query_concat.h5ad")
adata_combined.write(out_combined)
print("Saved concatenated file to:", out_combined)
EOF