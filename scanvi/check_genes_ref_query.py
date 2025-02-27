#!/bin/bash
#PBS -N gene_check     
#PBS -l select=1:ncpus=40:mem=480gb
#PBS -l walltime=02:00:00
#PBS -j oe 

# Load or activate your conda environment (adjust to your system)
eval "$(~/anaconda3/bin/conda shell.bash hook)" 
source activate scanvi_update_env 

# Make sure we are in the directory from which you submitted the job
cd /rds/general/user/pr422/home/PhD/scanvi_tutorial/human_epilep_scanvi_test/scripts

# Python inline script:
python <<EOF
import scanpy as sc
import numpy as np
import pandas as pd 
import mygene 

###############################################################################
# 1. Load the Reference Dataset with Ensembl IDs
###############################################################################
ref_path = "/rds/general/user/pr422/home/PhD/ref_data_humanbrain/Combined.h5ad"
query_path = "/rds/general/user/pr422/home/PhD/nfcore/scdownstream_epilep_human/output/run3/combine/merge/merged_inner.h5ad"
save_path = "/rds/general/user/pr422/projects/puklandmarkproject/ephemeral/Parisa/scanvi_output/Combined_symbol.h5ad"

print("Loading reference dataset:", ref_path)
adata_ref = sc.read_h5ad(ref_path)
print(f"Reference shape: {adata_ref.shape}")
print("Sample of reference var_names:", adata_ref.var_names[:10].tolist(), "\\n")

###############################################################################
# 2. Convert Ensembl IDs to Gene Symbols
###############################################################################
# We will:
#  - Strip version numbers (e.g., ENSG00000000419.13 -> ENSG00000000419)
#  - Use 'mygene' to map ENSG IDs to gene symbols
#  - Fallback to original name if no mapping found

mg = mygene.MyGeneInfo()

# Extract unversioned Ensembl IDs where appropriate
ref_var_unversioned = []
for g in adata_ref.var_names:
    if g.startswith("ENSG"):
        ref_var_unversioned.append(g.split(".")[0])
    else:
        # Genes like Cas9 or EGFP remain as-is
        ref_var_unversioned.append(g)

# Query MyGene to map ENSG IDs -> symbols
print("Querying mygene to map Ensembl IDs to gene symbols...")
gene_info = mg.querymany(
    ref_var_unversioned,
    scopes="ensembl.gene",
    fields="symbol",
    species="human",
    as_dataframe=True
)

# Build dictionary: { unversioned_ensg : gene_symbol }
ensg_to_symbol = {}
for idx, row in gene_info.iterrows():
    query_id = idx
    if pd.isnull(row.get("notfound", False)) and isinstance(row.get("symbol"), str):
        ensg_to_symbol[query_id] = row["symbol"]

# Replace var_names in adata_ref with gene symbols if found
new_var_names = []
for original, unversioned in zip(adata_ref.var_names, ref_var_unversioned):
    # If we have a symbol for that Ensembl ID, use it. Otherwise, keep original.
    if unversioned in ensg_to_symbol:
        new_var_names.append(ensg_to_symbol[unversioned])
    else:
        new_var_names.append(original)

adata_ref.var_names = new_var_names
adata_ref.var_names_make_unique()

print("After conversion, sample var_names:", adata_ref.var_names[:10].tolist())

###############################################################################
# 3. Save the Updated Reference Dataset
###############################################################################
print(f"Saving updated reference dataset with gene symbols to: {save_path}")
adata_ref.write(save_path)
print("Saved successfully.\\n")

###############################################################################
# 4. Load the Query Dataset (Gene Symbols) and Check Intersection
###############################################################################
print("Loading query dataset:", query_path)
adata_query = sc.read_h5ad(query_path)
print(f"Query shape: {adata_query.shape}")
print("Sample of query var_names:", adata_query.var_names[:10].tolist(), "\\n")

# Compute intersection
common_genes = np.intersect1d(adata_ref.var_names, adata_query.var_names)
print(f"Number of common genes: {len(common_genes)}")
if len(common_genes) > 0:
    print("Sample of intersecting genes:", common_genes[:20])
EOF