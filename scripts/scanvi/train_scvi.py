#!/bin/bash
#PBS -N scvi_train    
#PBS -l select=1:ncpus=40:mem=480gb
#PBS -l walltime=12:00:00
#PBS -j oe 

# Load or activate your conda environment (adjust to your system)
eval "$(~/anaconda3/bin/conda shell.bash hook)"
source activate scanvi_update_env  

# Make sure we are in the directory from which you submitted the job
cd /rds/general/user/pr422/home/PhD/scanvi_tutorial/human_epilep_scanvi_test/scripts

python <<EOF
import scanpy as sc
import scvi
import os

# -------------------------------------------------------------------
# 1. PATHS
# -------------------------------------------------------------------
input_file = "/rds/general/user/pr422/projects/puklandmarkproject/ephemeral/Parisa/scanvi_output/Combined_symbol_query_concat_2000hvg.h5ad"
out_dir = "/rds/general/user/pr422/projects/puklandmarkproject/ephemeral/Parisa/scanvi_output"
os.makedirs(out_dir, exist_ok=True)

# final integrated AnnData path:
output_h5ad = os.path.join(out_dir, "scvi_integrated_2000hvg.h5ad")

# UMAP plot path:
sc.settings.figdir = out_dir
umap_plot_file_suffix = "_scvi_umap_integration_2.png"

# -------------------------------------------------------------------
# 2. LOAD THE HVG-SUBSET DATA
# -------------------------------------------------------------------
print(f"Loading HVG-subset AnnData from: {input_file}")
adata = sc.read_h5ad(input_file)
print(f"Data shape: {adata.shape}")
print(f"Observations: {adata.n_obs}, Genes: {adata.n_vars}")

# -------------------------------------------------------------------
# 3. SETUP ANNDATA FOR SCVI
# -------------------------------------------------------------------
# We assume the batch key is 'dataset' (reference vs query)
print("Setting up anndata for scVI with batch_key='dataset'...")
scvi.model.SCVI.setup_anndata(adata, batch_key="dataset")

# -------------------------------------------------------------------
# 4. TRAIN SCVI
# -------------------------------------------------------------------
print("Initializing SCVI model (n_latent=30, n_layers=2)...")
scvi_model = scvi.model.SCVI(adata, n_latent=30, n_layers=2)

print("Training scVI model for up to 50 epochs...")
scvi_model.train(max_epochs=50)

# Get the scVI latent representation
adata.obsm["X_scVI"] = scvi_model.get_latent_representation()

# -------------------------------------------------------------------
# 5. COMPUTE NEIGHBORS AND UMAP FOR VISUALIZATION
# -------------------------------------------------------------------
print("Computing neighbors on 'X_scVI'...")
sc.pp.neighbors(adata, use_rep="X_scVI")

print("Computing UMAP (min_dist=0.3)...")
sc.tl.umap(adata, min_dist=0.3)

# -------------------------------------------------------------------
# 6. SAVE UMAP PLOT
# -------------------------------------------------------------------
print(f"Saving UMAP plot to directory: {out_dir}")
sc.pl.umap(
    adata,
    color=["dataset"],  
    save=umap_plot_file_suffix,
    show=False
)

# -------------------------------------------------------------------
# 7. SAVE THE UPDATED AnnData
# -------------------------------------------------------------------
print(f"Saving final integrated AnnData to: {output_h5ad}")
adata.write(output_h5ad)
print("Done.")
EOF