import scvi
import scanpy as sc
import numpy as np

print("Loading AnnData")
adata = sc.read_h5ad("/rds/general/user/pr422/home/PhD/ref_data_humanbrain/Combined.h5ad")

print("Subsampling to 10,000 cells for quick testing")
adata = adata[np.random.RandomState(seed=0).choice(adata.shape[0], 10000, replace=False)].copy()

print("Checking for 'cell_type' in adata.obs")
assert "cell_type" in adata.obs.columns

print("Setting up AnnData for scVI")
scvi.model.SCVI.setup_anndata(adata, batch_key=None)

print("Training scVI...")
model = scvi.model.SCVI(adata)
model.train()

print("Transferring to scANVI...")
scanvi_model = scvi.model.SCANVI.from_scvi_model(
    model,
    unlabeled_category="unlabeled_dummy",  # Dummy category (not in real labels)
    labels_key="cell_type"
)
scanvi_model.train()

print("Saving model...")
scanvi_model.save(
    "/rds/general/Users/Parisa/ephemeral_files/Parisa_2/Parisa/train_ref_scdown/scanvi_model/",
    save_anndata=True
)

print("✅ Model training complete.")