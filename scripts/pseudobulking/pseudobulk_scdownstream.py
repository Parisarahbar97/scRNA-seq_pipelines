#!/usr/bin/env python3

import os
import scanpy as sc
import pandas as pd
from scipy import sparse

def pseudo_bulk(adata: sc.AnnData,
                celltype_column: str,
                sample_column: str) -> dict:
    """
    Pseudobulk a Scanpy AnnData from scdownstream pipeline by summing counts for each gene,
    per sample, within each cell type.
    """
    # 1) pull the count matrix: raw if available, otherwise X
    if adata.raw is not None and adata.raw.X is not None:
        X = adata.raw.X
        genes = adata.raw.var_names
    else:
        X = adata.X
        genes = adata.var_names

    # 2) To dense
    if sparse.issparse(X):
        X = X.toarray()

    # 3) Build cells×genes DataFrame
    df = pd.DataFrame(X, index=adata.obs_names, columns=genes)

    # 4) Group by (cell type, sample) and sum
    grouped = df.groupby(
        [adata.obs[celltype_column], adata.obs[sample_column]]
    ).sum()

    # 5) Split into one genes×samples DataFrame per cell type
    out = {}
    for ct in adata.obs[celltype_column].unique():
        if ct in grouped.index.get_level_values(0):
            # samples × genes → transpose → genes × samples
            out[ct] = grouped.loc[ct].T.astype(int)
        else:
            out[ct] = pd.DataFrame(index=genes, columns=[])
    return out

if __name__ == '__main__':
    # --- Config ---
    h5ad_path    = '/mnt/epilep/merged.h5ad'
    celltype_col = 'celltypist:Human_AdultAged_Hippocampus'
    sample_col   = 'sample'
    out_dir      = '/mnt/epilep/pseudobulk_outputs'

    # --- Run ---
    adata = sc.read_h5ad(h5ad_path)
    pseudo_dict = pseudo_bulk(adata, celltype_col, sample_col)

    os.makedirs(out_dir, exist_ok=True)
    for ct, df in pseudo_dict.items():
        # move gene names into a column, reorder so geneid is last
        df_out = (
            df.reset_index()
              .rename(columns={'index': 'geneid'})
        )
        samples = [c for c in df_out.columns if c != 'geneid']
        df_out = df_out[samples + ['geneid']]

        fname = os.path.join(out_dir, f"{ct.replace('/', '_')}_pseudobulk.csv")
        df_out.to_csv(fname, index=False)
        print(f"Wrote {fname}: {df_out.shape[0]} genes × {df_out.shape[1]-1} samples")