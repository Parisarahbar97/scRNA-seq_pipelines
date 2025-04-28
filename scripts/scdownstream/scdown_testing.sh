#for removing the last workdir: 
rm -rf /rds/general/user/pr422/ephemeral/NEXTFLOW/*


export TMPDIR="/rds/general/user/pr422/projects/puklandmarkproject/ephemeral/Parisa_scdownstream/scdown_old_samples/tmp"
mkdir -p "$TMPDIR"
chmod 777 "$TMPDIR"

nextflow run nf-core/scdownstream \
    -profile imperial \
    -c /rds/general/user/pr422/projects/puklandmarkproject/ephemeral/Parisa_scdownstream/scdown_old_samples/test2.config \
    --input /rds/general/user/pr422/projects/puklandmarkproject/ephemeral/Parisa_scdownstream/scdown_old_samples/input/samplesheet_converted.csv \
    --outdir /rds/general/user/pr422/projects/puklandmarkproject/ephemeral/Parisa_scdownstream/scdown_old_samples/output \
    -w /rds/general/user/pr422/ephemeral/NEXTFLOW/ \
    --scanvi_model /rds/general/user/pr422/projects/puklandmarkproject/ephemeral/Parisa/train_ref_scdown/scanvi_model/model.pt \
	  --base_adata /rds/general/user/pr422/projects/puklandmarkproject/ephemeral/Parisa/train_ref_scdown/scanvi_model/adata.h5ad \
	  --base_label_col cell_type \
	  --integration_methods scanvi \
    --ambient_removal cellbender \
    --doublet_detection scrublet \
    --clustering_resolutions "0.5,1.0" \
    --integration_hvgs 10000 \
    --cluster_per_label \
    --cluster_global \
    --memory_scale 1 \
    --skip_liana \
    --save_intermediates false \
    > /rds/general/user/pr422/projects/puklandmarkproject/ephemeral/Parisa_scdownstream/scdown_old_samples/output/test.log 2>&1

