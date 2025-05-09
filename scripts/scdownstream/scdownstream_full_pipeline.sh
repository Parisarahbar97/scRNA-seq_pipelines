# Set a dedicated temporary directory with strict permissions
export TMPDIR="/rds/general/user/pr422/projects/puklandmarkproject/ephemeral/Parisa_scdownstream/tmp" 
mkdir -p "$TMPDIR"
chmod 777 "$TMPDIR" 

# Define inputs and outputs
samplesheet="/rds/general/user/pr422/projects/puklandmarkproject/ephemeral/Parisa/converted_h5ad_files/32human_epilep_scdownstream.csv"
out="/rds/general/user/pr422/projects/puklandmarkproject/ephemeral/Parisa_scdownstream/output_full/run3_celltypist_new_model"
config_file="/rds/general/user/pr422/projects/puklandmarkproject/ephemeral/Parisa_scdownstream/scRNA-seq_pipelines/scripts/scdownstream/test_full.config"

# Run the pipeline with resume capability
nextflow run nf-core/scdownstream \
  -profile imperial > /rds/general/user/pr422/projects/puklandmarkproject/ephemeral/Parisa_scdownstream/output_full/run3_celltypist_new_model/test.log \
  -c "$config_file" \
  --input "$samplesheet" \
  --outdir "$out" \
  -w /rds/general/user/${USER}/ephemeral/NEXTFLOW/ \
  --ambient_removal cellbender \
  --doublet_detection scrublet \
  --clustering_resolutions "0.5,1.0" \
  --integration_methods "scvi" \
  --integration_hvgs 10000 \
  --cluster_per_label \
  --cluster_global \
  --celltypist_model "Human_AdultAged_Hippocampus" \
  --memory_scale 1 \
  --skip_liana true \
  --save_intermediates false