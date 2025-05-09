
main_dir=/rds/general/user/ah3918/projects/puklandmarkproject/
export NXF_LOG_FILE=/rds/general/user/${USER}/ephemeral/NEXTFLOW/scdownstream.log
mkdir -p /rds/general/user/ah3918/ephemeral/tmp

# Directory containing the mapping outputs
MAPPING_DIR="/rds/general/user/ah3918/projects/puklandmarkproject/ephemeral/alex/mapping_outs"
# Output CSV file path
OUTPUT_CSV="$MAPPING_DIR/samplesheet_h5_files.csv"

# Write header
echo "sample,filtered,unfiltered" > "$OUTPUT_CSV"

# Find all directories ending with "_mapped" and extract sample names and paths
for dir in "$MAPPING_DIR"/*_mapped; do
  if [ -d "$dir" ]; then
    # Extract sample name (remove "_mapped" suffix)
    sample_name=$(basename "$dir" | sed 's/_mapped$//')
    
    # Check if required files exist
    filtered_path="$dir/outs/filtered_feature_bc_matrix.h5"
    raw_path="$dir/outs/raw_feature_bc_matrix.h5"
    
    if [[ -f "$filtered_path" && -f "$raw_path" ]]; then
      # Write to CSV
      echo "$sample_name,$filtered_path,$raw_path" >> "$OUTPUT_CSV"
    else
      echo "Warning: Missing h5 files for sample $sample_name" >&2
    fi
  fi
done



# Define inputs and outputs
out="/rds/general/user/ah3918/ephemeral/scdownstream_outs"
config_file="/rds/general/user/ah3918/home/scRNA-seq_pipelines/scripts/scdownstream/alex.config"

# Run the pipeline with resume capability
nextflow run nf-core/scdownstream \
  -r dev \
  -profile imperial \
  -c "$config_file" \
  --input "$OUTPUT_CSV" \
  --outdir "$out" \
  --email ah3918@ic.ac.uk \
  -w /rds/general/user/${USER}/ephemeral/NEXTFLOW/ \
  --ambient_removal cellbender \
  --doublet_detection scrublet \
  --clustering_resolutions "0.5,1.0" \
  --integration_methods "scvi" \
  --integration_hvgs 10000 \
  --cluster_per_label \
  --cluster_global \
  --celltypist_model "Adult_Human_PrefrontalCortex" \
  --memory_scale 1 \
  --skip_liana true \
  --save_intermediates true