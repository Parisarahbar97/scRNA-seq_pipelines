# Set main paths
MAPPING_DIR="/rds/general/user/pr422/projects/puklandmarkproject/ephemeral/alex/mapping_outs"
OUTPUT_CSV="/rds/general/user/pr422/projects/puklandmarkproject/live/Users/Parisa/parisa_scdownstream/input/samplesheet_with_filtering_parameters.csv"

# Create directory if it doesn't exist
mkdir -p $(dirname "$OUTPUT_CSV")

# Write CSV header with extra columns
echo "sample,filtered,unfiltered,min_genes,min_counts_cell,max_mito_percentage" > "$OUTPUT_CSV"

# Loop through sample directories
for dir in "$MAPPING_DIR"/*_mapped; do
  if [ -d "$dir" ]; then
    sample_name=$(basename "$dir" | sed 's/_mapped$//')
    filtered_path="$dir/outs/filtered_feature_bc_matrix.h5"
    raw_path="$dir/outs/raw_feature_bc_matrix.h5"

    if [[ -f "$filtered_path" && -f "$raw_path" ]]; then
      # Add sample with default QC values
      echo "$sample_name,$filtered_path,$raw_path,300,500,5" >> "$OUTPUT_CSV"
    else
      echo "Warning: Missing h5 files for sample $sample_name" >&2
    fi
  fi
done

#Verify the CSV file:
head -n 5 "$OUTPUT_CSV"


export NXF_LOG_FILE=/rds/general/user/pr422/ephemeral/NEXTFLOW/scdownstream.log

nextflow run nf-core/scdownstream \
  -r 85a3ca5c56 \
  -profile imperial \
  -c "/rds/general/user/pr422/projects/puklandmarkproject/ephemeral/Parisa_scdownstream/scRNA-seq_pipelines/scripts/scdownstream/final.config" \
  --input "/rds/general/user/pr422/projects/puklandmarkproject/live/Users/Parisa/parisa_scdownstream/input/samplesheet_with_filtering_parameters.csv" \
  --outdir "/rds/general/user/pr422/projects/puklandmarkproject/live/Users/Parisa/parisa_scdownstream/output/with_filtering_parameters" \
  --email pr422@ic.ac.uk \
  -w /rds/general/user/pr422/ephemeral/NEXTFLOW/ \
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
  --save_intermediates true 