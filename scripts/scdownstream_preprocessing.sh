#!/bin/bash
#PBS -l walltime=12:00:00
#PBS -l select=1:ncpus=40:mem=480GB
#PBS -N scdownstream_analysis
#PBS -o scdownstream_analysis.out
#PBS -e scdownstream_analysis.err

# Load necessary modules and set environment variables
module load Java/11.0.20

export NXF_OPTS='-Xms1g -Xmx4g'
export JAVA_HOME=$(dirname $(dirname $(which java)))
export PATH=$JAVA_HOME/bin:$PATH
export SINGULARITY_CACHEDIR="/rds/general/user/pr422/ephemeral/singularity_cache"

# Set a dedicated temporary directory with strict permissions
export TMPDIR="/rds/general/user/pr422/projects/puklandmarkproject/ephemeral/Parisa_scdownstream/tmp" 
mkdir -p "$TMPDIR"
chmod 777 "$TMPDIR" 

# Define inputs and outputs
samplesheet="/rds/general/user/pr422/projects/puklandmarkproject/ephemeral/Parisa/converted_h5ad_files/32human_epilep_scdownstream.csv"
out="/rds/general/user/pr422/projects/puklandmarkproject/ephemeral/Parisa_scdownstream/output_preprocess/run1"
config_file="/rds/general/user/pr422/projects/puklandmarkproject/ephemeral/Parisa_scdownstream/template/test.config"

# Run the pipeline with resume capability
nextflow run nf-core/scdownstream \
  -profile imperial \
  -c "$config_file" \
  --input "$samplesheet" \
  --outdir "$out" \
  -w /rds/general/user/${USER}/ephemeral/NEXTFLOW/ \
  --ambient_removal cellbender \
  --doublet_detection scrublet \
  --memory_scale 1 \
  --skip_liana true \
  --save_intermediates false 