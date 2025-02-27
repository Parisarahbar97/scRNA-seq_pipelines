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
export SINGULARITY_CACHEDIR="/rds/general/user/$USER/ephemeral/singularity_cache"

# Set a dedicated temporary directory with strict permissions
export TMPDIR="/rds/general/user/$USER/projects/puklandmarkproject/ephemeral/Parisa_scdownstream/tmp" 
mkdir -p "$TMPDIR"
chmod 777 "$TMPDIR" 

# Define inputs and outputs
#sed -i 's/pr422/ah3918/g' //rds/general/user/ah3918/ephemeral/32human_epilep_scdownstream.csv
# samplesheet="/rds/general/user/$USER/projects/puklandmarkproject/ephemeral/Parisa/converted_h5ad_files/32human_epilep_scdownstream.csv"
samplesheet="/rds/general/user/ah3918/ephemeral/32human_epilep_scdownstream.csv"
out="/rds/general/user/$USER/projects/puklandmarkproject/ephemeral/outs/alex_test/"
config_file="/rds/general/user/$USER/ephemeral/PROJECTS/PARISA/scRNA-seq_pipelines/alex_test.config"

# Run the pipeline with resume capability
nextflow run nf-core/scdownstream \
  -r dev \
  -c $config_file \
  --input "$samplesheet" \
  --outdir "$out" \
  -w /rds/general/user/${USER}/ephemeral/NEXTFLOW/ \
  --ambient_removal cellbender \
  --doublet_detection scrublet \
  --skip_liana true \
  --save_intermediates false 