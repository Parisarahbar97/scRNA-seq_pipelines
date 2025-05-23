export NXF_LOG_FILE="/rds/general/user/pr422/ephemeral/NEXTFLOW_EQTL/nextflow.log"

# Set path to remove old cached pipeline repo
REPO_PATH="/rds/general/user/pr422/home/.nextflow/assets/johnsonlab-ic/sc-eQTL-pipeline"
if [ -d "$REPO_PATH" ]; then
    echo "Removing cached repository at $REPO_PATH"
    rm -rf "$REPO_PATH"
fi


outdir="/rds/general/user/pr422/projects/puklandmarkproject/live/Users/Parisa/parisa_eqtl/eQTL_output"
cov_file="/rds/general/user/pr422/projects/puklandmarkproject/live/Users/Parisa/parisa_eqtl/epilepsy_metadata_matrixeqtl.csv"
gds_file="/rds/general/user/pr422/projects/puklandmarkproject/ephemeral/alex/sc_eQTL_runs/BONN_post_imputation_QC.gds"
single_cell_file="/rds/general/user/pr422/projects/puklandmarkproject/live/Users/Parisa/parisa_eqtl/matched_seuratobj/matched_seurat_demuxlet.rds"


nextflow run johnsonlab-ic/sc-eQTL-pipeline \
-profile imperial \
--outdir $outdir \
--gds_file $gds_file \
--single_cell_file $single_cell_file \
--cov_file $cov_file \
--celltype_column "cell_type" \
--individual_column "Individual_ID" \
--counts_assay "originalexp" \
--counts_slot "counts" \
--with-report pipeline_report.html \
--report true \
--covariates_to_include "Sex,Pathology,AgeAtSurgery"