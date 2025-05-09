NXF_LOG_FILE="/rds/general/user/pr422/ephemeral/NEXTFLOW_EQTL/nextflow.log"


nextflow run HaglundA/eQTL_PIPELINE -r dev \
-profile imperial \
--outdir /rds/general/user/pr422/ephemeral/NEXTFLOW_EQTL/eQTL_output \
--gds_file /rds/general/user/pr422/projects/puklandmarkproject/ephemeral/alex/sc_eQTL_runs/BONN_post_imputation_QC.gds \
--single_cell_file /rds/general/user/pr422/projects/puklandmarkproject/live/Users/Parisa/parisa_eqtl/updated_seurat_with_individual_id.rds \
--celltype_column "cell_type" \
--individual_column "Individual_ID" \
--counts_assay "originalexp" \
--counts_slot "counts" \
--optimize_pcs true \
-with-report pipeline_report.html \
--report true