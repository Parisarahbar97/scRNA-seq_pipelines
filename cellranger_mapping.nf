#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Parameters
params.samples_csv = "/rds/general/user/pr422/home/samples_part2.csv" 
params.outdir = "/rds/general/user/pr422/home/PhD/DATA/epilep/mouse/hippo/output"
params.genome = "/rds/general/user/pr422/home/PhD/refrence_genome/mouse/refdata-gex-GRCm39-2024-A"
params.cellranger = "/rds/general/user/pr422/home/PhD/cellranger/cellranger-8.0.1/cellranger"

// Cell Ranger process
process runCellRanger {
    publishDir "${params.outdir}", mode: 'copy'

    conda "/rds/general/user/pr422/home/anaconda3/envs/cellranger_env"

    executor = "pbspro"
    clusterOptions = '-lselect=1:ncpus=40:mem=480gb -l walltime=08:00:00'

    input:
    tuple val(sample), val(fastq_dirs)

    output:
    path("${sample}_outs")

    script:
    """
    echo "Processing sample: ${sample}"
    echo "FASTQ directories: ${fastq_dirs}"

    ${params.cellranger} count \
        --id=${sample}_outs \
        --transcriptome=${params.genome} \
        --fastqs=${fastq_dirs.join(',')} \
        --sample=${sample} \
        --create-bam true
    """
}

// Workflow definition
workflow {

    println "Starting the mapping pipeline"
    samples = Channel
        .fromPath(params.samples_csv)
        .splitCsv(header: true)
        .map { row ->
            def sample = row['sample']
            // Remove extra quotes from the 'fastq_dirs' field
            def fastq_dirs_str = row['fastq_dirs'].replaceAll(/^"|"$/, '')
            def fastq_dirs = fastq_dirs_str.split(',').collect { it.replaceAll(/^"|"$/, '') }
            tuple(sample, fastq_dirs)
        }

    samples.view()

    samples | runCellRanger
}
