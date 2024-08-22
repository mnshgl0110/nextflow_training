params.reads = "$projectDir/../nf-training/data/ggal/*_{1,2}.fq"
params.transcriptome_file = "$projectDir/../nf-training/data/ggal/transcriptome.fa"
params.multiqc = "$projectDir/multiqc"
params.outdir = "$projectDir/results"

println "reads: $params.reads"



log.info """
    $params.reads
    $params.transcriptome_file
    $params.multiqc
    $params.outdir
"""

process INDEX {
    cpus 2

    input:
    path transcriptome

    output:
    path 'salmon_index'

    script:
    """
    salmon index --threads $task.cpus -t $transcriptome -i salmon_index
    #salmon index -t $transcriptome -i salmon_index
    #salmon --version
    """

}

process QUANTIFICATION {
    tag "Salmon on $sample_id"
    publishDir params.outdir, mode: 'copy'

    input:
    path salmon_index
    tuple val(sample_id), path(reads)



    output:
    path "$sample_id"

    script:
    """
    salmon quant --threads $task.cpus --libType=U -i $salmon_index -1 ${reads[0]} -2 ${reads[1]} -o $sample_id
    """

}

process FASTQC {
    tag "FASTQC on $sample_id"

    input:
    tuple val(sample_id), path(reads)

    output:
    path "fastqc_${sample_id}_logs"

    script:
    """
    mkdir fastqc_${sample_id}_logs
    fastqc -o fastqc_${sample_id}_logs -f fastq -q ${reads}
    """
}

process MULTIQC {
    publishDir params.outdir, mode:'copy'

    input:
    path '*'

    output:
    path 'multiqc_report.html'

    script:
    """
    multiqc .
    """
}


workflow {
    index_ch = INDEX(params.transcriptome_file)
    channel
        .fromFilePairs(params.reads, checkIfExists:true)
       .set{read_pairs_ch}
    quant_ch = QUANTIFICATION(index_ch, read_pairs_ch)
    fastqc_ch = FASTQC(read_pairs_ch)
    multiqc_ch = MULTIQC(quant_ch.mix(fastqc_ch))
}

workflow.onComplete{
    log.info ( workflow.success ? "\nDone! Open the following report in your browser --> $params.outdir/multiqc_report.html\n" : "Oops .. something went wrong" )
}