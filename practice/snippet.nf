reads = Channel.fromPath("../nf-training/data/ggal/*_1.fq")
transcriptome = Channel.fromPath("../nf-training/data/ggal/transcriptome.fa")

process PRINT{
    debug true

    input:
    path x
    path y

    script:
    """
    echo $x and $y
    """

}

workflow {
    PRINT(reads, transcriptome.first())
}