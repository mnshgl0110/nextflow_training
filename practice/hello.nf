#!/usr/bin/env nextflow
params.greeting = 'Hello world!'
greeting_ch = Channel.of(params.greeting)

process SPLITLETTERS{
    input:
    val x

    output:
    path 'chunk_*'

    script:
    """
    printf '$x' | split -b 6 - chunk_
    """
}

process COVERTTOUPPER {
    input:
    path y

    output:
    stdout

    script:
    """
    cat $y | tr '[a-z]' '[A-Z]'
    #rev $y
    """

}

workflow {
    letters_ch = SPLITLETTERS(greeting_ch)
    result_ch = COVERTTOUPPER(letters_ch.flatten())
    result_ch.view{it}
}