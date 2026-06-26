#!/usr/bin/env nextflow


process fetch_reference {

    publishDir "${params.out}/reference"

    conda "bioconda::entrez-direct=24.0"

    input:
        val accession

    output:
        path "${accession}.fasta"

    script:
    """
    esearch -db nucleotide -query "$accession" \
    | efetch -format fasta > "${accession}.fasta"
    """
}

process combine_fastas {

    publishDir "${params.out}/combined"

    input:
        path genomes
        path reference

    output:
        path "combined.fasta"

    script:
    """
    cat "$reference" ${genomes} > combined.fasta
    """
}

process mafft_align {

    publishDir "${params.out}/alignment"

    conda "bioconda::mafft=7.525"

    input:
        path fasta

    output:
        path "alignment.fasta"

    script:
    def args = task.ext.args ?: ""

    """
    mafft $args "$fasta" > alignment.fasta
    """
}

process trim_alignment {

    publishDir "${params.out}/results"

    conda "bioconda::trimal=1.5.0"

    input:
        path alignment

    output:
        path "alignment_trimmed.fasta"
        path "alignment_trimmed.html"

    script:
    def args = task.ext.args ?: ""

    """
    trimal \
        $args \
        -in "$alignment" \
        -out alignment_trimmed.fasta \
        -htmlout alignment_trimmed.html
    """
}

workflow {

    if (!params.in) {
        error "Please provide an input directory using --in parameter."
    }

    def ch_ref = fetch_reference(params.accession)

    def ch_genomes = channel.fromPath("${params.in}/*.fasta")

    def ch_combined = combine_fastas(ch_genomes.collect(), ch_ref)

    def ch_alignment = mafft_align(ch_combined)

    trim_alignment(ch_alignment)
}