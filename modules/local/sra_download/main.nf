process SRA_DOWNLOAD {
    conda "bioconda::sra-tools"
    label 'process_medium'
    
    publishDir "${params.outdir}/fastq", mode: 'copy'
    
    input:
    tuple val(accession), val(genus_species), val(antibiotic), val(metadata)
    
    output:
    tuple val(accession), val(genus_species), val(antibiotic), val(metadata), path("${accession}*.fastq.gz"), emit: fastq
    
    script:
    // Process only uses accession for download, but passes everything through
    """
    fasterq-dump "${accession}" --split-files --threads ${task.cpus} --outdir .

    gzip *.fastq
    
    ls
    """
}