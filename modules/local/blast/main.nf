process BLAST {
    conda "bioconda::blast"
    
    publishDir "${params.outdir}/blast", mode: 'copy'
    
    input:
    tuple val(accession), val(genus_species), val(antibiotic), val(metadata), path(assembly)
    
    output:
    tuple val(accession), val(genus_species), val(antibiotic), val(metadata), path(assembly), path("${accession}_hits.tsv"), emit: blast
    
    script:
    def sig_seqs = file("bin/resources/sig_seqs_5e-2/${genus_species}_sig_sequences.fasta")

    """
    echo "Running BLAST for ${accession}"
    makeblastdb -in ${assembly} -dbtype nucl -out ${accession}
    blastn -query ${sig_seqs} -db ${accession} \
    -max_target_seqs ${params.blast_max_target_seqs} \
    -perc_identity ${params.blast_percid} \
    -outfmt "6 qseqid sseqid pident length qstart qend sstart send sstrand bitscore" \
    -out ${accession}_hits.tsv 

    """
}