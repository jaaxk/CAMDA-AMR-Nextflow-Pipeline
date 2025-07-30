process FLANK {
    conda "bin/envs/flank_env.yml"
    publishDir "${params.outdir}/flank", mode: 'copy'
    
    input:
    tuple val(accession), val(genus_species), val(antibiotic), val(metadata), path(assembly), path(blast)
    path script_file
    
    output:
    tuple val(accession), val(genus_species), val(antibiotic), val(metadata), path("${accession}_flanks.fasta"), emit: flanks
    
    script:
    """
    echo "Getting flanking regions for ${accession}"   
    python ${script_file} ${assembly} ${params.input_seq_length} ${blast} ${accession}_flanks.fasta
    """
}