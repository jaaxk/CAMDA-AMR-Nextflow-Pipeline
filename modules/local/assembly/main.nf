process ASSEMBLY {
    label "process_cpu_high"
    conda "bioconda::spades"
    
    publishDir "${params.outdir}/assembly", mode: 'copy'
    input:
    tuple val(accession), val(genus_species), val(antibiotic), val(metadata), path(fastp)
    
    output:
    tuple val(accession), val(genus_species), val(antibiotic), val(metadata), path("${accession}/contigs.fasta"), emit: assembly
    
    script:
    def mode
    def fastp_cmd
    if (fastp instanceof List && fastp.size() == 2) {
        mode = "paired-end"
        def r1 = fastp[0]
        def r2 = fastp[1]
        spades_cmd = """
        spades.py \
            --isolate \
            -1 ${r1} -2 ${r2} \
            -o ${accession} \
            --only-assembler \
            --threads ${task.cpus}
        """
    } else {
        mode = "single-end"
        spades_cmd = """
        spades.py \
            --isolate \
            -s ${fastp} \
            -o ${accession} \
            --only-assembler \
            --threads ${task.cpus}
      """
    }
    """
    echo "Running assembly for ${accession} in ${mode} mode"
    ${spades_cmd}
    """
}
