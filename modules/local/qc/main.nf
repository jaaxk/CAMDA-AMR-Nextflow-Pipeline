process QC {
    label "process_cpu_high"
    conda "bioconda::fastp"
    
    publishDir "${params.outdir}/fastp", mode: 'copy'
    input: 
    tuple val(accession), val(genus_species), val(antibiotic), val(metadata), path(fastq)

    output:
    tuple val(accession), val(genus_species), val(antibiotic), val(metadata), path("${accession}*.trimmed.fastq.gz"), emit: fastp
    path("${accession}_fastp.html"), emit: report_html
    path("${accession}_fastp.json"), emit: report_json

    script:
    
    def mode
    def fastp_cmd

    if (fastq instanceof List && fastq.size() == 2) {
        mode = "paired-end"
        def r1 = fastq[0]
        def r2 = fastq[1]
        fastp_cmd = """
        fastp \
            --in1 ${r1} --in2 ${r2} \
            --out1 ${accession}_1.trimmed.fastq.gz \
            --out2 ${accession}_2.trimmed.fastq.gz \
            --thread ${task.cpus} \
            --html ${accession}_fastp.html \
            --json ${accession}_fastp.json
        """
    } else {
        mode = "single-end"
        fastp_cmd = """
        fastp \
            --in1 ${fastq} \
            --out1 ${accession}.trimmed.fastq.gz \
            --thread ${task.cpus} \
            --html ${accession}_fastp.html \
            --json ${accession}_fastp.json
        """
    }

    """
    echo "Running QC for ${accession} in ${mode} mode"
    ${fastp_cmd}
    """
}
