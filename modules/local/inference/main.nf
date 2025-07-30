process INFERENCE {
    label "process_a100"
    conda 'bin/envs/infer_env.yml'
    publishDir "${params.outdir}/inference", mode: 'copy'
    
    input:
    tuple val(antibiotic), val(samples)
    path script_file

    output:
    tuple path("${antibiotic}_results.tsv"), emit: inference

    script:
    def models_path = file("bin/resources/models/2class_5e-2")
    """
    echo "Running inference for ${antibiotic}"
    echo -e "${samples.collect{ it.join('\t') }.join('\n')}" > samples.tsv
    python ${script_file} ${antibiotic} samples.tsv ${models_path} ${params.model_max_length} ${params.num_hits_norm}
    """

}