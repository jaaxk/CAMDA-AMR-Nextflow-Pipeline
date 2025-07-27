#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Parameters
params.input = "data/samples.tsv"
params.outdir = "results"

// Include modules
include { SRA_DOWNLOAD } from './modules/local/sra_download/main'
include { QC } from './modules/local/QC/main'
include { ASSEMBLY } from './modules/local/assembly/main'
include { BLAST } from './modules/local/blast/main'
include { FLANK } from './modules/local/flank/main'
include { INFERENCE } from './modules/local/inference/main'

workflow {
    // Create a channel with ALL metadata for each sample
    samples_ch = Channel
        .fromPath(params.input, checkIfExists: true)
        .splitCsv(header: true, sep: '\t')
        .map { row -> 
            // Keep ALL the data together
            [row.accession, row.genus_species, row.antibiotic, row]
        }

   /*  fastq_ch = Channel.of(
        tuple(
            "SRR5386858",                        // accession
            "klebsiella_pneumoniae",             // genus_species
            "GEN",                                // antibiotic
            [:],                                  // metadata (empty map)
            [                                     
                file("results/fastq/SRR5386858_1.fastq.gz"),
                file("results/fastq/SRR5386858_2.fastq.gz")
            ]                                     // path list (paired-end)
        )
    ) */
    
    // Each process gets the full row and extracts what it needs
    SRA_DOWNLOAD(samples_ch)
    QC(SRA_DOWNLOAD.out.fastq)
    //QC(fastq_ch)
    ASSEMBLY(QC.out.fastp)
    BLAST(ASSEMBLY.out.assembly)
    flank_script = file("bin/flank/get_flanking_regions.py")
    FLANK(BLAST.out.blast, flank_script)

    FLANK.out.flanks
        .map { accession, genus_species, antibiotic, metadata, flanks -> 
            tuple(antibiotic, tuple(accession, genus_species, flanks))
        }
        .groupTuple()
        .set { grouped_flanks }

    inference_script = file("bin/inference/inference.py")
    INFERENCE(grouped_flanks, inference_script)
}