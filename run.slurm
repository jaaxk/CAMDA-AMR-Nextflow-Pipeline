#!/bin/bash
#SBATCH --partition=long-40core
#SBATCH --output=res.txt
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mail-user=jack.vaska@stonybrook.edu
#SBATCH --mail-type=ALL


nextflow run main.nf -profile slurm
