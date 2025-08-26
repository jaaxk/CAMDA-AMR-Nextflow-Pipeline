# Nextflow implementation of 2nd place [CAMDA AMR Prediction Challenge Pipeline](https://github.com/jaaxk/CAMDA_AMR)
## Running inference on WGS samples
1. Clone repository
2. Download [models]() and place in bin/resources/models (models to be released soon)
3. Generate `samples.tsv` with header `accession  genus_species  antibiotic`
   `Example:
   accession	genus_species	antibiotic
   SRR5827165	neisseria_gonorrhoeae	TET
   SRR5827255	neisseria_gonorrhoeae	TET`
4. Place `samples.tsv` in `data/`
5. Modify `nextflow.config` to your configuration
   Current implementation launches separate SLURM jobs for each sample
6. See per-antibiotic results in `results/inference`

For more information on this project, visit the [official repository](https://github.com/jaaxk/CAMDA_AMR)
Please note: This pipeline and model is currently under consideration for a patent
