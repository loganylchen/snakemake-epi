# Snakemake workflow: `snakemake-epitranscriptome`

[![Snakemake](https://img.shields.io/badge/snakemake-≥6.3.0-brightgreen.svg)](https://snakemake.github.io)
[![GitHub actions status](https://github.com/loganylchen/snakemake-epitranscriptome/workflows/Tests/badge.svg?branch=main)](https://github.com/loganylchen/snakemake-epitranscriptome/actions?query=branch%3Amain+workflow%3ATests)


A Snakemake workflow for `epitranscriptome analysis`


## Usage

``` bash
# Step1: Deploy the workflow (Don't forget the last .[dot])
snakedeploy deploy-workflow  --branch main  https://github.com/loganylchen/snakemake-epitranscriptome .

# Step2: Modify the samples.tsv file in the config directory
# For now, only the original gloritools was supported, so the ToolType can only be gloritools.
# Read2Fastq should be the path to the fastq file.
# As the gloritools does not need the read1 file, the Read1Fastq could be leaved as empty

# Step3: Modify the config.yaml in the config directory
# Maybe the only thing need to be modified are the references
# You'd better copy the files in the working directory, and input the relative path to these files in the config.yaml file.
# Step4: Running the workflow 
snakemake --cores THREADS --use-singularity --use-conda  

# Final results are in the `results/$SAMPLENAME/gloritools/` directory.
```

- [`samples.tsv`](config/samples.tsv) example.
- [`config.yaml`](config/config.yaml) example.

```yaml
reference:
  transcriptome_fa: transcriptome.fa
  gtf: transcript.gtf
  genome_fa: genome.fa
```

The usage of this workflow is described in the [Snakemake Workflow Catalog](https://snakemake.github.io/snakemake-workflow-catalog/?usage=loganylchen%2Fsnakemake-epitranscriptome).

If you use this workflow in a paper, don't forget to give credits to the authors by citing the URL of this (original) snakemake-epitranscriptomesitory and its DOI (see above).

