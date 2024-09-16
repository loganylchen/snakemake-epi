#!/usr/bin/env bash
set -x 
set -e

exec 1> "${snakemake_log[0]}"
exec 2> "${snakemake_log[1]}"  # send all stderr from this script to the log file

converted_transcriptome_reference="${snakemake_input[converted_transcriptome_reference]}"
threads="${snakemake[threads]}"

bowtie-build --threads ${threads} -q ${converted_transcriptome_reference} ${converted_transcriptome_reference} 


