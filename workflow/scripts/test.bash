#!/usr/bin/env bash

exec 2> "${snakemake_log[0]}"  # send all stderr from this script to the log file

reads=(${snakemake_input[reads]})  # don't double-quote this - we want word splitting

r1="${reads[0]}"
r2="${reads[1]}"

bwa index "${snakemake_input[reference]}"
bwa mem ${snakemake_params[opts]} -t ${snakemake[threads]} \
    "${snakemake_input[reference]}" "$r1" "$r2" > "${snakemake_output[0]}"