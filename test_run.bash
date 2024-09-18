#!/bin/bash

snakemake  --cores 16  --use-conda  --rerun-incomplete -k -p   --directory .test --snakefile workflow/Snakefile --show-failed-logs