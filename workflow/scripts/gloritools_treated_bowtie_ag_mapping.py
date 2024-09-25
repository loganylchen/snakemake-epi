from snakemake import shell
import os
import pysam
import re
import json
import gzip
from Bio.Seq import reverse_complement


log = snakemake.log_fmt_shell(stdout=True, stderr=True,append=True)





# inputs
fastq=snakemake.input.fastq
threads = snakemake.threads
transcriptome_index=snakemake.input.ag_transcriptome_reference


# prefix
output_prefix = snakemake.params.output_prefix

# outputs

transcriptome_unmapped_fastq=snakemake.output.ag_transcriptome_fastq
transcriptome_bowtie_bam = snakemake.output.ag_transcriptome_bowtie_bam

# temp
bowtie_raw_bam=f'{output_prefix}.rawBowtie.out.bam'



cmd6=f'''
bowtie -k 1 -m 1  -v 2  \
    --best --strata -p {threads}  \
    -x {transcriptome_index} {fastq}  \
    -S {bowtie_raw_bam}  --un  {transcriptome_unmapped_fastq} {log}
'''
# print(cmd6)
shell(cmd6)


cmd7=f'''
samtools view -F 2324 -@ {threads} -h {bowtie_raw_bam} | samtools sort - -n -@ {threads} -o {transcriptome_bowtie_bam}
'''
# print(cmd7)
shell(cmd7)
