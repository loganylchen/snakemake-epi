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
transcriptome_bowtie_bam = snakemake.output.readname_sorted_bam

# temp
bowtie_raw_bam=f'{output_prefix}.rawBowtie.out.bam'
unmapped_fastq=f'{output_prefix}.rawBowtie.unmapped.fq'


cmd6=f'''
bowtie -k 1 -m 1  -v 2  \
    --best --strata -p {threads}  \
    -x {transcriptome_index} {fastq}  \
    -S {bowtie_raw_bam}  --un  {unmapped_fastq} {log}
'''
# print(cmd6)
shell(cmd6)


cmd7=f'''
samtools view -F 2324 -@ {threads} -h {bowtie_raw_bam} | samtools sort - -n -@ {threads} -o {transcriptome_bowtie_bam}
'''
# print(cmd7)
shell(cmd7)


if transcriptome_unmapped_fastq.endswith('.gz'):
    cmd = 'gzip -c {unmapped_fastq} > {transcriptome_unmapped_fastq} ; rm {unmapped_fastq}'
else:
    cmd = 'mv  {unmapped_fastq}  {transcriptome_unmapped_fastq} '
shell(cmd)