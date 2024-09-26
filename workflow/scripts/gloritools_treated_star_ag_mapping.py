from snakemake import shell
import os
import pysam
import re
import json
import pgzip
from Bio.Seq import reverse_complement


log = snakemake.log_fmt_shell(stdout=True, stderr=True,append=True)
# inputs
ag_change_fastq=snakemake.input.ag_change_fastq
threads = snakemake.threads
ag_genome_dir = os.path.dirname(snakemake.input.ag_genome_indexes[0])

# outputs

ag_genome_unmapped_fastq = snakemake.output.ag_genome_unmapped_fastq
readname_sorted_star_bam=snakemake.output.ag_genome_readname_sorted_bam


# prefix
output_prefix = snakemake.params.output_prefix

# temp
raw_star_bam=f"{output_prefix}.Aligned.out.bam"


extra_para=' --readFilesCommand zcat ' if ag_change_fastq.endswith('.gz') else ''

cmd1=f'''
STAR --runThreadN {threads} \
    --genomeDir {ag_genome_dir} \
    --limitOutSJcollapsed 5000000 \
    --outFilterMismatchNmax 2 \
    --outFilterScoreMinOverLread 0.5 \
    --outFilterMatchNminOverLread 0.5 \
    --seedSearchStartLmax 30 \
    --outSAMattributes All --outSAMprimaryFlag AllBestScore --outMultimapperOrder Random --outSAMmultNmax 1 --outSAMtype BAM Unsorted \
    --outFilterMultimapNmax 1 {extra_para} \
    --outFileNamePrefix {output_prefix}.  --readFilesIn {ag_change_fastq} \
    --outSAMunmapped Within --outReadsUnmapped Fastx {log}
'''
# print(cmd1)
shell(cmd1)

shell('echo "`date`|STAR MAPPING DONE" {log}')
# read unmapped (0x4)
# read reverse strand (0x10)
# not primary alignment (0x100)
# supplementary alignment (0x800)
# 2324
cmd2 =f'''
samtools view -F 2324 -@ {threads} -h {raw_star_bam} | samtools sort -n -@ {threads} -o {readname_sorted_star_bam}
'''
# print(cmd2)
shell(cmd2)

shell('echo "`date`|FILTERING & SORTING  DONE" {log}')


if ag_genome_unmapped_fastq.endswith('.gz'):
    cmd = 'gzip -c {output_prefix}.Unmapped.out.mate1 > {ag_genome_unmapped_fastq} ; rm {output_prefix}.Unmapped.out.mate1'
else:
    cmd = 'mv {output_prefix}.Unmapped.out.mate1 {ag_genome_unmapped_fastq} '
shell(cmd)






