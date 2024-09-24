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
rvs_genome_dir = os.path.dirname(snakemake.input.rvs_genome_indexes[0])

# outputs

rvs_genome_unmapped_fastq = snakemake.output.rvs_genome_unmapped_fastq
readname_sorted_bam=snakemake.output.readname_sorted_bam


# prefix
output_prefix = snakemake.params.output_prefix

# temp
raw_star_bam=f"{output_prefix}.rvs.Aligned.out.bam"

# params
extra_para=' --readFilesCommand zcat ' if fastq.endswith('.gz') else ''

cmd1=f'''
STAR --runThreadN {threads} \
    --genomeDir {rvs_genome_dir} \
    --limitOutSJcollapsed 5000000 \
    --outFilterMismatchNmax 2 \
    --outFilterScoreMinOverLread 0.5 \
    --outFilterMatchNminOverLread 0.5 \
    --seedSearchStartLmax 30 \
    --outSAMattributes All --outSAMprimaryFlag AllBestScore --outMultimapperOrder Random --outSAMmultNmax 1 --outSAMtype BAM Unsorted \
    --outFilterMultimapNmax 1 {extra_para} \
    --outFileNamePrefix {output_prefix}.rvs.  --readFilesIn {fastq} \
    --outSAMunmapped Within --outReadsUnmapped Fastx {log}
'''
# print(cmd1)
shell(cmd1)
# read unmapped (0x4)
# read reverse strand (0x10)
# not primary alignment (0x100)
# supplementary alignment (0x800)
# 2324
cmd2 =f'''
samtools view -F 2324 -@ {threads} -h {raw_star_bam} | samtools sort -n -@ {threads} -o {readname_sorted_bam}
'''
# print(cmd2)
shell(cmd2)


