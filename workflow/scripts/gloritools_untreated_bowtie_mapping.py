from snakemake import shell
import os

log = snakemake.log_fmt_shell(stdout=True, stderr=True,append=True)


# inputs
fastq=snakemake.input.fastq
threads = snakemake.threads
transcriptome_index=snakemake.input.transcriptome_reference


# prefix
output_prefix = snakemake.params.output_prefix

# outputs

transcriptome_unmapped_fastq=snakemake.output.transcriptome_unmapped_fastq
transcriptome_bowtie_bam = snakemake.output.transcriptome_bowtie_bam

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
samtools view -F 4 -bS -@ {threads} -h {bowtie_raw_bam} | samtools sort - -o {transcriptome_bowtie_bam} 
'''
# print(cmd7)
shell(cmd7)

cmd8=f'''
samtools index {transcriptome_bowtie_bam}
'''
# print(cmd8)
shell(cmd8)

if transcriptome_unmapped_fastq.endswith('.gz'):
    cmd = 'gzip -c {unmapped_fastq} > {transcriptome_unmapped_fastq} ; rm {unmapped_fastq}'
else:
    cmd = 'mv  {unmapped_fastq}  {transcriptome_unmapped_fastq} '
shell(cmd)