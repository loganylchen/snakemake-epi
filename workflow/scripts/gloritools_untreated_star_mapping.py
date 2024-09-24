from snakemake import shell
import os

log = snakemake.log_fmt_shell(stdout=True, stderr=True,append=True)


# inputs
fastq=snakemake.input.fastq
threads = snakemake.threads
genome_dir = os.path.dirname(snakemake.input.genome_indexes[0])

# output_prefix
output_prefix = snakemake.params.output_prefix

# outputs
genome_unmapped_fastq=snakemake.output.genome_unmapped_fastq
genome_star_bam=f'{output_prefix}.star.bam'

# temp
star_raw_bam=f'{output_prefix}.Aligned.out.bam'


extra_para=' --readFilesCommand zcat ' if fastq.endswith('.gz') else ''
cmd1=f'''
STAR --runThreadN {threads} \
    --genomeDir {genome_dir} \
    --limitOutSJcollapsed 5000000 \
    --outFilterMismatchNmax 2 \
    --outFilterScoreMinOverLread 0.5 \
    --outFilterMatchNminOverLread 0.5 \
    --seedSearchStartLmax 30 \
    --outSAMattributes All --outSAMprimaryFlag AllBestScore \
    --outMultimapperOrder Random --outSAMmultNmax 1 --outSAMtype BAM Unsorted \
    --outFilterMultimapNmax 1 {extra_para} \
    --outFileNamePrefix {output_prefix}.  --readFilesIn {fastq} \
    --outSAMunmapped Within --outReadsUnmapped Fastx {log}
'''

shell(cmd1)

shell('echo "`date`|STAR MAPPING DONE" {log}')

cmd2 =f'''
samtools view -F 4 -@ {threads} -h {star_raw_bam} | samtools sort -@ {threads} -o {genome_star_bam} 
'''
shell(cmd2)

shell('echo "`date`|FILTERING & SORTING DONE" {log}')

cmd3=f'''
mv {output_prefix}.Unmapped.out.mate1 {genome_unmapped_fastq} {log}
'''
shell(cmd3)
cmd4=f'''
samtools index {genome_star_bam} {log}
'''
shell(cmd4)


