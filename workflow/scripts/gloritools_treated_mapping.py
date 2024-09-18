from snakemake import shell
import os
import pysam
import re
import json
import gzip

log = snakemake.log_fmt_shell(stdout=True, stderr=True,append=True)



def A2G_change_fastq(raw_fastq,out_fastq,info):
    change_info = dict()
    with pysam.FastxFile(raw_fastq) as fastq, gzip.open(out_fastq,'wt') as outfq:
        for entry in fastq:
            A_sites = [m.start() for m in re.finditer('A', entry.sequence)]
            entry.sequence=entry.sequence.replace('A','G')
            change_info[entry.name]=A_sites
            outfq.write(str(entry)+'\n')
    with open(info,'w') as f:
        f.write(json.dumps(change_info,indent=4))
        
# inputs
fastq=snakemake.input.fastq
threads = snakemake.threads
ag_genome_dir = os.path.dirname(snakemake.input.ag_genome_indexes[0])
ag_transcriptome_index=snakemake.input.ag_transcriptome_reference
# outputs
ag_change_fastq=snakemake.output.ag_change_fastq
info_json = snakemake.output.info_json

unmapped_fastq=snakemake.output.unmapped_fastq
transcriptome_unmapped_fastq=snakemake.output.unmapped_fastq2
output_prefix = snakemake.params.output_prefix
star_mapping_bam_step1=f'{output_prefix}.tmp.step1_sortbyname.star.bam'
star_mapping_bam_step2=f'{output_prefix}.tmp.step2_sortbyloc.star.bam'
bowtie_mapping_bam_step3=f'{output_prefix}.tmp.step3.bowtie.bam'
bowtie_mapping_bam_step4=f'{output_prefix}.tmp.step4_sortbyloc.bowtie.bam'


# prepare
A2G_change_fastq(fastq,ag_change_fastq,info_json)

# params
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

cmd2 =f'''
samtools view -F 20 -@ {threads} -h {output_prefix}.Aligned.out.bam | samtools sort -n -o {star_mapping_bam_step1}
'''
# print(cmd2)
shell(cmd2)

cmd3=f'''
mv {output_prefix}.Unmapped.out.mate1 {unmapped_fastq}
'''
# print(cmd3)
shell(cmd3)


cmd4=f'''
samtools view -F 4 -bS -@ {threads} -h {star_mapping_bam_step1} | samtools sort - -o {star_mapping_bam_step2} 
'''
# print(cmd4)
shell(cmd4)

cmd5=f'''
samtools index {star_mapping_bam_step2}
'''
# print(cmd5)
shell(cmd5)



# cmd6=f'''
# bowtie -k 1 -m 1  -v 2  \
#     --best --strata -p {threads}  \
#     -x {transcriptome_index} {unmapped_fastq}  \
#     -S {bowtie_mapping_bam_step3}  --un  {transcriptome_unmapped_fastq} {log}
# '''
# # print(cmd6)
# shell(cmd6)

# cmd7=f'''
# samtools view -F 4 -bS -@ {threads} -h {bowtie_mapping_bam_step3} | samtools sort - -o {bowtie_mapping_bam_step4} 
# '''
# # print(cmd7)
# shell(cmd7)

# cmd8=f'''
# samtools index {bowtie_mapping_bam_step4}
# '''
# # print(cmd8)
# shell(cmd8)

shell('touch {unmapped_fastq} {transcriptome_unmapped_fastq} {log}')