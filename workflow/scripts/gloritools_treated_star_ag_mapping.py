from snakemake import shell
import os
import pysam
import re
import json
import gzip
from Bio.Seq import reverse_complement


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


def _multi_replace(string,index_list,nucleotide='A'):
    string_list = [x for x in string]
    for idx in index_list:
        string_list[idx]=nucleotide
    return ''.join(string_list)

def recover_A(readname_sorted_bam,output_bam,index_json): 
    with open(index_json) as f:
        info_dict = json.load(f)
    
    with pysam.AlignmentFile(readname_sorted_bam) as input_bam, pysam.AlignmentFile(output_bam,'wb',template=input_bam) as output_bam:
        previous_read_name=None
        for read in input_bam:
            if read.query_name == previous_read_name:
                raise ValueError(f'{read.query_name} was duplicated in the {readname_sorted_bam}')
            else:
                qualities = read.query_qualities
                if read.is_forward:
                    read.query_sequence = _multi_replace(read.query_sequence,info_dict[read.query_name])
                else:
                    read.query_sequence = reverse_complement(_multi_replace(read.get_forward_sequence(),info_dict[read.query_name]))
                read.query_qualities = qualities
                previous_read_name =read.query_name
                output_bam.write(read)
    print('Done the convertion')

# inputs
fastq=snakemake.input.fastq
threads = snakemake.threads
ag_genome_dir = os.path.dirname(snakemake.input.ag_genome_indexes[0])

# outputs
ag_change_fastq=snakemake.output.ag_change_fastq
info_json = snakemake.output.info_json
ag_genome_unmapped_fastq = snakemake.output.ag_genome_unmapped_fastq
ag_genome_star_bam=snakemake.output.ag_genome_star_bam


# prefix
output_prefix = snakemake.params.output_prefix

# temp
raw_star_bam=f"{output_prefix}.Aligned.out.bam"
readname_sorted_bam=f"{output_prefix}.Aligned.ReadnameSorted.bam"
convered_bam=f"{output_prefix}.Aligned.Converted.bam"


# prepare
shell('echo "`date`|A2G_change_fastq START" {log}')
A2G_change_fastq(fastq,ag_change_fastq,info_json)
shell('echo "`date`|A2G_change_fastq DONE" {log}')
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

shell('echo "`date`|STAR MAPPING DONE" {log}')
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

shell('echo "`date`|FILTERING & SORTING  DONE" {log}')


recover_A(readname_sorted_bam,convered_bam,info_json)

shell('echo "`date`|RECOVER A  DONE" {log}')

cmd3=f'''
mv {output_prefix}.Unmapped.out.mate1 {ag_genome_unmapped_fastq}
'''
# print(cmd3)
shell(cmd3)


cmd4=f'''
samtools view -F 4 -bS -@ {threads} -h {convered_bam} | samtools sort -@ {threads} - -o {ag_genome_star_bam} 
'''
# print(cmd4)
shell(cmd4)

cmd5=f'''
samtools index {ag_genome_star_bam}
'''
# print(cmd5)
shell(cmd5)


