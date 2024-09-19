from snakemake import shell
import os
import pysam
import re
import json
import gzip
from Bio.Seq import reverse_complement

log = snakemake.log_fmt_shell(stdout=True, stderr=True,append=True)



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
rvs_genome_dir = os.path.dirname(snakemake.input.rvs_genome_indexes[0])
info_json=snakemake.input.info_json
# outputs

rvs_genome_unmapped_fastq = snakemake.output.rvs_genome_unmapped_fastq
rvs_genome_star_bam=snakemake.output.rvs_genome_star_bam


# prefix
output_prefix = snakemake.params.output_prefix

# temp
raw_star_bam=f"{output_prefix}.rvs.Aligned.out.bam"
readname_sorted_bam=f"{output_prefix}.rvs.Aligned.ReadnameSorted.bam"
convered_bam=f"{output_prefix}.rvs.Aligned.Converted.bam"




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


recover_A(readname_sorted_bam,convered_bam,info_json)


cmd3=f'''
mv {output_prefix}.rvs.Unmapped.out.mate1 {rvs_genome_unmapped_fastq}
'''
# print(cmd3)
shell(cmd3)


cmd4=f'''
samtools view -F 4 -bS -@ {threads} -h {convered_bam} | samtools sort -@ {threads} - -o {rvs_genome_star_bam} 
'''
# print(cmd4)
shell(cmd4)

cmd5=f'''
samtools index {rvs_genome_star_bam}
'''
# print(cmd5)
shell(cmd5)


