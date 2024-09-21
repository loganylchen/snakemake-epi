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
                read.query_sequence = _multi_replace(read.query_sequence,info_dict[read.query_name])
                read.query_qualities = qualities
                previous_read_name =read.query_name
                output_bam.write(read)
    print('Done the convertion')





# inputs
fastq=snakemake.input.fastq
threads = snakemake.threads
ag_genome_dir = os.path.dirname(snakemake.input.ag_genome_indexes[0])
ag_rvs_genome_dir = os.path.dirname(snakemake.input.ag_rvs_genome_indexes[0])
ag_transcriptome_index=snakemake.input.ag_transcriptome_reference


# outputs
ag_change_fastq=snakemake.output.ag_change_fastq
info_json = snakemake.output.info_json
ag_genome_unmapped_fastq = snakemake.output.ag_genome_unmapped_fastq
ag_rvs_genome_unmapped_fastq = snakemake.output.ag_rvs_genome_unmapped_fastq
ag_transcriptome_unmapped_fastq=snakemake.output.ag_transcriptome_unmapped_fastq



star_mapping_bam=f"{output_prefix}.Aligned.out.bam"
step1_star_mapping_namesorted_bam=f'{output_prefix}.tmp.step1_sortbyname.star.bam'
step2_star_mapping_ag_converted_bam=f'{output_prefix}.tmp.step2_agconverted.star.bam'
step3_star_mapping_locsorted_bam=f'{output_prefix}.tmp.step3_sortbyloc.star.bam'

star_mapping_rvs_bam=f"{output_prefix}.rvs.Aligned.out.bam"

step4_star_mappingrvs_namesorted_bam=f'{output_prefix}.tmp.step4_sortbyname.starrvs.bam'
step5_star_mappingrvs_ag_converted_bam=f'{output_prefix}.tmp.step5_agconverted.starrvs.bam'
step6_star_mappingrvs_locsorted_bam=f'{output_prefix}.tmp.step6_sortbyloc.starrvs.bam'



step7_bowtie_mapping_bam=f'{output_prefix}.tmp.step4.bowtie.bam'
step8_bowtie_mapping_bam=f'{output_prefix}.tmp.step5_sortbyloc.bowtie.bam'

output_prefix = snakemake.params.output_prefix


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
# read unmapped (0x4)
# read reverse strand (0x10)
# not primary alignment (0x100)
# supplementary alignment (0x800)
# 2324
cmd2 =f'''
samtools view -F 2324 -@ {threads} -h {star_mapping_bam} | samtools sort -n -o {step1_star_mapping_namesorted_bam}
'''
# print(cmd2)
shell(cmd2)


recover_A(step1_star_mapping_namesorted_bam,step2_star_mapping_ag_converted_bam,info_json)


cmd3=f'''
mv {output_prefix}.Unmapped.out.mate1 {ag_genome_unmapped_fastq}
'''
# print(cmd3)
shell(cmd3)


cmd4=f'''
samtools view -F 4 -bS -@ {threads} -h {step2_star_mapping_ag_converted_bam} | samtools sort - -o {step3_star_mapping_locsorted_bam} 
'''
# print(cmd4)
shell(cmd4)

cmd5=f'''
samtools index {step3_star_mapping_locsorted_bam}
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