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
transcriptome_index=snakemake.input.ag_transcriptome_reference
info_json=snakemake.input.info_json

# prefix
output_prefix = snakemake.params.output_prefix

# outputs

transcriptome_unmapped_fastq=snakemake.output.ag_transcriptome_fastq
transcriptome_bowtie_bam = snakemake.output.ag_transcriptome_bowtie_bam

# temp
bowtie_raw_bam=f'{output_prefix}.rawBowtie.out.bam'
bowtie_sortname_bam=f'{output_prefix}.rawBowtie.sortname.bam'
bowtie_converted_bam=f'{output_prefix}.rawBowtie.converted.bam'

cmd6=f'''
bowtie -k 1 -m 1  -v 2  \
    --best --strata -p {threads}  \
    -x {transcriptome_index} {fastq}  \
    -S {bowtie_raw_bam}  --un  {transcriptome_unmapped_fastq} {log}
'''
# print(cmd6)
shell(cmd6)




cmd7=f'''
samtools view -F 2324 -@ {threads} -h {bowtie_raw_bam} | samtools sort - -n -@ {threads} -o {bowtie_sortname_bam}
'''
# print(cmd7)
shell(cmd7)
recover_A(bowtie_sortname_bam,bowtie_converted_bam,info_json)

cmd8=f'''
samtools view -F 2324 -@ {threads} -h {bowtie_converted_bam} | samtools sort -@ {threads} - -o {transcriptome_bowtie_bam}
'''
shell(cmd8)
cmd9=f'''
samtools index {transcriptome_bowtie_bam}
'''
# print(cmd8)
shell(cmd9)

