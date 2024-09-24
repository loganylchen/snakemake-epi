import os
import re
import json
import pysam
from Bio.Seq import reverse_complement
import time


sys.stdout = open(snakemake.log.log,'w') 
sys.stderr= open(snakemake.log.err, 'w')

# inputs
readname_sorted_bam=snakemake.input.readname_sorted_bam
info_json = snakemake.input.info_json

# outputs
output_bam=snakemake.output.output_bam


def _format_seconds(sec):
    hours = sec // 3600  
    minutes = (sec % 3600) // 60  
    seconds = sec % 60  
    return hours, minutes, seconds


def _multi_replace(string,index_list,nucleotide='A'):
    string_list = [x for x in string]
    for idx in index_list:
        string_list[idx]=nucleotide
    return ''.join(string_list)

def recover_A(readname_sorted_bam,output_bam,index_json): 
    print("Start loading json file:", time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
    with open(index_json) as f:
        info_dict = json.load(f)
    print("Done loading json file:", time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
    
    n=0
    start_time = time.time()
    with pysam.AlignmentFile(readname_sorted_bam) as input_bam, pysam.AlignmentFile(output_bam,'wb',template=input_bam) as output_bam:
        previous_read_name=None
        for read in input_bam:
            n+=1
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
            if n %10000 ==0:
                end_time = time.time()
                h,m,s = _format_seconds(end_time-start_time)
                print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
                print(f'{n} reads has been reconverted, the last 10000 reads used: {h}:{m}:{s}')
                start_time=time.time()
    print("End:", time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))





START_TIME= time.time()

recover_A(readname_sorted_bam,output_bam,info_json)

END_TIME = time.time()
h,m,s=_format_seconds(END_TIME-START_TIME)
print(f"ALL Program Takes: {h}:{m}:{s} ")
