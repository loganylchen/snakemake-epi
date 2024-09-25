import os
import re
import json
import pysam
import pickle
import sqlite3
from Bio.Seq import reverse_complement
import time


sys.stdout = open(snakemake.log.log,'w') 
sys.stderr= open(snakemake.log.err, 'w')

# inputs
readname_sorted_bam=snakemake.input.readname_sorted_bam
info_db = snakemake.input.info_db

# outputs
output_bam=snakemake.output.output_bam
BATCH_SIZE=5000000

def _format_seconds(sec):
    hours = int(sec // 3600 ) 
    minutes = int((sec % 3600) // 60)  
    seconds = int(sec % 60 ) 
    return hours, minutes, seconds


def _multi_replace(string,index_list,nucleotide='A'):
    string_list = [x for x in string]
    for idx in index_list:
        string_list[idx]=nucleotide
    return ''.join(string_list)

def recover_A(readname_sorted_bam,output_bam,info_db): 
    print("Start loading json file:", time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
    conn = sqlite3.connect(info_db)
    c = conn.cursor()
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
                c.execute('SELECT value FROM reads WHERE name = ? LIMIT 1', (read,query_name,))
                result = c.fetchone()
                index_list = pickle.loads(result[0])
                if read.is_forward:
                    read.query_sequence = _multi_replace(read.query_sequence,index_list)
                else:
                    read.query_sequence = reverse_complement(_multi_replace(read.get_forward_sequence(),index_list))
                read.query_qualities = qualities
                previous_read_name =read.query_name
                output_bam.write(read)
            if n %BATCH_SIZE ==0:
                end_time = time.time()
                h,m,s = _format_seconds(end_time-start_time)
                print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
                print(f'{n} reads has been reconverted, the last {BATCH_SIZE} reads used: {h}:{m}:{s}')
                start_time=time.time()
    print("End:", time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))






START_TIME= time.time()

recover_A(readname_sorted_bam,output_bam,info_db)

END_TIME = time.time()
h,m,s=_format_seconds(END_TIME-START_TIME)
print(f"ALL Program Takes: {h}:{m}:{s} ")
