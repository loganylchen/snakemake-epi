import os
import re
import json
import pgzip
from Bio.Seq import reverse_complement
import time
import pysam


sys.stdout = open(snakemake.log.log,'w') 
sys.stderr= open(snakemake.log.err, 'w')

# inputs
fastq=snakemake.input.fastq
threads = snakemake.threads

# outputs
ag_change_fastq=snakemake.output.ag_change_fastq
info_json = snakemake.output.info_json


BATCH_SIZE=5000000

def _format_seconds(sec):
    hours = sec // 3600  
    minutes = (sec % 3600) // 60  
    seconds = sec % 60  
    return hours, minutes, seconds


def A2G_change_fastq(raw_fastq,out_fastq,info,threads=threads):
    start_time = time.time()
    print("Start:", time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
    change_info = dict()
    n=0
    with pysam.FastxFile(raw_fastq) as fastq, pgzip.open(out_fastq,'wt',thread=threads, blocksize=2*10**8) as outfq:
        for entry in fastq:
            n+=1
            A_sites = [m.start() for m in re.finditer('A', entry.sequence)]
            entry.sequence=entry.sequence.replace('A','G')
            change_info[entry.name]=A_sites
            outfq.write(str(entry)+'\n')
            if n %BATCH_SIZE ==0:
                end_time = time.time()
                h,m,s = _format_seconds(end_time-start_time)
                print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
                print(f'{n} reads has been converted, the last {BATCH_SIZE} reads used: {h}:{m}:{s}')
                start_time=time.time()

    end_time = time.time()
    h,m,s = _format_seconds(end_time-start_time)
    print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
    print(f'{n} reads has been converted, the last {BATCH_SIZE} reads used: {h}:{m}:{s}')
    start_time=time.time()
    
    with open(info,'w') as f:
        f.write(json.dumps(change_info,indent=4))
    print("End:", time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))


START_TIME= time.time()

A2G_change_fastq(fastq,ag_change_fastq,info_json)

END_TIME = time.time()
h,m,s=_format_seconds(END_TIME-START_TIME)
print(f"ALL Program Takes: {h}:{m}:{s} ")

