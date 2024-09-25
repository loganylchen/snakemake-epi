import os
import re
import pickle
import sqlite3
from Bio.Seq import reverse_complement
import time
import pysam
import pgzip

sys.stdout = open(snakemake.log.log,'w') 
sys.stderr= open(snakemake.log.err, 'w')

# inputs
fastq=snakemake.input.fastq
threads = snakemake.threads

# outputs
ag_change_fastq=snakemake.output.ag_change_fastq
info_db = snakemake.output.info_db






BATCH_SIZE=5000000

def _format_seconds(sec):
    hours = int(sec // 3600 ) 
    minutes = int((sec % 3600) // 60 ) 
    seconds = int(sec % 60)
    return hours, minutes, seconds







def A2G_change_fastq(raw_fastq,out_fastq,info,threads=threads):
    start_time = time.time()
    print("Start:", time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
    if os.path.exists(info):
        print(f'Existing {info}: removing it!!')
        os.remove(info)
    os.makedirs(os.path.dirname(info),exist_ok=True)
    conn = sqlite3.connect(info)
    c = conn.cursor()
    c.execute('CREATE TABLE IF NOT EXISTS reads (name TEXT, value BLOB)')
    data_to_insert = []
    n=0
    with pysam.FastxFile(raw_fastq) as fastq, pgzip.open(out_fastq,'wt',thread=threads, blocksize=2*10**8)  as outfq:
        for entry in fastq:
            n+=1
            A_sites = [m.start() for m in re.finditer('A', entry.sequence)]
            entry.sequence=entry.sequence.replace('A','G')
            # change_info[entry.name]=A_sites
            data_to_insert.append((entry.name,pickle.dumps(A_sites)))
            outfq.write(str(entry)+'\n')
            if n %BATCH_SIZE ==0:
                conn.executemany('INSERT INTO reads (name, value) VALUES (?, ?)', data_to_insert)
                conn.commit()
                data_to_insert=[]
                end_time = time.time()
                h,m,s = _format_seconds(end_time-start_time)
                print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
                print(f'{n} reads has been converted, the last {BATCH_SIZE} reads used: {h}:{m}:{s}')
                start_time=time.time()
                # print(sys.getsizeof(change_info))

    end_time = time.time()
    h,m,s = _format_seconds(end_time-start_time)
    print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
    print(f'{n} reads has been converted, the last {BATCH_SIZE} reads used: {h}:{m}:{s}')
    conn.executemany('INSERT INTO reads (name, value) VALUES (?, ?)', data_to_insert)
    conn.commit()
    conn.execute('PRAGMA synchronous = OFF')
    conn.execute('PRAGMA journal_mode = WAL')
    conn.execute('PRAGMA cache_size = 10000')
    c.execute('CREATE INDEX IF NOT EXISTS idx_name ON reads(name)')
    c.close()
    conn.close()


START_TIME= time.time()

A2G_change_fastq(fastq,ag_change_fastq,info_db)

END_TIME = time.time()
h,m,s=_format_seconds(END_TIME-START_TIME)
print(f"ALL Program Takes: {h}:{m}:{s} ")

