import pysam
import os
import sys
from snakemake import shell
import numpy as np

log = snakemake.log_fmt_shell(stdout=False, stderr=True)

sys.stdout = open(snakemake.log.log,'w')
threads=snakemake.threads

convert_genome_reference=snakemake.input.convert_genome_reference
rev_convert_genome_reference=snakemake.input.rev_convert_genome_reference
raw_genome_reference=snakemake.input.raw_genome_reference

convert_genome_prefix=snakemake.params.convert_genome_prefix
rev_convert_genome_prefix=snakemake.params.rev_convert_genome_prefix
raw_genome_prefix=snakemake.params.raw_genome_prefix



cmd = f'''
STAR --runMode genomeGenerate -runThreadN {threads} \
        --genomeDir {rev_convert_genome_prefix}  \
        --genomeFastaFiles {rev_convert_genome_reference} \
        --genomeSAindexNbases {int(round(min(14, np.log2(sum(pysam.FastaFile(rev_convert_genome_reference).lengths))/2 - 1)))} \
        --limitGenomeGenerateRAM 84807429045 {log}
'''
print(cmd)
shell(cmd)

cmd = f'''
STAR --runMode genomeGenerate -runThreadN {threads} \
        --genomeDir {convert_genome_prefix}  \
        --genomeFastaFiles {convert_genome_reference} \
        --genomeSAindexNbases {int(round(min(14, np.log2(sum(pysam.FastaFile(convert_genome_reference).lengths))/2 - 1)))} \
        --limitGenomeGenerateRAM 84807429045 {log}
'''
print(cmd)
shell(cmd)

cmd = f'''
STAR --runMode genomeGenerate -runThreadN {threads} \
        --genomeDir {raw_genome_prefix}  \
        --genomeFastaFiles {raw_genome_reference} \
        --genomeSAindexNbases {int(round(min(14, np.log2(sum(pysam.FastaFile(raw_genome_reference).lengths))/2 - 1)))} \
        --limitGenomeGenerateRAM 84807429045 {log}
'''
print(cmd)
shell(cmd)


        
