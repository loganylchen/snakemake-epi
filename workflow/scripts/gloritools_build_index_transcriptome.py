import pysam
import os
import sys
from snakemake import shell

log = snakemake.log_fmt_shell(stdout=False, stderr=True)

sys.stdout = open(snakemake.log.log,'w')
threads=snakemake.threads

converted_transcriptome_reference=snakemake.input.select_transcriptome_reference
raw_transcript_reference=snakemake.input.raw_transcript_reference
threads=snakemake.threads





cmd = f'''
bowtie-build --threads {threads} -q {converted_transcriptome_reference} {converted_transcriptome_reference} {log}
'''
print(cmd)
shell(cmd)

cmd = f'''
bowtie-build --threads {threads} -q {raw_transcript_reference} {raw_transcript_reference} {log}
'''
print(cmd)
shell(cmd)




        





