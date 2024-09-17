import pysam
import os
import sys
from snakemake import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True,append=True)


threads=snakemake.threads

converted_transcriptome_reference=snakemake.input.select_transcriptome_reference
raw_transcript_reference=snakemake.input.raw_transcript_reference
threads=snakemake.threads





cmd = f'''
bowtie-build --threads {threads} -q {converted_transcriptome_reference} {converted_transcriptome_reference} {log}
'''
shell(cmd)

cmd = f'''
bowtie-build --threads {threads} -q {raw_transcript_reference} {raw_transcript_reference} {log}
'''
shell(cmd)




        





