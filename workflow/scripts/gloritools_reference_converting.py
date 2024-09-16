import sys
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Seq import reverse_complement
import numpy as np
import pysam
import re

sys.stdout = open(snakemake.log.log,'w') 
sys.stderr= open(snakemake.log.err, 'w')

genome_output_dir = os.path.dirname(snakemake.output.genome_reference)
transcriptome_output_dir = os.path.dirname(snakemake.output.transcriptome_reference)

os.makedirs(transcriptome_output_dir,exist_ok=True)
os.makedirs(genome_output_dir,exist_ok=True)

raw_genome_reference = snakemake.input.genome
raw_transcriptome_reference = snakemake.input.transcriptome

print('Coverting genome')
with open(snakemake.output.genome_reference,'w') as out:
    for record in SeqIO.parse(raw_genome_reference, "fasta"):
        record.seq = Seq(re.sub('A','G',str(record.seq).upper()))
        reversed_seq = reverse_complement(record.seq)
        print('Genome',record.id)
        record.id = record.id + "_AG_converted"
        SeqIO.write(record, out, "fasta")
        
print('Coverting genome reverse complement')
with open(snakemake.output.genome_reference_complement,'w') as out:
    for record in SeqIO.parse(raw_genome_reference, "fasta"):
        reversed_seq = reverse_complement(record.seq)
        converted_seq = re.sub('A', 'G', str(reversed_seq).upper())
        reversed_seq2 = reverse_complement(converted_seq)
        record.seq = Seq(reversed_seq2)
        record.id = record.id
        print('reversecomplement Genome',record.id)
        record.id = record.id + "_AG_converted"
        SeqIO.write(record, out, "fasta")
        


print('Coverting Transcritome')
with open(snakemake.output.transcriptome_reference,'w') as out:
        for record in SeqIO.parse(raw_transcriptome_reference, "fasta"):
            record.seq = Seq(re.sub('A','G',str(record.seq).upper()))
            record.id = record.id + "_AG_converted"
            SeqIO.write(record, out, "fasta")




