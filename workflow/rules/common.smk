import pandas as pd
import sys


def get_output_list_for_one_sample(sample):
    return [
        # f"data/{sample}/fastq/pass.fq.gz",
        f"results/{sample}/gloritools/glori/{sample}.totalm6A.FDR.csv",
        f"results/{sample}/gloritools/glori_asControl/{sample}.totalm6A.FDR.csv",
    ]


samples = {"SRR21356250": ""}


def get_final_output():
    final_output = []
    for sample in samples.keys():
        final_output += get_output_list_for_one_sample(sample)
    return final_output
