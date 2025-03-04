import pandas as pd
import sys

samples = (
    pd.read_csv(config["samples"], sep="\t", dtype={"SampleName": str}, comment="#")
    .set_index("SampleName", drop=False)
    .sort_index()
    .T.to_dict()
)


def get_raw_fastq(wildcards):
    if samples[wildcards.sample]["ToolType"] == "gloritools":
        return {"fastq": samples[wildcards.sample]["Read2Fastq"]}
    else:
        raise ValueError(f'{samples[wildcards.sample]["ToolType"]} was not supported')


def get_output_list_for_one_sample(sample):
    if samples[sample]["ToolType"] == "gloritools":
        return [
            f"results/{sample}/gloritools/glori/{sample}.totalm6A.FDR.csv",
            # f"results/{sample}/gloritools/glori_asControl/{sample}.totalm6A.FDR.csv",
        ]
    else:
        raise ValueError(f'{samples[sample]["ToolType"]} was not supported')


def get_final_output():
    final_output = []
    for sample in samples.keys():
        final_output += get_output_list_for_one_sample(sample)
    return final_output
