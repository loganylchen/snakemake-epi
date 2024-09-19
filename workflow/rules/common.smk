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
    elif samples[wildcards.sample]["ToolType"] == "patch-gloritools":
        return {"fastq": samples[wildcards.sample]["Read2Fastq"]}
    else:
        raise ValueError(f'{samples[wildcards.sample]["ToolType"]} was not supported')


def get_output_list_for_one_sample(sample):
    if samples[sample]["ToolType"] == "gloritools":
        return [
            f"results/{sample}/gloritools/glori/{sample}.totalm6A.FDR.csv",
            # f"results/{sample}/gloritools/glori_asControl/{sample}.totalm6A.FDR.csv",
        ]
    elif samples[sample]["ToolType"] == "patch-gloritools":
        return [
            f"results/{sample}/gloritools/cleandata/{sample}_rmdup.json",
            f"results/{sample}/gloritools/untreated/{sample}.star.bam",
            f"results/{sample}/gloritools/untreated/{sample}.bowtie.bam",
            f"results/{sample}/gloritools/treated/{sample}.star.ag.bam",
            f"results/{sample}/gloritools/treated/{sample}.star.rvs.bam",
            f"results/{sample}/gloritools/treated/{sample}.bowtie.ag.bam",
        ]
    else:
        raise ValueError(
            f'|{samples[sample]["ToolType"]}| was not supported, {samples[sample]["ToolType"] in["gloritools", "patch-gloritools"]}'
        )


def get_final_output():
    final_output = [
        "references/gloritools/annotation.tbl",
        "references/gloritools/transcriptome_AG.fa",
        "references/gloritools/genome_AG.fa",
        "references/gloritools/genome_rc_AG.fa",
        "references/gloritools/genome/Log.out",
        "references/gloritools/transcriptome_AG.fa.rev.1.ebwt",
    ]
    for sample in samples.keys():
        final_output += get_output_list_for_one_sample(sample)
    return final_output
