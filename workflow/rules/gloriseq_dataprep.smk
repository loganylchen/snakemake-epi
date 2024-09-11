# rule glori_trim:
#     input:
#         fastq="data/{sample}/{sample}.read2.fq.gz",
#     output:
#         clean_fastq="results/{sample}/gloritools/cleandata/{sample}_trimmed.fq.gz",
#     params:
#         params=config["gloritools"]["trim_galore"],
#         outdir=lambda w, output: os.path.dirname(output.clean_fastq),
#     threads: config["threads"]["trim_galore"]
#     conda:
#         "../envs/trim_galore.yaml"
#     log:
#         log="logs/gloritools/{sample}_trim_galore.log",
#         err="logs/gloritools/{sample}_trim_galore.err",
#     benchmark:
#         "benchmarks/gloritools/{sample}_trim_galore.txt"
#     shell:
#         " trim_galore {params.params} "
#         " -j {threads} -o {params.outdir} "
#         " --basename {wildcards.sample} "
#         " {input.fastq} 1>{log.log} 2>{log.err} "


# rule glori_rmdup:
#     input:
#         trimmed_fastq="results/{sample}/gloritools/cleandata/{sample}_trimmed.fq.gz",
#     output:
#         rmdup_fastq="results/{sample}/gloritools/cleandata/{sample}_rmdup.fq.gz",
#         dupname_f="results/{sample}/gloritools/cleandata/{sample}_rmdup.txt",
#     threads: config["threads"]["seqkit_rmdup"]
#     conda:
#         "../envs/seqkit.yaml"
#     log:
#         log="logs/gloritools/{sample}_seqkit_rmdup.log",
#         err="logs/gloritools/{sample}_seqkit_rmdup.err",
#     benchmark:
#         "benchmarks/gloritools/{sample}_seqkit_rmdup.txt"
#     shell:
#         " seqkit rmdup -j {threads} -s "
#         " -D {output.dupname_f} {input.trimmed_fastq} "
#         " -o {output.rmdup_fastq} 1>{log.log} 2>{log.err} "


rule glori_trim_dedup:
    input:
        fastq="data/{sample}/{sample}.read2.fq.gz",
    output:
        clean_fastq="results/{sample}/gloritools/cleandata/{sample}_rmdup.fq.gz",
        fastp_html="results/{sample}/gloritools/cleandata/{sample}_rmdup.html",
        fastp_json="results/{sample}/gloritools/cleandata/{sample}_rmdup.json",
    params:
        params=config["gloritools"]["fastp"],
    threads: config["threads"]["fastp"]
    conda:
        "../envs/fastp.yaml"
    log:
        log="logs/gloritools/{sample}_fastp.log",
        err="logs/gloritools/{sample}_fastp.err",
    benchmark:
        "benchmarks/gloritools/{sample}_fastp.txt"
    shell:
        " fastp --in1 {input.fastq} "
        " --out1 {output.clean_fastq} "
        " --json {output.fastp_json} "
        " --html {output.fastp_html} "
        " --thread {threads} 1>{log.log} 2>{log.err} "


rule glori_trim_umi:
    input:
        rmdup_fastq="results/{sample}/gloritools/cleandata/{sample}_rmdup.fq.gz",
    output:
        rmumi_fastq="results/{sample}/gloritools/cleandata/{sample}_rmumi.fq",
    params:
        params=config["gloritools"]["fastx_trimmer"],
    threads: config["threads"]["fastx_trimmer"]
    conda:
        "../envs/fastx_toolkit.yaml"
    log:
        log="logs/gloritools/{sample}_fastx_trimmer.log",
        err="logs/gloritools/{sample}_fastx_trimmer.err",
    benchmark:
        "benchmarks/gloritools/{sample}_fastx_trimmer.txt"
    shell:
        "fastx_trimmer {params.params} -i {input.rmdup_fastq} -o {output.rmumi_fastq} 1>{log.log} 2>{log.err}"
