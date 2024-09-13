rule glori_trim_dedup_umi:
    input:
        unpack(get_raw_fastq),
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
