rule gloritools_genome_mapping:
    input:
        fastq="results/{sample}/gloritools/cleandata/{sample}_rmdup.fq.gz",
        ag_genome_dir=multiext(
            "references/gloritools/genome_AG/",
            "chrLength.txt",
            "chrName.txt",
            "Genome",
            "Log.out",
            "SAindex",
            "chrNameLength.txt",
            "chrStart.txt",
            "genomeParameters.txt",
            "SA",
        ),
    output:
        unmapped_fastq="results/{sample}/gloritools/{sample}_unmapped.fq",
    params:
        output_prefix=lambda w, output: output.unmapped_fastq.replace(
            "_unmapped.fq", ""
        ),
    threads: config["threads"]["gloritools_star_mapping"]
    conda:
        "../envs/star.yaml"
    log:
        "logs/gloritools/{sample}_gloritools_genome_mapping.log",
        "logs/gloritools/{sample}_gloritools_genome_mapping.err",
    benchmark:
        "benchmarks/gloritools/{sample}_gloritools_genome_mapping.txt"
    script:
        "../scripts/gloritools_star_mapping.py"
