rule gloritools_genome_mapping_treated:
    input:
        fastq="results/{sample}/gloritools/cleandata/{sample}_rmdup.fq.gz",
        ag_genome_indexes=multiext(
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
        ag_transcriptome_indexes=multiext(
            "references/gloritools/transcriptome_AG.fa",
            ".1.ebwt",
            ".2.ebwt",
            ".3.ebwt",
            ".4.ebwt",
            ".rev.1.ebwt",
            ".rev.2.ebwt",
        ),
        ag_transcriptome_reference="references/gloritools/transcriptome_AG.fa",
    output:
        unmapped_fastq="results/{sample}/gloritools/{sample}_1_unmapped.fq",
        unmapped_fastq2="results/{sample}/gloritools/{sample}_2_unmapped.fq",
    params:
        output_prefix=lambda w, output: output.unmapped_fastq.replace(
            "_1_unmapped.fq", ""
        ),
    threads: config["threads"]["gloritools_star_mapping"]
    conda:
        "../envs/star.yaml"
    log:
        log="logs/gloritools/{sample}_gloritools_genome_mapping.log",
        err="logs/gloritools/{sample}_gloritools_genome_mapping.err",
    benchmark:
        "benchmarks/gloritools/{sample}_gloritools_genome_mapping.txt"
    script:
        "../scripts/gloritools_treated_mapping.py"
