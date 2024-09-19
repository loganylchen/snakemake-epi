rule gloritools_star_mapping_untreated:
    input:
        fastq="results/{sample}/gloritools/cleandata/{sample}_rmdup.fq.gz",
        genome_indexes=multiext(
            "references/gloritools/genome/",
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
        genome_unmapped_fastq="results/{sample}/gloritools/untreated/{sample}_genome_unmapped.fq",
        genome_star_bam="results/{sample}/gloritools/untreated/{sample}.star.bam",
    params:
        output_prefix=lambda w, output: output.genome_unmapped_fastq.replace(
            "_genome_unmapped.fq", ""
        ),
    threads: config["threads"]["gloritools_star_mapping"]
    conda:
        "../envs/mapping.yaml"
    log:
        "logs/gloritools/{sample}_untreated_gloritools_star_mapping.log",
    benchmark:
        "benchmarks/gloritools/{sample}_untreated_gloritools_star_mapping.txt"
    script:
        "../scripts/gloritools_untreated_star_mapping.py"


rule gloritools_bowtie_mapping_untreated:
    input:
        fastq="results/{sample}/gloritools/untreated/{sample}_genome_unmapped.fq",
        transcriptome_indexes=multiext(
            "references/gloritools/selected_transcriptome.fa",
            ".1.ebwt",
            ".2.ebwt",
            ".3.ebwt",
            ".4.ebwt",
            ".rev.1.ebwt",
            ".rev.2.ebwt",
        ),
        transcriptome_reference="references/gloritools/selected_transcriptome.fa",
    output:
        transcriptome_unmapped_fastq="results/{sample}/gloritools/untreated/{sample}_transcriptome_unmapped.fq",
        transcriptome_bowtie_bam="results/{sample}/gloritools/untreated/{sample}.bowtie.bam",
    params:
        output_prefix=lambda w, output: output.transcriptome_unmapped_fastq.replace(
            "_transcriptome_unmapped.fq", ""
        ),
    threads: config["threads"]["gloritools_star_mapping"]
    conda:
        "../envs/mapping.yaml"
    log:
        "logs/gloritools/{sample}_untreated_gloritools_bowtie_mapping.log",
    benchmark:
        "benchmarks/gloritools/{sample}_untreated_gloritools_bowtie_mapping.txt"
    script:
        "../scripts/gloritools_untreated_bowtie_mapping.py"


rule gloritools_star_ag_mapping_treated:
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
    output:
        ag_genome_unmapped_fastq="results/{sample}/gloritools/treated/{sample}_star_ag_unmapped.fq",
        ag_change_fastq="results/{sample}/gloritools/treated/{sample}_AG.fq.gz",
        info_json="results/{sample}/gloritools/treated/{sample}_AG_changed_info.json",
        ag_genome_star_bam="results/{sample}/gloritools/treated/{sample}.star.ag.bam",
    params:
        output_prefix=lambda w, output: output.ag_genome_unmapped_fastq.replace(
            "_star_ag_unmapped.fq", ""
        ),
    threads: config["threads"]["gloritools_star_mapping"]
    conda:
        "../envs/mapping.yaml"
    log:
        "logs/gloritools/{sample}_treated_gloritools_star_ag_mapping.log",
    benchmark:
        "benchmarks/gloritools/{sample}_treated_gloritools_star_ag_mapping.txt"
    script:
        "../scripts/gloritools_treated_star_ag_mapping.py"


rule gloritools_star_rvs_mapping_treated:
    input:
        fastq="results/{sample}/gloritools/treated/{sample}_star_ag_unmapped.fq",
        info_json="results/{sample}/gloritools/treated/{sample}_AG_changed_info.json",
        rvs_genome_indexes=multiext(
            "references/gloritools/genome_rc_AG/",
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
        rvs_genome_unmapped_fastq="results/{sample}/gloritools/treated/{sample}_star_rvs_unmapped.fq",
        rvs_genome_star_bam="results/{sample}/gloritools/treated/{sample}.star.rvs.bam",
    params:
        output_prefix=lambda w, output: output.rvs_genome_unmapped_fastq.replace(
            "_star_rvs_unmapped.fq", ""
        ),
    threads: config["threads"]["gloritools_star_mapping"]
    conda:
        "../envs/mapping.yaml"
    log:
        "logs/gloritools/{sample}_treated_gloritools_star_rvs_mapping.log",
    benchmark:
        "benchmarks/gloritools/{sample}_treated_gloritools_star_rvs_mapping.txt"
    script:
        "../scripts/gloritools_treated_star_rvs_mapping.py"


rule gloritools_bowtie_ag_mapping_treated:
    input:
        fastq="results/{sample}/gloritools/treated/{sample}_star_rvs_unmapped.fq",
        info_json="results/{sample}/gloritools/treated/{sample}_AG_changed_info.json",
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
        info_json="results/{sample}/gloritools/treated/{sample}_AG_changed_info.json",
    output:
        ag_transcriptome_fastq="results/{sample}/gloritools/treated/{sample}_bowtie_ag_unmapped.fq",
        ag_transcriptome_bowtie_bam="results/{sample}/gloritools/treated/{sample}.bowtie.ag.bam",
    params:
        output_prefix=lambda w, output: output.ag_transcriptome_fastq.replace(
            "_bowtie_ag_unmapped.fq", ""
        ),
    threads: config["threads"]["gloritools_star_mapping"]
    conda:
        "../envs/mapping.yaml"
    log:
        "logs/gloritools/{sample}_treated_gloritools_bowtie_ag_mapping.log",
    benchmark:
        "benchmarks/gloritools/{sample}_treated_gloritools_bowtie_ag_mapping.txt"
    script:
        "../scripts/gloritools_treated_bowtie_ag_mapping.py"
