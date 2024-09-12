rule glori_preptbl:
    input:
        gtf=config["reference"]["gtf"],
    output:
        annotation_tbl=temp("references/gloritools/annotation.tbl"),
    threads: 1
    container:
        "docker://btrspg/gloritools:latest"
    log:
        log="logs/gloritools/indexing_gtf2anno.log",
        err="logs/gloritools/indexing_gtf2anno.err",
    benchmark:
        "benchmarks/gloritools/indexing_gtf2anno.txt"
    shell:
        " python3 /opt/GLORI-tools/get_anno/gtf2anno.py "
        " -i {input.gtf} "
        " -o {output.annotation_tbl} 1>{log.log} 2>{log.err} "


rule glori_preptbl2:
    input:
        annotation_tbl="references/gloritools/annotation.tbl",
    output:
        annotation_tbl2="references/gloritools/annotation2.tbl",
    threads: 1
    container:
        "docker://btrspg/gloritools:latest"
    log:
        log="logs/gloritools/indexing_preptbl2.log",
    benchmark:
        "benchmarks/gloritools/indexing_gtf2anno2.txt"
    shell:
        """awk '$3!~/_/&&$3!="na"' {input.annotation_tbl} | grep -v 'unknown_transcript'  > {output.annotation_tbl2} && echo date > {log.log}"""


rule glori_select_longest_transcripts:
    input:
        annotation_tbl="references/gloritools/annotation2.tbl",
        transcriptome_reference=config["reference"]["transcriptome_fa"],
    output:
        select_transcriptome_reference="references/gloritools/selected_transcriptome.fa",
    threads: 1
    container:
        "docker://btrspg/gloritools:latest"
    log:
        log="logs/gloritools/indexing_select_longest_transcripts.log",
        err="logs/gloritools/indexing_select_longest_transcripts.err",
    benchmark:
        "benchmarks/gloritools/indexing_select_longest_transcripts.txt"
    shell:
        " python3 /opt/GLORI-tools/get_anno/selected_longest_transcrpts_fa.py "
        " -anno {input.annotation_tbl} "
        " -fafile {input.transcriptome_reference} "
        " --outname_prx {output.select_transcriptome_reference} 1>{log.log} 2>{log.err} "


rule glori_build_index_transcriptome:
    input:
        select_transcriptome_reference="references/gloritools/selected_transcriptome.fa",
    output:
        transcript_fa="references/gloritools/transcriptome_index/transcriptome.AG_conversion.fa",
        transcript_indexes=multiext(
            "references/gloritools/transcriptome_index/transcriptome.AG_conversion.fa",
            ".1.ebwt",
            ".2.ebwt",
            ".3.ebwt",
            ".4.ebwt",
            ".rev.1.ebwt",
            ".rev.2.ebwt",
        ),
        outdir=directory("references/gloritools/transcriptome_index/"),
    threads: config["threads"]["gloritools_build_index"]
    params:
        prefix="transcriptome",
    container:
        "docker://btrspg/gloritools:latest"
    log:
        log="logs/gloritools/indexing_build_index_transcriptome.log",
        err="logs/gloritools/indexing_build_index_transcriptome.err",
    benchmark:
        "benchmarks/gloritools/indexing_build_index_transcriptome.txt"
    shell:
        " python3 /opt/GLORI-tools/pipelines/build_transcriptome_index.py "
        " -f {input.select_transcriptome_reference} "
        " -p {threads} "
        " -o {output.outdir} "
        " -pre {params.prefix} 1>{log.log} 2>{log.err} "


rule glori_build_index_genome:
    input:
        genome_reference=config["reference"]["genome_fa"],
    output:
        genome_indexes=multiext(
            "references/gloritools/genome_index/genome",
            ".AG_conversion.fa",
            ".AG_conversion.fa.fai",
            ".rvsCom.fa",
            ".rvsCom.fa.fai",
        ),
        outdir=directory("references/gloritools/genome_index/"),
        raw_index=directory("references/gloritools/genome_index/genome.AG_conversion"),
        rvs_index=directory("references/gloritools/genome_index/genome.rvsCom"),
    threads: config["threads"]["gloritools_build_index"]
    container:
        "docker://btrspg/gloritools:latest"
    log:
        log="logs/gloritools/indexing_build_index_genome.log",
        err="logs/gloritools/indexing_build_index_genome.err",
    benchmark:
        "benchmarks/gloritools/indexing_build_index_genome.txt"
    shell:
        " python3 /opt/GLORI-tools/pipelines/build_genome_index.py "
        " -f {input.genome_reference} "
        " -p {threads} "
        " -o {output.outdir}/ "
        " -pre genome 1>{log.log} 2>{log.err} "


rule glori_index_baseanno:
    input:
        annotation_tbl="references/gloritools/annotation2.tbl",
    output:
        baseanno_tbl="references/gloritools/baseanno.tbl",
    threads: 1
    container:
        "docker://btrspg/gloritools:latest"
    log:
        log="logs/gloritools/indexing_index_baseanno.log",
        err="logs/gloritools/indexing_index_baseanno.err",
    benchmark:
        "benchmarks/gloritools/indexing_index_baseanno.txt"
    shell:
        " python3 /opt/GLORI-tools/get_anno/anno_to_base.py "
        " -i {input.annotation_tbl} "
        " -o {output.baseanno_tbl} 1>{log.log} 2>{log.err} "


rule glori_gtf2genelist:
    input:
        gtf=config["reference"]["gtf"],
        annotation_tbl="references/gloritools/annotation2.tbl",
        transcriptome_reference=config["reference"]["transcriptome_fa"],
    output:
        gene_list=temp("references/gloritools/genelist.tmp.tbl"),
    threads: 1
    container:
        "docker://btrspg/gloritools:latest"
    log:
        log="logs/gloritools/indexing_gtf2genelist.log",
        err="logs/gloritools/indexing_gtf2genelist.err",
    benchmark:
        "benchmarks/gloritools/indexing_gtf2genelist.txt"
    shell:
        " python3 /opt/GLORI-tools/get_anno/gtf2genelist.py "
        " -i {input.annotation_tbl} "
        " -f {input.transcriptome_reference} "
        " -o {output.gene_list} 1>{log.log} 2>{log.err} "


rule glori_filtered_genelist:
    input:
        gene_list="references/gloritools/genelist.tmp.tbl",
    output:
        gene_list2="references/gloritools/genelist.tbl",
    threads: 1
    container:
        "docker://btrspg/gloritools:latest"
    log:
        log="logs/gloritools/indexing_filtered_genelist.log",
    benchmark:
        "benchmarks/gloritools/indexing_filtered_genelist.txt"
    shell:
        """awk '$6!~/_/&&$6!="na"' {input.genelist} > {output.genelist2} && echo date > {log.log}"""


rule glori_remove_redundance:
    input:
        baseanno=config["reference"]["gtf"],
        genelist="references/gloritools/genelist.tbl",
    output:
        noredundance_base="references/gloritools/noredundance.base",
    threads: 1
    container:
        "docker://btrspg/gloritools:latest"
    log:
        log="logs/gloritools/indexing_remove_redundance.log",
        err="logs/gloritools/indexing_remove_redundance.err",
    benchmark:
        "benchmarks/gloritools/indexing_remove_redundance.txt"
    shell:
        " python3 /opt/GLORI-tools/get_anno/anno_to_base_remove_redundance_v1.0.py "
        " -i {input.baseanno} "
        " -g {input.genelist} "
        " -o {output.noredundance_base} 1>{log.log} 2>{log.err} "
