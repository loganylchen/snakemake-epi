# rule run_gloritools:
#     input:
#         rmumi_fastq="results/{sample}/gloritools/cleandata/{sample}_rmumi.fq",
#         genome_ag_reference="references/gloritools/genome_index/genome.AG_conversion.fa",
#         genome_fa=config["reference"]["genome_fa"],
#         genome_rev_reference="references/gloritools/genome_index/genome.rvsCom.fa",
#         transcriptome_ag_reference="references/gloritools/transcriptome_index/transcriptome.AG_conversion.fa",
#         annotation_tbl2="references/gloritools/annotation2.tbl",
#         baseanno_tbl="references/gloritools/baseanno.tbl",
#     output:
#         m6a_results="results/{sample}/gloritools/glori/{sample}.totalm6A.FDR.csv",
#     threads: config["threads"]["gloritools_run"]
#     container:
#         "docker://btrspg/gloritools:latest"
#     params:
#         params=config["gloritools"]["gloritools_run"],
#         outdir=lambda x, output: os.path.dirname(output.m6a_results),
#         prefix="{sample}",
#     log:
#         log="logs/gloritools/{sample}_gloritools_m6Arun.log",
#         err="logs/gloritools/{sample}_gloritools_m6Arun.err",
#     benchmark:
#         "benchmarks/gloritools/{sample}_gloritools_m6Arun.txt"
#     shell:
#         " python3 /opt/GLORI-tools/run_GLORI.py "
#         " -i /opt/GLORI-tools "
#         " -q {input.rmumi_fastq} "
#         " -T {threads} -f {input.genome_ag_reference} "
#         " -f2 {input.genome_fa} "
#         " -rvs {input.genome_rev_reference} "
#         " -Tf {input.transcriptome_ag_reference} "
#         " -a {input.annotation_tbl2} "
#         " -b {input.baseanno_tbl} "
#         " -pre {params.prefix} -o {params.outdir} "
#         " {params.params} 1>{log.log} 2>{log.err} "
# rule run_gloritools_control:
#     input:
#         rmumi_fastq="results/{sample}/gloritools/cleandata/{sample}_rmumi.fq",
#         genome_ag_reference="references/gloritools/genome_index/genome.AG_conversion.fa",
#         genome_fa=config["reference"]["genome_fa"],
#         genome_rev_reference="references/gloritools/genome_index/genome.rvsCom.fa",
#         transcriptome_ag_reference="references/gloritools/transcriptome_index/transcriptome.AG_conversion.fa",
#         annotation_tbl2="references/gloritools/annotation2.tbl",
#         baseanno_tbl="references/gloritools/baseanno.tbl",
#     output:
#         m6a_results="results/{sample}/gloritools/glori_asControl/{sample}.totalm6A.FDR.csv",
#     threads: config["threads"]["gloritools_run"]
#     container:
#         "docker://btrspg/gloritools:latest"
#     params:
#         params=config["gloritools"]["gloritools_run_control"],
#         outdir=lambda x, output: os.path.dirname(output.m6a_results),
#         prefix="{sample}",
#     log:
#         log="logs/gloritools/{sample}_gloritools_control_m6Arun.log",
#         err="logs/gloritools/{sample}_gloritools_control_m6Arun.err",
#     benchmark:
#         "benchmarks/gloritools/{sample}_gloritools_control_m6Arun.txt"
#     shell:
#         " python3 /opt/GLORI-tools/run_GLORI.py "
#         " -i /opt/GLORI-tools "
#         " -q {input.rmumi_fastq} "
#         " -T {threads} "
#         " -f {input.genome_ag_reference} "
#         " -f2 {input.genome_fa} "
#         " -rvs {input.genome_rev_reference} "
#         " -Tf {input.transcriptome_ag_reference} "
#         " -a {input.annotation_tbl2} "
#         " -b {input.baseanno_tbl} "
#         " -pre {params.prefix} -o {params.outdir} "
#         " {params.params} 1>{log.log} 2>{log.err} "
