# define the config file
configfile: "/capstor/scratch/cscs/ymorsy/workflow/config.yaml"

import os

def get_ids(filename):
    with open(filename, "r") as file:
        id_column = []
        for line in file:
            columns = line.rstrip("\n").split(",")
            id_column.append(columns[0])
    return id_column

ids = get_ids(config["ids_file_name"])

rule targets:
    input:
        expand(config["output_path"] + "/01_non_human/{sample}_unmapped.fastq.1.gz", sample=ids),
        expand(config["output_path"] + "/01_non_human/{sample}_unmapped.fastq.2.gz", sample=ids),
        expand(config["output_path"] + "/02_non_human_concat/{sample}_unmapped.fastq.gz", sample=ids),
        expand(config["output_path"] + "/03_Metaphlan_output/{sample}_profile.tsv", sample=ids),
        expand(config["output_path"] + "/04_Humann_output/{sample}", sample=ids),
        expand(config["output_path"] + "/kraken2/k2_{sample}_report.txt", sample=ids),
        expand(config["output_path"] + "/kraken2/k2_nonhuman_{sample}_report.txt", sample=ids),


rule host_remove:
    input:
        R1=config["input_fastq_path"] + "/{sample}" + config["R1_extension"],
        R2=config["input_fastq_path"] + "/{sample}" + config["R2_extension"],
    params:
        prefix=config["output_path"] + "/01_non_human/{sample}_unmapped.fastq.gz",
        db_path_bowtie2=config["db_path"] + "/01_bowtie2_human/human_index",
    output:
        Out_R1=config["output_path"] + "/01_non_human/{sample}_unmapped.fastq.1.gz",
        Out_R2=config["output_path"] + "/01_non_human/{sample}_unmapped.fastq.2.gz",
    threads: 20,
    shell:
        """
        bowtie2 \
            -x {params.db_path_bowtie2} \
            -p {threads} \
            -1 {input.R1} \
            -2 {input.R2}\
            --un-conc-gz {params.prefix}
        """

rule concatenate:
    input:
        R1=rules.host_remove.output.Out_R1,
        R2=rules.host_remove.output.Out_R2,
    output:
        fastq=config["output_path"] + "/02_non_human_concat/{sample}_unmapped.fastq.gz",
    shell:
        """
        cat {input.R1} {input.R2} > {output.fastq}
        """

rule metaphlan:
    input:
        concat_fastq=rules.concatenate.output.fastq,
    output:
        tsv=config["output_path"] + "/03_Metaphlan_output/{sample}_profile.tsv",
    params:
        db_path_bowtie2=config["db_path"] + "/02_metaphlan3_db",
    threads: 20,
    shell:
        """
        metaphlan \
            {input.concat_fastq} \
            --input_type fastq \
            --bowtie2db {params.db_path_bowtie2} \
            --ignore_eukaryotes \
            --ignore_archaea \
            --nproc {threads} \
            -o {output.tsv}
        """

rule humann:
    input:
        concat_fastq=rules.concatenate.output.fastq,
        tsv=rules.metaphlan.output.tsv,
    output:
        output=config["output_path"] +"/04_Humann_output/{sample}",
    params:
        prefix="{sample}",
        db_path_choco=config["db_path"] + "/04_chocophlan",
        db_path_uniref=config["db_path"] + "/05_uniref",
        db_path_bowtie2=config["db_path"] + "/02_metaphlan3_db",
    threads: 20,
    shell:
        """
        humann \
            --input {input.fastq} \
            --taxonomic-profile {input.tsv} \
            --output {output.output} \
            --output-basename {params.prefix} \
            --threads {threads} \
            --memory-use maximum \
            --prescreen-threshold 0.0001 \
            --nucleotide-database {params.db_path_choco} \
            --protein-database {params.db_path_uniref} \
            --metaphlan-options "--bowtie2db {params.db_path_bowtie2} --nproc {threads}"
        """

rule kraken2:
    input:
        R1=config["input_fastq_path"] + "/{sample}" + config["R1_extension"],
        R2=config["input_fastq_path"] + "/{sample}" + config["R2_extension"],
    output:
        res=config["output_path"] + "/kraken2/k2_{sample}_report.txt",
    params:
        db_path_kra=config["db_path"] + "/kradb25",
    resources:
        threads=20,
        mem_mb=7700,
    shell:
        """
        kraken2 \
        --db {params.db_path_kra} \
        --threads {resources.threads} \
        --paired \
        --use-names \
        --gzip-compressed \
        --report-zero-counts \
        --use-mpa-style \
        --report {output.res} \
        {input.R1} \
        {input.R2}
        """

rule kraken2_nonhuman:
    input:
        R1=rules.host_remove.output.Out_R1,
        R2=rules.host_remove.output.Out_R2,
    output:
        res_nonhuman=config["output_path"] + "/kraken2/k2_nonhuman_{sample}_report.txt",
    params:
        db_path_kra=config["db_path"] + "/kradb25",
    resources:
        threads=20,
        mem_mb=7700,
    shell:
        """
        kraken2 \
        --db {params.db_path_kra} \
        --threads {resources.threads} \
        --paired \
        --use-names \
        --gzip-compressed \
        --report-zero-counts \
        --use-mpa-style \
        --report {output.res} \
        {input.R1} \
        {input.R2}
        """