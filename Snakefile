import os
import pandas as pd

samples_df = pd.read_csv("samples.tab", sep="\t", comment="#")
SAMPLES = samples_df["sample_id"].tolist()

species = "homo_sapiens"
genome_version = None

if os.path.exists("splicekit.config"):
    with open("splicekit.config") as config_file:
        for line in config_file:
            exec(line.strip())

genome_version = f"--genome_version {genome_version}" if genome_version else ""

rule all:
    input:
        expand("bam/{sample}.bam", sample=SAMPLES),
        expand("data/sample_junctions_data/sample_{sample}_raw.tab.gz", sample=SAMPLES)

rule map_fastq:
    input:
        fastq1="fastq/{sample}_1.fastq.gz",
        fastq2="fastq/{sample}_2.fastq.gz"
    output:
        "bam/{sample}.bam"
    shell:
        f"""
        echo mapping {{input.fastq1}} {{input.fastq2}} {{output}}
        echo species = {species}
        echo genome_version = {genome_version}
        echo pybio star {species} {{input.fastq1}} {{input.fastq2}} {{output}} -t 7 {genome_version}
        pybio star {species} {{input.fastq1}} {{input.fastq2}} {{output}} -t 7 {genome_version}
        """

rule junctions:
    container:
        "docker://ghcr.io/bedapub/splicekit:main"
    input:
        "bam/{sample}.bam"
    output:
        fname_comp="data/sample_junctions_data/sample_{sample}_raw.tab.gz",
    params:
        fname=lambda wildcards: f"data/sample_junctions_data/sample_{wildcards.sample}"
    shell:
        """
        python /usr/splicekit/splicekit/core/junctions.py {input} {params.fname}
        """
