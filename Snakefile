import os
import pandas as pd

container: "docker://ghcr.io/bedapub/splicekit:main"
threads = config.get("threads", 1)

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
        # bam files
        expand("bam/{sample}.bam", sample=SAMPLES),
        # sample junction counts from bam files
        expand("data/sample_junctions_data/sample_{sample}_raw.tab.gz", sample=SAMPLES),
        # GTF files
        "reference/junctions.tab.gz",
        "reference/acceptor_anchors.gtf.gz",
        "reference/donor_anchors.gtf.gz",
        # anchors
        expand("data/sample_acceptor_anchors_data/sample_{sample}.tab.gz", sample=SAMPLES),
        expand("logs/count_acceptor_anchors/sample_{sample}.tab.summary", sample=SAMPLES),
        expand("data/sample_donor_anchors_data/sample_{sample}.tab.gz", sample=SAMPLES),
        expand("logs/count_donor_anchors/sample_{sample}.tab.summary", sample=SAMPLES)

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
        echo pybio star {species} {{input.fastq1}} {{input.fastq2}} {{output}} -t {threads} {genome_version}
        pybio star {species} {{input.fastq1}} {{input.fastq2}} {{output}} -t {threads} {genome_version}
        """

rule junctions:
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

rule junctions_make_master:
    input:
        expand("data/sample_junctions_data/sample_{sample}_raw.tab.gz", sample=SAMPLES)
    output:
        "reference/junctions.tab.gz"
    shell:
        """
        python -c 'import splicekit; splicekit.core.junctions.make_master();'
        """

rule anchors_gtf:
    input:
        "reference/junctions.tab.gz",
        expand("data/sample_junctions_data/sample_{sample}_raw.tab.gz", sample=SAMPLES)
    output:
        "reference/acceptor_anchors.gtf.gz",
        "reference/donor_anchors.gtf.gz"
    shell:
        """
        python -c 'import splicekit; splicekit.core.anchors.write_anchor_gtf();'
        """

rule anchors:
    params:
        library_type_insert = {"single-end":"", "paired-end":"-p "}[library_type],
        library_strand_insert = {"FIRST_READ_TRANSCRIPTION_STRAND":1, "SINGLE_STRAND":1, "SINGLE_REVERSE":1, "SECOND_READ_TRANSCRIPTION_STRAND":2, "NONE":0}[library_strand],
        tab_fname = lambda wildcards: f"data/sample_{wildcards.anchor_type}_anchors_data/sample_{wildcards.sample}.tab",
        logs_dir = lambda wildcards: f'logs/count_{wildcards.anchor_type}_anchors'
    input:
        gtf_fname = "reference/{anchor_type}_anchors.gtf.gz",
        bam_fname = "bam/{sample}.bam"
    output:
        "data/sample_{anchor_type}_anchors_data/sample_{sample}.tab.gz"
    wildcard_constraints:
        anchor_type="acceptor|donor"
    shell:
        """
        featureCounts {params.library_type_insert}-s {params.library_strand_insert} -M -O -T {threads} -F GTF -f -t anchor -g {wildcards.anchor_type}_anchor_id -a {input.gtf_fname} -o {params.tab_fname} {input.bam_fname}
        cp {params.tab_fname} {params.tab_fname}_temp
        echo "anchor_id\tcount" >| {params.tab_fname}
        tail -n +3 {params.tab_fname}_temp| cut -f1,7 >> {params.tab_fname}
        rm {params.tab_fname}_temp
        mv {params.tab_fname}.summary {params.logs_dir}
        gzip -f {params.tab_fname}
        """
