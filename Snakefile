import os
import pandas as pd
import splicekit

DEFAULT_CORES = 1
DEFAULT_MEM = 4 # MB
DEFAULT_TIME = "01:00:00" # 1h

if not os.path.exists("results"):
    splicekit.setup()

if not os.path.exists("annotation/comparisons.tab"):
    splicekit.annotation()

splicekit_folder = os.path.dirname(splicekit.__file__)

container: "docker://ghcr.io/bedapub/splicekit:main"

samples_df = pd.read_csv("samples.tab", sep="\t", comment="#")
SAMPLES = samples_df["sample_id"].tolist()

species = splicekit.config.species
genome_version = splicekit.config.genome_version

if os.path.exists("splicekit.config"):
    with open("splicekit.config") as config_file:
        for line in config_file:
            exec(line.strip())

genome_version = f"--genome_version {genome_version}" if genome_version else ""

if not os.path.exists("annotation/comparisons.tab"):
    splicekit.annotation()

comparisons_df = pd.read_csv("annotation/comparisons.tab", sep="\t", comment="#")
COMPARISONS = comparisons_df["comparison"].tolist()
comparisons_df.set_index('comparison', inplace=True)
comps = comparisons_df.to_dict(orient='index')

rule all:
    input:
        # annotation
        "annotation/comparisons.tab",

        # bam files
        expand("bam/{sample}.bam", sample=SAMPLES),

        # sample junction counts from bam files
        expand("data/sample_junctions_data/sample_{sample}_raw.tab.gz", sample=SAMPLES),

        # GTF files
        "reference/junctions.tab.gz",
        "reference/acceptor_anchors.gtf.gz",
        "reference/donor_anchors.gtf.gz",
        "reference/exons.gtf.gz",
        "reference/genes.gtf.gz",

        # anchors
        expand("data/sample_acceptor_anchors_data/sample_{sample}.tab.gz", sample=SAMPLES),
        expand("data/sample_donor_anchors_data/sample_{sample}.tab.gz", sample=SAMPLES),
        # anchor counts
        "data/samples_acceptor_anchors_counts.tab.gz",
        "data/samples_donor_anchors_counts.tab.gz",

        # exons counts
        expand("data/sample_exons_data/sample_{sample}.tab.gz", sample=SAMPLES),
        "data/samples_exons_counts.tab.gz",

        # junctions counts
        expand("data/sample_junctions_data/sample_{sample}.tab.gz", sample=SAMPLES),
        "data/samples_junctions_counts.tab.gz",

        # genes counts
        expand("data/sample_genes_data/sample_{sample}.tab.gz", sample=SAMPLES),
        "data/samples_genes_counts.tab.gz",

        expand("results/edgeR/{feature_type}/{comparison}_altsplice.tab.gz", feature_type=["junctions", "exons", "donor_anchors", "acceptor_anchors"], comparison=COMPARISONS),
        expand("results/edgeR/{feature_type}/{comparison}_difffeature.tab.gz", feature_type=["junctions", "exons", "donor_anchors", "acceptor_anchors", "genes"], comparison=COMPARISONS),

        # juDGE
        "results/judge/scored.tab.gz",

        # juan
        "results/edgeR/juan.done"

rule setup:
    input:
        "samples.tab"
    output:
        "annotation/comparisons.tab"
    shell:
        """
        splicekit setup
        splicekit annotation
        """

rule map_fastq:
    resources:
        mem = 40, # GB
        time = "02:00:00",
        cores = 8
    input:
        fastq1="fastq/{sample}_1.fastq.gz",
        fastq2="fastq/{sample}_2.fastq.gz"
    output:
        "bam/{sample}.bam",
    shell:
        f"""
        echo mapping {{input.fastq1}} {{input.fastq2}} {{output}}
        echo species = {species}
        echo genome_version = {genome_version}
        echo pybio star {species} {{input.fastq1}} {{input.fastq2}} {{output}} -t {{resources.cores}} {genome_version}
        pybio star {species} {{input.fastq1}} {{input.fastq2}} {{output}} -t {{resources.cores}} {genome_version}
        """

rule junctions:
    input:
        "bam/{sample}.bam"
    output:
        fname_comp="data/sample_junctions_data/sample_{sample}_raw.tab.gz",
    params:
        fname=lambda wildcards: f"data/sample_junctions_data/sample_{wildcards.sample}",
        junctions_path = os.path.join(splicekit_folder, "core", "junctions.py")
    resources:
        mem = DEFAULT_MEM,
        time = DEFAULT_TIME,
        cores = DEFAULT_CORES
    shell:
        """
        python {params.junctions_path} {input} {params.fname}
        """

rule junctions_per_sample:
    input:
        expand("data/sample_junctions_data/sample_{sample}_raw.tab.gz", sample=SAMPLES),
        "reference/junctions.tab.gz"
    output:
        expand("data/sample_junctions_data/sample_{sample}.tab.gz", sample=SAMPLES),
    resources:
        mem = DEFAULT_MEM,
        time = DEFAULT_TIME,
        cores = DEFAULT_CORES
    run:
        import splicekit
        splicekit.core.annotation.make_comparisons()
        splicekit.core.junctions.junctions_per_sample()

rule junctions_make_master:
    input:
        expand("data/sample_junctions_data/sample_{sample}_raw.tab.gz", sample=SAMPLES)
    output:
        "reference/junctions.tab.gz"
    resources:
        mem = DEFAULT_MEM,
        time = DEFAULT_TIME,
        cores = DEFAULT_CORES
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
    resources:
        mem = DEFAULT_MEM,
        time = DEFAULT_TIME,
        cores = DEFAULT_CORES
    shell:
        """
        python -c 'import splicekit; splicekit.core.anchors.write_anchor_gtf();'
        """

rule exons_gtf:
    output:
        "reference/exons.gtf.gz",
    resources:
        mem = DEFAULT_MEM,
        time = DEFAULT_TIME,
        cores = DEFAULT_CORES
    shell:
        """
        python -c 'import splicekit; splicekit.core.exons.write_exons_gtf()'
        """

rule genes_gtf:
    output:
        "reference/genes.gtf.gz",
    resources:
        mem = DEFAULT_MEM,
        time = DEFAULT_TIME,
        cores = DEFAULT_CORES
    shell:
        """
        python -c 'import splicekit; splicekit.core.genes.write_genes_gtf()'
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
    resources:
        mem = DEFAULT_MEM,
        time = DEFAULT_TIME,
        cores = DEFAULT_CORES
    shell:
        """
        featureCounts {params.library_type_insert}-s {params.library_strand_insert} -M -O -T {resources.cores} -F GTF -f -t anchor -g {wildcards.anchor_type}_anchor_id -a {input.gtf_fname} -o {params.tab_fname} {input.bam_fname}
        cp {params.tab_fname} {params.tab_fname}_temp
        echo "anchor_id\tcount" >| {params.tab_fname}
        tail -n +3 {params.tab_fname}_temp| cut -f1,7 >> {params.tab_fname}
        rm {params.tab_fname}_temp
        mv {params.tab_fname}.summary {params.logs_dir}
        gzip -f {params.tab_fname}
        """

rule exons:
    params:
        library_type_insert = {"single-end":"", "paired-end":"-p "}[library_type],
        library_strand_insert = {"FIRST_READ_TRANSCRIPTION_STRAND":1, "SINGLE_STRAND":1, "SINGLE_REVERSE":1, "SECOND_READ_TRANSCRIPTION_STRAND":2, "NONE":0}[library_strand],
        tab_fname = lambda wildcards: f"data/sample_exons_data/sample_{wildcards.sample}.tab",
        logs_dir = lambda wildcards: f'logs/count_exons'
    input:
        gtf_fname = "reference/exons.gtf.gz",
        bam_fname = "bam/{sample}.bam"
    output:
        "data/sample_exons_data/sample_{sample}.tab.gz"
    resources:
        mem = DEFAULT_MEM,
        time = DEFAULT_TIME,
        cores = DEFAULT_CORES
    shell:
        """
        featureCounts {params.library_type_insert}-s {params.library_strand_insert} -M -O -T {resources.cores} -F GTF -f -t exon -g exon_id -a {input.gtf_fname} -o {params.tab_fname} {input.bam_fname}
        cp {params.tab_fname} {params.tab_fname}_temp
        echo "anchor_id\tcount" >| {params.tab_fname}
        tail -n +3 {params.tab_fname}_temp| cut -f1,7 >> {params.tab_fname}
        rm {params.tab_fname}_temp
        mv {params.tab_fname}.summary {params.logs_dir}
        gzip -f {params.tab_fname}
        """

rule genes:
    params:
        library_type_insert = {"single-end":"", "paired-end":"-p "}[library_type],
        library_strand_insert = {"FIRST_READ_TRANSCRIPTION_STRAND":1, "SINGLE_STRAND":1, "SINGLE_REVERSE":1, "SECOND_READ_TRANSCRIPTION_STRAND":2, "NONE":0}[library_strand],
        tab_fname = lambda wildcards: f"data/sample_genes_data/sample_{wildcards.sample}.tab",
        logs_dir = lambda wildcards: f'logs/count_genes'
    input:
        gtf_fname = "reference/genes.gtf.gz",
        bam_fname = "bam/{sample}.bam"
    output:
        "data/sample_genes_data/sample_{sample}.tab.gz"
    resources:
        mem = DEFAULT_MEM,
        time = DEFAULT_TIME,
        cores = DEFAULT_CORES
    shell:
        """
        featureCounts {params.library_type_insert}-s {params.library_strand_insert} -M -O -T {resources.cores} -F GTF -f -t exon -g gene_id -a {input.gtf_fname} -o {params.tab_fname} {input.bam_fname}
        cp {params.tab_fname} {params.tab_fname}_temp
        echo "anchor_id\tcount" >| {params.tab_fname}
        tail -n +3 {params.tab_fname}_temp| cut -f1,7 >> {params.tab_fname}
        rm {params.tab_fname}_temp
        mv {params.tab_fname}.summary {params.logs_dir}
        gzip -f {params.tab_fname}
        """

rule anchors_counts:
    input:
        expand("data/sample_donor_anchors_data/sample_{sample}.tab.gz", sample=SAMPLES),
        expand("data/sample_acceptor_anchors_data/sample_{sample}.tab.gz", sample=SAMPLES),
        "reference/junctions.tab.gz"
    output:
        "data/samples_acceptor_anchors_counts.tab.gz",
        "data/samples_donor_anchors_counts.tab.gz"
    resources:
        mem = DEFAULT_MEM,
        time = DEFAULT_TIME,
        cores = DEFAULT_CORES
    run:
        import splicekit
        splicekit.core.annotation.make_comparisons()
        splicekit.core.features.load_genes()
        splicekit.core.features.read_anchors("donor")
        splicekit.core.features.read_anchors("acceptor")
        os.system("rm -f data/comparison_donor_anchors_data/*.tab.gz > /dev/null 2>&1")
        os.system("rm -f data/comparison_acceptor_anchors_data/*.tab.gz > /dev/null 2>&1")
        splicekit.core.features.make_counts_table("donor_anchors")
        splicekit.core.features.make_counts_table("acceptor_anchors")

rule junctions_count:
    input:
        expand("data/sample_junctions_data/sample_{sample}.tab.gz", sample=SAMPLES),
        "reference/junctions.tab.gz"
    output:
        "data/samples_junctions_counts.tab.gz"
    resources:
        mem = DEFAULT_MEM,
        time = DEFAULT_TIME,
        cores = DEFAULT_CORES
    run:
        import splicekit
        splicekit.core.annotation.make_comparisons()
        splicekit.core.junctions.junctions_per_sample()
        splicekit.core.features.load_genes()
        splicekit.core.features.read_junctions()
        splicekit.core.features.make_counts_table("junctions")

rule exons_count:
    input:
        expand("data/sample_exons_data/sample_{sample}.tab.gz", sample=SAMPLES),
        "reference/junctions.tab.gz"
    output:
        "data/samples_exons_counts.tab.gz"
    resources:
        mem = DEFAULT_MEM,
        time = DEFAULT_TIME,
        cores = DEFAULT_CORES
    run:
        import splicekit
        splicekit.core.annotation.make_comparisons()
        splicekit.core.features.load_genes()
        splicekit.core.features.read_exons()
        splicekit.core.features.make_counts_table("exons")

rule genes_count:
    input:
        expand("data/sample_genes_data/sample_{sample}.tab.gz", sample=SAMPLES),
        "reference/junctions.tab.gz"
    output:
        "data/samples_genes_counts.tab.gz"
    resources:
        mem = DEFAULT_MEM,
        time = DEFAULT_TIME,
        cores = DEFAULT_CORES
    run:
        import splicekit
        splicekit.core.annotation.make_comparisons()
        splicekit.core.features.load_genes()
        splicekit.core.features.read_genes()
        splicekit.core.features.make_counts_table("genes")

rule edgeR:
    input:
        "data/samples_{feature_type}_counts.tab.gz"
    output:
        "results/edgeR/{feature_type}/{comparison}_altsplice.tab.gz",
        "results/edgeR/{feature_type}/{comparison}_difffeature.tab.gz",
    wildcard_constraints:
        feature_type="junctions|exons|donor_anchors|acceptor_anchors"
    resources:
        mem = DEFAULT_MEM,
        time = DEFAULT_TIME,
        cores = DEFAULT_CORES
    run:
        splicekit.core.features.load_genes()
        splicekit.edgeR(wildcards.feature_type)

rule edgeR_genes:
    input:
        "data/samples_{feature_type}_counts.tab.gz"
    output:
        "results/edgeR/{feature_type}/{comparison}_difffeature.tab.gz",
    wildcard_constraints:
        feature_type="genes"
    resources:
        mem = DEFAULT_MEM,
        time = DEFAULT_TIME,
        cores = DEFAULT_CORES
    run:
        splicekit.core.features.load_genes()
        splicekit.edgeR(wildcards.feature_type)

rule juDGE:
    input:
        "results/edgeR/genes_results_complete.tab.gz",
        "results/edgeR/junctions_results_complete.tab.gz"
    output:
        "results/judge/scored.tab.gz"
    resources:
        mem = DEFAULT_MEM,
        time = DEFAULT_TIME,
        cores = DEFAULT_CORES
    shell:
        "splicekit judge"

rule juan:
    input:
        expand("results/edgeR/{feature_type}/{comparison}_altsplice.tab.gz", feature_type=["donor_anchors", "acceptor_anchors"], comparison = COMPARISONS)
    output:
        "results/edgeR/juan.done"
    resources:
        mem = DEFAULT_MEM,
        time = DEFAULT_TIME,
        cores = DEFAULT_CORES
    shell:
        "splicekit juan"
