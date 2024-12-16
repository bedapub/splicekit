import os
import pandas as pd
import splicekit

container: "docker://ghcr.io/bedapub/splicekit:main"
available_threads = workflow.cores

samples_df = pd.read_csv("samples.tab", sep="\t", comment="#")
SAMPLES = samples_df["sample_id"].tolist()

species = "homo_sapiens"
genome_version = None

if os.path.exists("splicekit.config"):
    with open("splicekit.config") as config_file:
        for line in config_file:
            exec(line.strip())

genome_version = f"--genome_version {genome_version}" if genome_version else ""

if not os.path.exists("annotation/comparisons.tab"):
    splicekit.annotation()

comparisons_df = pd.read_csv("annotation/comparisons.tab", sep="\t", comment="#")
COMPARISONS = comparisons_df["comparison"].tolist()

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

        # edgeR
        expand("results/edgeR/{feature_type}_results_complete.tab.gz", feature_type=["junctions", "exons", "donor_anchors", "acceptor_anchors", "genes"]),
        expand("results/edgeR/{feature_type}_results_fdr005.tab.gz", feature_type=["junctions", "exons", "donor_anchors", "acceptor_anchors", "genes"])

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
    threads: available_threads
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
        echo pybio star {species} {{input.fastq1}} {{input.fastq2}} {{output}} -t {{threads}} {genome_version}
        pybio star {species} {{input.fastq1}} {{input.fastq2}} {{output}} -t {{threads}} {genome_version}
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

rule junctions_per_sample:
    input:
        expand("data/sample_junctions_data/sample_{sample}_raw.tab.gz", sample=SAMPLES),
        "reference/junctions.tab.gz"
    output:
        expand("data/sample_junctions_data/sample_{sample}.tab.gz", sample=SAMPLES),
    run:
        import splicekit
        splicekit.core.annotation.make_comparisons()
        splicekit.core.junctions.junctions_per_sample()

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

rule exons_gtf:
    output:
        "reference/exons.gtf.gz",
    shell:
        """
        python -c 'import splicekit; splicekit.core.exons.write_exons_gtf()'
        """

rule genes_gtf:
    output:
        "reference/genes.gtf.gz",
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
    shell:
        """
        featureCounts {params.library_type_insert}-s {params.library_strand_insert} -M -O -T {threads} -F GTF -f -t exon -g exon_id -a {input.gtf_fname} -o {params.tab_fname} {input.bam_fname}
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
    shell:
        """
        featureCounts {params.library_type_insert}-s {params.library_strand_insert} -M -O -T {threads} -F GTF -f -t exon -g gene_id -a {input.gtf_fname} -o {params.tab_fname} {input.bam_fname}
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
    run:
        import splicekit
        splicekit.core.annotation.make_comparisons()
        splicekit.core.features.load_genes()
        splicekit.core.features.read_exons()
        splicekit.core.features.make_counts_table("exons")

rule genes_count:
    input:
        expand("data/sample_genes_data/sample_{sample}.tab.gz", sample=SAMPLES),
    output:
        "data/samples_genes_counts.tab.gz"
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
        feature_type="junctions"
    run:
        job_sh_edgeR="""R --no-save --args {input_folder} {atype} {control_name} {test_name} {comparison_name} {sample_membership} {filter_low} < {core_path}/comps_edgeR.R"""
        comparison = comps[wildcards.comparison]
        comparison_name, control_set, test_set, control_group_id, test_group_id = comps[wildcards.comparison]
        control_name = control_group_id
        test_name = test_group_id
        control_ids = []
        test_ids = []
        for (sample_id, compound, rep, _) in control_set:
            control_ids.append(sample_id)
        for (sample_id, compound, rep, _) in test_set:
            test_ids.append(sample_id)
        try:
            control_ids = sort_readout_id(control_ids)
        except:
            pass
        try:
            test_ids = sort_readout_id(test_ids)
        except:
            pass

        try:
            filter_low
        except:
            filter_low = "filter_low"

        job_sh_junctions = job_sh_edgeR.format(filter_low=filter_low, core_path=os.path.dirname(splicekit.core.__file__), comparison_name=comparison_name, input_folder=os.getcwd(), atype=wildcards.feature_type, control_name=control_name, test_name=test_name, sample_membership=",".join(str(el) for el in sample_membership))
        if workflow.use_singularity:
            job_sh_junctions = f"{container} {job_sh_junctions}"
        shell(job_sh_junctions)

rule edgeR_compile:
    input:
        expand("results/edgeR/{feature_type}/{comparison}_altsplice.tab.gz", comparison=COMPARISONS, feature_type=["junctions", "exons", "donor_anchors", "acceptor_anchors", "genes"]),
        expand("results/edgeR/{feature_type}/{comparison}_difffeature.tab.gz", comparison=COMPARISONS, feature_type=["junctions", "exons", "donor_anchors", "acceptor_anchors", "genes"])
    output:
        "results/edgeR/{feature_type}_results_complete.tab.gz",
        "results/edgeR/{feature_type}_results_fdr005.tab.gz"
    wildcard_constraints:
        feature_type="junctions"
    run:
        import splicekit
        splicekit.core.report.edgeR_feature(wildcards.feature_type)

"""

# alpha numeric sort
def break_readout_id(item):
    if type(item) is tuple:
        item_process = item[0]
    else:
        item_process = item
    for splitter in [" ", "_", "-"]:
        item_temp = item_process.split(splitter)
        # separates the string by splitter
        # find first element that casts to int and return (int, str_without_element)
        if len(item_temp)>1:
            for item_index_test in range(0, len(item_temp)):
                if item_temp[item_index_test].isdigit():
                    return (int(item_temp[item_index_test]), splitter.join(item_temp[:item_index_test]+item_temp[item_index_test+1:]))
    if item_process.isdigit():
        return int(item_process)
    else:
        return item_process

def sort_readout_id(data):
    return sorted(data, key=break_readout_id)

def make_short_names(text):
    result = text
    try:
        short_names
    except:
        return result
    for from_text, to_text, replace_type in short_names:
        if replace_type=="complete":
            if text==from_text:
                result = to_text
        else:
            result = text.replace(from_text, to_text)
    return result

def make_comparisons():
    comparisons = []
    comps = {}
    comparisons_ids = []
    treatments = {}
    samples = set()
    f = open("samples.tab")
    r = f.readline()
    while r.startswith("#"):
        r = f.readline()
    header = r.replace("\r", "").replace("\n", "").split("\t")
    r = f.readline()
    separates = set()
    while r:
        if r.startswith("#"):
            r = f.readline()
            continue
        r = r.replace("\r", "").replace("\n", "").split("\t")
        data = dict(zip(header, r))
        sample_id = data[sample_column]
        samples.add(sample_id)
        treatment_id = data[treatment_column]
        separates.add(data.get(separate_column, ""))
        treatments.setdefault(treatment_id, []).append((sample_id, treatment_id, data.get(separate_column, ""), data.get(group_column, "")))
        r = f.readline()
    f.close()
    for treatment, data in treatments.items(): # sort by sample ID
        try:
            treatments[treatment] = sort_readout_id(data)
        except:
            pass
    dmso_hash = {}
    dmso_letter = "A"
    separates = sorted(separates) # sometimes the separates set was reshufled, very difficult to pinpoint/debug this
    for sep in separates:
        for test_sample, test_data in treatments.items():
            if test_sample == control_name:
                continue
            temp_test = [(a,b,c,d) for (a,b,c,d) in test_data if c==sep]
            group_test = set([d for (a,b,c,d) in temp_test])
            if len(temp_test)==0:
                continue
            for control_sample, control_data in treatments.items():
                if control_sample != control_name:
                    continue
                temp_control = [(a,b,c,d) for (a,b,c,d) in control_data if c==sep]
                if group_column!="":
                    temp_control = [(a,b,c,d) for (a,b,c,d) in temp_control if d in group_test]
                if len(temp_control)==0:
                    continue
                # add comparison to the list
                # also add a unique control_group_id name, not only with _sep
                # this is then useful for the contrast and design matrix for DGE analysis
                if dmso_hash.get(tuple(temp_control), None)==None:
                    dmso_hash[tuple(temp_control)] = dmso_letter
                    dmso_letter = chr(ord(dmso_letter) + 1)
                sep_temp = sep.replace(" ", "_").lower()
                if sep_temp!="":
                    test_group_id = f"{test_sample}_{sep_temp}"
                    control_group_id = f"{control_sample}_{dmso_hash[tuple(temp_control)]}_{sep_temp}"
                    comparison_name = f"{test_sample}_{control_sample}{dmso_hash[tuple(temp_control)]}_{sep_temp}"
                else:
                    test_group_id = f"{test_sample}"
                    control_group_id = f"{control_sample}{dmso_hash[tuple(temp_control)]}"
                    comparison_name = f"{test_sample}_{control_sample}{dmso_hash[tuple(temp_control)]}"
                comparisons.append((make_short_names(comparison_name), temp_control, temp_test, make_short_names(control_group_id), make_short_names(test_group_id)))
                comps[comparisons[-1][0]] = comparisons[-1]
                comparisons_ids.append(make_short_names(comparison_name))
    return comps, comparisons, comparisons_ids

comps, comparisons, comparisons_ids = make_comparisons()

sample_membership = {}
for (comparison_name, control_set, test_set, control_group_id, test_group_id) in comparisons:
    for (sample_id, _, _, _) in control_set:
        # in some rare cases, the same sample can be part of diverse control groups
        sample_membership[sample_id] = control_group_id
    for (sample_id, _, _, _) in test_set:
        # in some rare cases, the same sample can be part of diverse test groups
        sample_membership[sample_id] = test_group_id
sample_membership = [sample_membership[sample_id] for sample_id in SAMPLES]

"""