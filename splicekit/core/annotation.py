# splicekit annotation module

import requests
import os
import subprocess
import sys
import splicekit.config as config
import splicekit.core.annotation as annotation
import splicekit.core as core
import splicekit


# dictionary of last nucleotide of first exon of transcripts
# used to annotate junctions hitting 5'UTR / alternative 5'UTR
# -> promoter changes
# annotation.first_exons[(chr, strand, pos)] = (transcript_id, gene_id)

def short_names(text):
    result = text
    if getattr(config, "short_names", None) == None:
        return result
    for from_text, to_text, replace_type in config.short_names:
        if replace_type=="complete":
            if text==from_text:
                result = to_text
        else:
            result = text.replace(from_text, to_text)
    return result

# alpha numeric sort
def break_readout_id(item):
    if type(item) is tuple:
        item_process = item[0]
    else:
        item_process = item
    for splitter in [" ", "_", "-"]:
        item_temp = item_process.split(splitter)
        if len(item_temp)>1:
            if item_temp[0].isdigit():
                return (int(item_temp[0]), splitter.join(item_temp[1:]))
    if item_process.isdigit():
        return int(item_process)
    else:
        return item_process

def sort_readout_id(data):
    return sorted(data, key=break_readout_id)

def to_int(el):
    try:
        return int(el)
    except:
        return el

class Cmd():
    # Interface for running commands from Python.
    def __init__(self, command):
        self.command = command
        self.returncode = None
        self.process = subprocess.Popen(['/bin/sh', '-c', command], stdout=subprocess.PIPE, stderr=subprocess.PIPE) # -cl and bash was original
        self.pid = self.process.pid

    def run(self):
        output, error = self.process.communicate()
        self.returncode = self.process.returncode
        return output.decode("utf-8"), error.decode("utf-8")

def read_comparisons():
    if not os.path.exists("samples.tab"):
        return
    annotation.comparisons = []
    annotation.treatments = {}
    samples = set()
    f = open("samples.tab")
    header = f.readline().replace("\r", "").replace("\n", "").split("\t")
    r = f.readline()
    separates = set()
    while r:
        if r.startswith("#"):
            r = f.readline()
            continue
        r = r.replace("\r", "").replace("\n", "").split("\t")
        data = dict(zip(header, r))
        sample_id = data[config.sample_column]
        samples.add(sample_id)
        treatment_id = data[config.treatment_column]
        separates.add(data.get(config.separate_column, ""))
        annotation.treatments.setdefault(treatment_id, []).append((sample_id, treatment_id, data.get(config.separate_column, ""), data.get(config.group_column, "")))
        r = f.readline()
    f.close()
    # sort by sample ID
    for treatment, data in annotation.treatments.items():
        annotation.treatments[treatment] = sort_readout_id(data)
    dmso_hash = {}
    dmso_letter = "A"
    # sometimes the separates set was reshufled
    # very difficult to pinpoint/debug this
    # sort it
    separates = sorted(separates)
    for sep in separates:
        for test_sample, test_data in annotation.treatments.items():
            if test_sample == config.control_name:
                continue
            temp_test = [(a,b,c,d) for (a,b,c,d) in test_data if c==sep]
            group_test = set([d for (a,b,c,d) in temp_test])
            if len(temp_test)==0:
                continue
            for control_sample, control_data in annotation.treatments.items():
                if control_sample != config.control_name:
                    continue
                temp_control = [(a,b,c,d) for (a,b,c,d) in control_data if c==sep]
                if config.group_column!="":
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
                    comparison_name = f"{test_sample}_{control_sample}_{dmso_hash[tuple(temp_control)]}_{sep_temp}"
                else:
                    test_group_id = f"{test_sample}"
                    control_group_id = f"{control_sample}_{dmso_hash[tuple(temp_control)]}"
                    comparison_name = f"{test_sample}_{control_sample}_{dmso_hash[tuple(temp_control)]}"
                annotation.comparisons.append((short_names(comparison_name), temp_control, temp_test, short_names(control_group_id), short_names(test_group_id)))
    # sort and cast sample_ids to string
    annotation.samples = list(samples)
    annotation.samples = sort_readout_id(annotation.samples)
    annotation.samples = [str(el) for el in annotation.samples]

def write_comparisons():
    job_edgeR="""
#!/bin/bash
#BSUB -J {job_name}                              # job name
#BSUB -n 4                                       # number of tasks
#BSUB -R "span[hosts=1]"                         # allocate hosts
#BSUB -M 16GB                                    # allocate memory
#BSUB -q short                                   # select queue
#BSUB -o logs/logs_edgeR_{atype}/{comp_name}.out # output file
#BSUB -e logs/logs_edgeR_{atype}/{comp_name}.err # error file

ml R
{container} R --no-save --args {input_folder} {data_folder} {atype} {control_name} {test_name} {comp_name} {control_list} {test_list} {filter_low} < {core_path}/comps_edgeR.R
"""

    job_sh_edgeR="""{container} R --no-save --args {input_folder} {data_folder} {atype} {control_name} {test_name} {comp_name} {control_list} {test_list} {filter_low} < {core_path}/comps_edgeR.R"""

    job_rmats="""
#!/bin/bash
#BSUB -J {job_name}                        # job name
#BSUB -n 4                                 # number of tasks
#BSUB -M 16GB                              # allocate memory
#BSUB -R "span[hosts=1]"                   # allocate hosts
#BSUB -q short                             # select queue
#BSUB -o logs/logs_rmats/{comp_name}.out   # output file
#BSUB -e logs/logs_rmats/{comp_name}.err   # error file

{container} run_rmats --b1 results/rmats/{comp_name}_test.tab --b2 results/rmats/{comp_name}_control.tab --gtf {gtf_path} -t paired --readLength 150 --variable-read-length --allow-clipping --nthread 4 --od results/rmats/{comp_name}_results --tmp results/rmats/{comp_name}_temp
"""

    comps_table = open("annotation/comparisons.tab", "wt")
    header = ["comparison", "control_samples", "test_samples", "control_group_id", "test_group_id"]
    comps_table.write("\t".join(header) + "\n")
    fsh_exons = open(f"jobs/jobs_edgeR_exons/process.sh", "wt")
    fsh_junctions = open(f"jobs/jobs_edgeR_junctions/process.sh", "wt")
    fsh_donor_anchors = open(f"jobs/jobs_edgeR_donor_anchors/process.sh", "wt")
    fsh_acceptor_anchors = open(f"jobs/jobs_edgeR_acceptor_anchors/process.sh", "wt")
    fsh_genes = open(f"jobs/jobs_edgeR_genes/process.sh", "wt")
    fsh_rmats = open(f"jobs/rmats/process.sh", "wt")
    for (comp_name, comp1, comp2, dmso_group_id, compound_group_id) in annotation.comparisons:
        comp1_compound = comp1[0][1]
        comp2_compound = comp2[0][1]
        fout_exons = open("jobs/jobs_edgeR_exons/{comp_name}.job".format(comp_name=comp_name), "wt")
        fout_junctions = open("jobs/jobs_edgeR_junctions/{comp_name}.job".format(comp_name=comp_name), "wt")
        fout_donor_anchors = open("jobs/jobs_edgeR_donor_anchors/{comp_name}.job".format(comp_name=comp_name), "wt")
        fout_acceptor_anchors = open("jobs/jobs_edgeR_acceptor_anchors/{comp_name}.job".format(comp_name=comp_name), "wt")
        fout_genes = open("jobs/jobs_edgeR_genes/{comp_name}.job".format(comp_name=comp_name), "wt")
        control_ids = []
        test_ids = []
        for (sample_id, compound, rep, _) in comp1:
            control_ids.append(sample_id)
        for (sample_id, compound, rep, _) in comp2:
            test_ids.append(sample_id)
        control_ids = sort_readout_id(control_ids)
        test_ids = sort_readout_id(test_ids)

        try:
            filter_low = splicekit.config.filter_low
        except:
            filter_low = "filter_low"

        # write rMATS {comp_name}_control.tab, {comp_name}_test.tab and {comp_name}_run.sh
        for rtype, rfile in [("control", control_ids), ("test", test_ids)]:
            bams = []
            for sample_id in rfile:
                bam_fname = os.path.abspath(os.path.join(f"{splicekit.config.bam_path}", f"{sample_id}.bam"))
                bams.append(bam_fname)
            f_rmats = open(f"results/rmats/{comp_name}_{rtype}.tab", "wt")
            f_rmats.write(",".join(bams))
            f_rmats.close()
        fsh_rmats.write(f"{config.container} run_rmats --b1 results/rmats/{comp_name}_test.tab --b2 results/rmats/{comp_name}_control.tab --gtf {splicekit.config.gtf_path[:-3]} -t paired --readLength 150 --variable-read-length --allow-clipping --nthread 4 --od results/rmats/{comp_name}_results --tmp results/rmats/{comp_name}_temp\n")
        f_rmats = open(f"jobs/rmats/{comp_name}.job", "wt")
        job_rmats_instance = job_rmats.format(container=splicekit.config.container, core_path=os.path.dirname(core.__file__), comp_name=comp_name, job_name="rmats_"+comp_name, gtf_path=splicekit.config.gtf_path[:-3])
        f_rmats.write(job_rmats_instance)
        f_rmats.close()

        # edgeR exons
        job_exons = job_edgeR.format(filter_low=filter_low, container=splicekit.config.container, core_path=os.path.dirname(core.__file__), comp_name=comp_name, input_folder=os.getcwd(), data_folder="data/comparison_exons_data", atype="exons", job_name="edgeR_exons_"+comp2_compound, control_name=comp1_compound, test_name=comp2_compound, control_list=",".join(str(el) for el in control_ids), test_list=",".join(str(el) for el in test_ids))
        fout_exons.write(job_exons)
        fout_exons.close()        
        job_sh_exons = job_sh_edgeR.format(filter_low=filter_low, container=splicekit.config.container, core_path=os.path.dirname(core.__file__), comp_name=comp_name, input_folder=os.getcwd(), data_folder="data/comparison_exons_data", atype="exons", control_name=comp1_compound, test_name=comp2_compound, control_list=",".join(str(el) for el in control_ids), test_list=",".join(str(el) for el in test_ids))
        fsh_exons.write(job_sh_exons+"\n")

        # edgeR junctions
        job_junctions = job_edgeR.format(filter_low=filter_low, container=splicekit.config.container, core_path=os.path.dirname(core.__file__), comp_name=comp_name, input_folder=os.getcwd(), data_folder="data/comparison_junctions_data", atype="junctions", job_name="edgeR_junctions_"+comp2_compound, control_name=comp1_compound, test_name=comp2_compound, control_list=",".join(str(el) for el in control_ids), test_list=",".join(str(el) for el in test_ids))
        fout_junctions.write(job_junctions)
        fout_junctions.close()
        job_sh_junctions = job_sh_edgeR.format(filter_low=filter_low, container=splicekit.config.container, core_path=os.path.dirname(core.__file__), comp_name=comp_name, input_folder=os.getcwd(), data_folder="data/comparison_junctions_data", atype="junctions", control_name=comp1_compound, test_name=comp2_compound, control_list=",".join(str(el) for el in control_ids), test_list=",".join(str(el) for el in test_ids))
        fsh_junctions.write(job_sh_junctions+"\n")

        # edgeR donor anchors
        job_donor_anchors = job_edgeR.format(filter_low=filter_low, container=splicekit.config.container, core_path=os.path.dirname(core.__file__), comp_name=comp_name, input_folder=os.getcwd(), data_folder="data/comparison_donor_anchors_data", atype="donor_anchors", job_name="edgeR_donor_anchors_"+comp2_compound, control_name=comp1_compound, test_name=comp2_compound, control_list=",".join(str(el) for el in control_ids), test_list=",".join(str(el) for el in test_ids))
        fout_donor_anchors.write(job_donor_anchors)
        fout_donor_anchors.close()
        job_sh_donor_anchors = job_sh_edgeR.format(filter_low=filter_low, container=splicekit.config.container, core_path=os.path.dirname(core.__file__), comp_name=comp_name, input_folder=os.getcwd(), data_folder="data/comparison_donor_anchors_data", atype="donor_anchors", control_name=comp1_compound, test_name=comp2_compound, control_list=",".join(str(el) for el in control_ids), test_list=",".join(str(el) for el in test_ids))
        fsh_donor_anchors.write(job_sh_donor_anchors+"\n")

        # edgeR acceptor anchors
        job_acceptor_anchors = job_edgeR.format(filter_low=filter_low, container=splicekit.config.container, core_path=os.path.dirname(core.__file__), comp_name=comp_name, input_folder=os.getcwd(), data_folder="data/comparison_acceptor_anchors_data", atype="acceptor_anchors", job_name="edgeR_acceptor_anchors_"+comp2_compound, control_name=comp1_compound, test_name=comp2_compound, control_list=",".join(str(el) for el in control_ids), test_list=",".join(str(el) for el in test_ids))
        fout_acceptor_anchors.write(job_acceptor_anchors)
        fout_acceptor_anchors.close()
        job_sh_acceptor_anchors = job_sh_edgeR.format(filter_low=filter_low, container=splicekit.config.container, core_path=os.path.dirname(core.__file__), comp_name=comp_name, input_folder=os.getcwd(), data_folder="data/comparison_acceptor_anchors_data", atype="acceptor_anchors", control_name=comp1_compound, test_name=comp2_compound, control_list=",".join(str(el) for el in control_ids), test_list=",".join(str(el) for el in test_ids))
        fsh_acceptor_anchors.write(job_sh_acceptor_anchors+"\n")

        # edgeR genes
        job_genes = job_edgeR.format(filter_low=filter_low, container=splicekit.config.container, core_path=os.path.dirname(core.__file__), comp_name=comp_name, input_folder=os.getcwd(), data_folder="data/comparison_genes_data", atype="genes", job_name="edgeR_genes_"+comp2_compound, control_name=comp1_compound, test_name=comp2_compound, control_list=",".join(str(el) for el in control_ids), test_list=",".join(str(el) for el in test_ids))
        fout_genes.write(job_genes)
        fout_genes.close()
        job_sh_genes = job_sh_edgeR.format(filter_low=filter_low, container=splicekit.config.container, core_path=os.path.dirname(core.__file__), comp_name=comp_name, input_folder=os.getcwd(), data_folder="data/comparison_genes_data", atype="genes", control_name=comp1_compound, test_name=comp2_compound, control_list=",".join(str(el) for el in control_ids), test_list=",".join(str(el) for el in test_ids))
        fsh_genes.write(job_sh_genes+"\n")

        row = [comp_name, ",".join(str(el) for el in test_ids), ",".join(str(el) for el in control_ids), compound_group_id, dmso_group_id]
        comps_table.write("\t".join(row) + "\n")
    comps_table.close()
    fsh_exons.close()
    fsh_junctions.close()
    fsh_donor_anchors.close()
    fsh_acceptor_anchors.close()
    fsh_rmats.close()

def write_comparisons_edgeR2():
    job_edgeR="""
#!/bin/bash
#BSUB -J {job_name}                                     # job name
#BSUB -n 4                                              # number of tasks
#BSUB -R "span[hosts=1]"                                # allocate hosts
#BSUB -M 16GB                                           # allocate memory
#BSUB -q short                                          # select queue
#BSUB -o logs/logs_edgeR2_{atype}/{comparison_name}.out # output file
#BSUB -e logs/logs_edgeR2_{atype}/{comparison_name}.err # error file

ml R
{container} R --no-save --args {input_folder} {atype} {control_name} {test_name} {comparison_name} {sample_membership} {filter_low} < {core_path}/comps_edgeR2.R
"""

    job_sh_edgeR="""{container} R --no-save --args {input_folder} {atype} {control_name} {test_name} {comparison_name} {sample_membership} {filter_low} < {core_path}/comps_edgeR2.R"""

    fsh_exons = open(f"jobs/jobs_edgeR2_exons/process.sh", "wt")
    fsh_junctions = open(f"jobs/jobs_edgeR2_junctions/process.sh", "wt")
    fsh_donor_anchors = open(f"jobs/jobs_edgeR2_donor_anchors/process.sh", "wt")
    fsh_acceptor_anchors = open(f"jobs/jobs_edgeR2_acceptor_anchors/process.sh", "wt")
    fsh_genes = open(f"jobs/jobs_edgeR2_genes/process.sh", "wt")
    sample_membership = {}
    for (comp_name, control_set, test_set, control_group_id, test_group_id) in annotation.comparisons:
        for (sample_id, _, _, _) in control_set:
            assert(sample_membership.get(sample_id, None)==None)
            sample_membership[sample_id] = control_group_id
        for (sample_id, _, _, _) in test_set:
            assert(sample_membership.get(sample_id, None)==None)
            sample_membership[sample_id] = test_group_id
    sample_membership = [sample_membership[sample_id] for sample_id in annotation.samples]
    print(sample_membership)
    for (comparison_name, control_set, test_set, test_group_id, control_group_id) in annotation.comparisons:
        #control_name = control_set[0][1]
        #test_name = test_set[0][1]
        control_name = control_group_id
        test_name = test_group_id
        fout_exons = open(f"jobs/jobs_edgeR2_exons/{comparison_name}.job", "wt")
        fout_junctions = open(f"jobs/jobs_edgeR2_junctions/{comparison_name}.job", "wt")
        fout_donor_anchors = open(f"jobs/jobs_edgeR2_donor_anchors/{comparison_name}.job", "wt")
        fout_acceptor_anchors = open(f"jobs/jobs_edgeR2_acceptor_anchors/{comparison_name}.job", "wt")
        fout_genes = open(f"jobs/jobs_edgeR2_genes/{comparison_name}.job", "wt")
        control_ids = []
        test_ids = []
        for (sample_id, compound, rep, _) in control_set:
            control_ids.append(sample_id)
        for (sample_id, compound, rep, _) in test_set:
            test_ids.append(sample_id)
        control_ids = sort_readout_id(control_ids)
        test_ids = sort_readout_id(test_ids)

        try:
            filter_low = splicekit.config.filter_low
        except:
            filter_low = "filter_low"

        # edgeR exons
        job_exons = job_edgeR.format(filter_low=filter_low, container=splicekit.config.container, core_path=os.path.dirname(core.__file__), comparison_name=comparison_name, input_folder=os.getcwd(), atype="exons", job_name="edgeR_exons_"+test_name, control_name=control_name, test_name=test_name, sample_membership=",".join(str(el) for el in sample_membership))
        fout_exons.write(job_exons)
        fout_exons.close()        
        job_sh_exons = job_sh_edgeR.format(filter_low=filter_low, container=splicekit.config.container, core_path=os.path.dirname(core.__file__), comparison_name=comparison_name, input_folder=os.getcwd(), atype="exons", control_name=control_name, test_name=test_name, sample_membership=",".join(str(el) for el in sample_membership))
        fsh_exons.write(job_sh_exons+"\n")

        # edgeR junctions
        job_junctions = job_edgeR.format(filter_low=filter_low, container=splicekit.config.container, core_path=os.path.dirname(core.__file__), comparison_name=comparison_name, input_folder=os.getcwd(), atype="junctions", job_name="edgeR_junctions_"+test_name, control_name=control_name, test_name=test_name, sample_membership=",".join(str(el) for el in sample_membership))
        fout_junctions.write(job_junctions)
        fout_junctions.close()
        job_sh_junctions = job_sh_edgeR.format(filter_low=filter_low, container=splicekit.config.container, core_path=os.path.dirname(core.__file__), comparison_name=comparison_name, input_folder=os.getcwd(), atype="junctions", control_name=control_name, test_name=test_name, sample_membership=",".join(str(el) for el in sample_membership))
        fsh_junctions.write(job_sh_junctions+"\n")

        # edgeR donor anchors
        job_donor_anchors = job_edgeR.format(filter_low=filter_low, container=splicekit.config.container, core_path=os.path.dirname(core.__file__), comparison_name=comparison_name, input_folder=os.getcwd(), atype="donor_anchors", job_name="edgeR_donor_anchors_"+test_name, control_name=control_name, test_name=test_name, sample_membership=",".join(str(el) for el in sample_membership))
        fout_donor_anchors.write(job_donor_anchors)
        fout_donor_anchors.close()
        job_sh_donor_anchors = job_sh_edgeR.format(filter_low=filter_low, container=splicekit.config.container, core_path=os.path.dirname(core.__file__), comparison_name=comparison_name, input_folder=os.getcwd(), atype="donor_anchors", control_name=control_name, test_name=test_name, sample_membership=",".join(str(el) for el in sample_membership))
        fsh_donor_anchors.write(job_sh_donor_anchors+"\n")

        # edgeR acceptor anchors
        job_acceptor_anchors = job_edgeR.format(filter_low=filter_low, container=splicekit.config.container, core_path=os.path.dirname(core.__file__), comparison_name=comparison_name, input_folder=os.getcwd(), atype="acceptor_anchors", job_name="edgeR_acceptor_anchors_"+test_name, control_name=control_name, test_name=test_name, sample_membership=",".join(str(el) for el in sample_membership))
        fout_acceptor_anchors.write(job_acceptor_anchors)
        fout_acceptor_anchors.close()
        job_sh_acceptor_anchors = job_sh_edgeR.format(filter_low=filter_low, container=splicekit.config.container, core_path=os.path.dirname(core.__file__), comparison_name=comparison_name, input_folder=os.getcwd(), atype="acceptor_anchors", control_name=control_name, test_name=test_name, sample_membership=",".join(str(el) for el in sample_membership))
        fsh_acceptor_anchors.write(job_sh_acceptor_anchors+"\n")

        # edgeR genes
        job_genes = job_edgeR.format(filter_low=filter_low, container=splicekit.config.container, core_path=os.path.dirname(core.__file__), comparison_name=comparison_name, input_folder=os.getcwd(), atype="genes", job_name="edgeR_genes_"+test_name, control_name=control_name, test_name=test_name, sample_membership=",".join(str(el) for el in sample_membership))
        fout_genes.write(job_genes)
        fout_genes.close()
        job_sh_genes = job_sh_edgeR.format(filter_low=filter_low, container=splicekit.config.container, core_path=os.path.dirname(core.__file__), comparison_name=comparison_name, input_folder=os.getcwd(), atype="genes", control_name=control_name, test_name=test_name, sample_membership=",".join(str(el) for el in sample_membership))
        fsh_genes.write(job_sh_genes+"\n")

    fsh_exons.close()
    fsh_junctions.close()
    fsh_donor_anchors.close()
    fsh_acceptor_anchors.close()

def make_design_contrast():
    # write design matrix
    f = open("annotation/design.tab", "wt")
    comparison_data = {}
    dmso_data = {}
    for (comp_id, dmso_list, test_list, compound_group_id, dmso_group_id) in splicekit.core.annotation.comparisons:
        comparison_data[compound_group_id] = test_list
        dmso_data[dmso_group_id] = dmso_list
    comparisons_list = list(comparison_data.keys())
    dmso_list = list(dmso_data.keys())
    header = comparisons_list+dmso_list
    master_order = header # contrast and design matrices rows and columns must be in the same order (they are not matched by name in the DGE script)
    f.write("\t".join([""]+header) + "\n")
    for sample_id in splicekit.core.annotation.samples:
        row = [sample_id]
        for comp_id, sample_list in comparison_data.items():
            sample_present = 0
            for (comparison_sample_id, _, _, _) in sample_list:
                if str(sample_id)==str(comparison_sample_id):
                    sample_present = 1
            row.append(sample_present)
        for dmso_id, dmso_list in dmso_data.items():
            sample_present = 0
            for (dmso_sample_id, _, _, _) in dmso_list:
                if str(sample_id)==str(dmso_sample_id):
                    sample_present = 1
            row.append(sample_present)
        f.write("\t".join([str(el) for el in row]) + "\n")
    f.close()

    # write contrast matrix
    f = open("annotation/contrast.tab", "wt")
    compound_groups = set()
    dmso_groups = set()
    for (comp_id, dmso_list, test_list, compound_group_id, dmso_group_id) in splicekit.core.annotation.comparisons:
        compound_groups.add(compound_group_id)
        dmso_groups.add(dmso_group_id)
    compound_groups = list(compound_groups)
    dmso_groups = list(dmso_groups)

    header = []
    for (comp_id, dmso_list, test_list, compound_group_id, dmso_group_id) in splicekit.core.annotation.comparisons:
        header.append(comp_id)

    f.write("\t".join([""]+header) + "\n")
    temp = {}
    for group_id in compound_groups + dmso_groups:
        row = [group_id]
        for (comp_id, dmso_list, test_list, compound_group_id, dmso_group_id) in splicekit.core.annotation.comparisons:
            if compound_group_id==group_id:
                row.append(1)
            elif dmso_group_id==group_id:
                row.append(-1)
            else:
                row.append(0)
        temp[group_id] = row
    for group_id in master_order:
        row = temp[group_id]
        f.write("\t".join([str(el) for el in row])+"\n")
    f.close()

def bam_count():
    return
    fout = open("annotation/bam_counts.tab", "wt")
    fout.write("\t".join([config.sample_column, "bam_count"]) + "\n")
    f = open("samples.tab")
    header = f.readline().replace("\r", "").replace("\n", "").split("\t")
    r = f.readline()
    while r:
        r = r.replace("\r", "").replace("\n", "").split("\t")
        data = dict(zip(header, r))
        bam_fname = os.path.join(config.bam_path, data[config.sample_column]+".bam")
        output = subprocess.check_output(f"samtools view -@ 4 -c {bam_fname}", shell=True)
        count = int(output.decode().split("\n")[0])
        print(f"splicekit | bam read counts: {data['readout_id']}, bam file {bam_fname}, counts = {count}")
        row = [data[config.sample_column], count]
        fout.write("\t".join([str(x) for x in row]) + "\n")
        r = f.readline()
    f.close()
    fout.close()
