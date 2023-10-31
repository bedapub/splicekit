"""
[motifs] module of splicekit

* retrieve sequence patterns (donor sites, classified by junction type) and construct sequence logos
* prepare fasta files for scanRBP and dreme
* generate html reports
"""

import pybio
import math
import random
random.seed(42)
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import logomaker as lm
from scipy.spatial import distance_matrix
from numpy import linalg as LA
import os
import Levenshtein as lev
import itertools
import splicekit
import splicekit.config as config
import sys
import pandas as pd
import seaborn as sns
import copy

logFC_thresh = 1 # threshold for events to be considered by effect size
motif_FDR = 0.05 # we use this threshold for motifs

# scan / motif areas
scanRBP_area = 100 # site-scanRBP_area .. site+-scanRBP_area
splice_sites_area = (-2, 6) # take -2 .. 6 around splice site (donor or acceptor)

# criteria is split by feature type, data read-in from "results/results_edgeR_{feature_type}_all.tab"
motif_criteria = {}

# (criteria_name, donor/acceptor, actual criteria)
# note that criteria_name must be unique across all feature types
motif_criteria["junctions"] = []
motif_criteria["junctions"].append(("jup", "donor", f"float(data['fdr'])<{motif_FDR} and float(data['donor_anchor_logfc'])>{logFC_thresh}"))
motif_criteria["junctions"].append(("jdown", "donor", f"float(data['fdr'])<{motif_FDR} and float(data['donor_anchor_logfc'])<-{logFC_thresh}"))
motif_criteria["junctions"].append(("jirr", "donor", f"float(data['fdr'])<{motif_FDR} and abs(float(data['donor_anchor_logfc']))<{logFC_thresh}"))
motif_criteria["junctions"].append(("jcontrols", "donor", f"float(data['fdr'])>{logFC_thresh}"))

# junctions_donor sites
motif_criteria["junctions"].append(("junctions_up_donor", "donor", f"float(data['fdr'])<{motif_FDR} and float(data['logFC'])>{logFC_thresh}"))
motif_criteria["junctions"].append(( "junctions_down_donor", "donor", f"float(data['fdr'])<{motif_FDR} and float(data['logFC'])<-{logFC_thresh}"))
motif_criteria["junctions"].append(("junctions_up_donor_controls", "donor", f"float(data['fdr'])>0.5 and float(data['logFC'])>0"))
motif_criteria["junctions"].append(( "junctions_down_donor_controls", "donor", f"float(data['fdr'])>0.5 and float(data['logFC'])<0"))

# junctions_acceptor sites
motif_criteria["junctions"].append(("junctions_up_acceptor", "acceptor", f"float(data['fdr'])<{motif_FDR} and float(data['logFC'])>{logFC_thresh}"))
motif_criteria["junctions"].append(( "junctions_down_acceptor", "acceptor", f"float(data['fdr'])<{motif_FDR} and float(data['logFC'])<-{logFC_thresh}"))
motif_criteria["junctions"].append(("junctions_up_acceptor_controls", "acceptor", f"float(data['fdr'])>0.5 and float(data['logFC'])>0"))
motif_criteria["junctions"].append(( "junctions_down_acceptor_controls", "acceptor", f"float(data['fdr'])>0.5 and float(data['logFC'])<0"))

# exons_donors sites
motif_criteria["exons"] = []
motif_criteria["exons"].append(("exons_up", "donor", f"float(data['fdr'])<{motif_FDR} and float(data['logFC'])>{logFC_thresh}"))
motif_criteria["exons"].append(("exons_down", "donor", f"float(data['fdr'])<{motif_FDR} and float(data['logFC'])<-{logFC_thresh}"))
motif_criteria["exons"].append(("exons_up_controls", "donor", f"float(data['fdr'])>0.5 and float(data['logFC'])>0"))
motif_criteria["exons"].append(("exons_down_controls", "donor", f"float(data['fdr'])>0.5 and float(data['logFC'])<0"))

# dreme positive sequences = criteria_name1
# dreme negative sequences = criteria_name2
# format of tuples (name is used for output folder name):
# (name, criteria_name1, criteria_name2)
dreme_criteria = {}
dreme_criteria = [("junctions_up_donor", "junctions_up_donor", "junction_up_donor_controls")]
dreme_criteria = [("junctions_down_donor", "junctions_down_donor", "junction_down_donor_controls")]
dreme_criteria = [("jup", "jup", "jcontrols")]
dreme_criteria = [("jdown", "jdown", "jcontrols")]
dreme_criteria = [("jirr", "jirr", "jcontrols")]

# pairs for scanRBP plots
# criteria_* refers to criteria_name
# (custom_name, criteria_up, criteria_down, criteria_up_control, criteria_down_control)
scanRBP_pairs = []
scanRBP_pairs.append(("junctions_donor", "junctions_up_donor", "junctions_down_donor", "junctions_up_donor_controls", "junctions_down_donor_controls")) # pair_name, motif_criteria_name_up, motif_criteria_name_down, motif_criteria_control_name_up, motif_criteria_control_name_down
scanRBP_pairs.append(("junctions_acceptor", "junctions_up_acceptor", "junctions_down_acceptor", "junctions_up_acceptor_controls", "junctions_down_acceptor_controls")) # pair_name, motif_criteria_name_up, motif_criteria_name_down, motif_criteria_control_name_up, motif_criteria_control_name_down

# check that the criteria names are unique
criteria_name_list = motif_criteria.values()
criteria_name_list = [item for sublist in criteria_name_list for item in sublist]
assert(len(set(criteria_name_list))==len(criteria_name_list))

def compute_mean_levdist(seq_list):   
    """
    given list of sequences, create list of all pairwise combinations (tuples)
    iterate over list and compute levensthein distance
    return average levensthein distance as float
    """
    # makes no biological sense to judge accuracy of seq_list of len > 3000 
    seq_list = seq_list[:3000]
    if len(seq_list)>1:
        seq_pairs = list(itertools.combinations(seq_list, 2))
        seq_dists = [lev.distance(pair[0], pair[1]) for pair in seq_pairs] 
        # normalize by lentgh of patters
        seqs_lens = [len(i) for i in seq_list]
        seqs_len_avg = sum(seqs_lens)/len(seqs_lens)
        seq_dist_avg = sum(seq_dists)/len(seq_dists)
        seq_dist_avg_norm = seq_dist_avg/seqs_len_avg
    else:
        seq_dist_avg_norm = None
    return seq_dist_avg_norm

def get_background():
    # equal background for all
    seqs = ["A" * (splice_sites_area[1]-splice_sites_area[0]+1)]
    seqs.append("T" * (splice_sites_area[1]-splice_sites_area[0]+1))
    seqs.append("C" * (splice_sites_area[1]-splice_sites_area[0]+1))
    seqs.append("G" * (splice_sites_area[1]-splice_sites_area[0]+1))
    counts_mat = lm.alignment_to_matrix(seqs)
    prob_mat = lm.transform_matrix(counts_mat, from_type='counts', to_type='probability')
    return len(seqs), prob_mat

def fix_matrix(mat):
    # make sure A, T, C, G is present in all matrices
    columns = ["A", "C", "G", "T"]
    for col in columns:
        if not col in mat:
            mat.insert({"A":0, "C":1, "G":2, "T":3}[col], col, [0]*len(mat))
    return mat

def mat_distance(matA, matB):
    matA = fix_matrix(matA)
    matB = fix_matrix(matB)
    return matA-matB

html_code = """
<html>
    <title>Donor site patterns by compound</title>
    <link rel="preconnect" href="https://fonts.googleapis.com">
    <link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
    <link href="https://fonts.googleapis.com/css2?family=Open+Sans:wght@300;400;700&display=swap" rel="stylesheet">
    <link href="fontawesome-free-6.1.2-web/css/all.css" rel="stylesheet">
<body>

<style>
    * {{
        font-family: "Open Sans", sans-serif;
        font-size: 11px;
        text-color: #f1f1f1;
    }}

</style>

Background ({num_background_sites}K donor sites)<br>
<img src='{background}' style='width:300px'>
<br><br>

<table border=0>
{image_blocks}
</table>

</body>
</html>
"""

scanRBP_data = {}
def make_scanRBP():
    print("[splicekit.core.motifs.make_scanRBP] start")

    scanRBP_to_process = [[a,b,c,d] for (_,a,b,c,d) in scanRBP_pairs]
    scanRBP_to_process = [item for sublist in scanRBP_to_process for item in sublist]
    comparison_counts = {}
    for feature_type, mcriteria in motif_criteria.items():
        print(f"[splicekit.core.motifs.make_scanRBP] processing {feature_type} criteria")
        scanRBP_data = {}
        f = open(f"results/results_edgeR_{feature_type}_all.tab", "rt")
        header = f.readline().replace("\r", "").replace("\n", "").split("\t")
        r = f.readline()
        while r:
            r = r.replace("\r", "").replace("\n", "").split("\t")
            data = dict(zip(header, r))
            fdr = float(data["fdr"])
            for criteria_name, seq_field, criteria in mcriteria:
                if criteria_name not in scanRBP_to_process:
                    continue
                try: # because some anchors are not reported (too low expression or other reasons), data["fdr_anchor"] is sometimes empty for a junction
                    criteria_met = eval(criteria)
                except:
                    criteria_met = False
                if criteria_met:
                    coords = data["feature_id"].split('_')
                    start = int(coords[-2])
                    stop = int(coords[-1])
                    strand = coords[-3][-1]
                    chr = '_'.join(coords[:-2])[:-1]
                    if comparison_counts.get((criteria_name, data["comparison"]), 0)<1000: # no more than 1000 sequences per criteria per comparison
                        comparison_counts[(criteria_name, data["comparison"])] = comparison_counts.setdefault((criteria_name, data["comparison"]), 0) + 1
                        if seq_field=="donor":
                            donor_site = start if strand=="+" and feature_type=="junctions" else stop
                            start_pos = donor_site - scanRBP_area
                            stop_pos = donor_site + scanRBP_area
                            donor_seq = pybio.core.genomes.seq_direct(config.species, chr, strand, start_pos, stop_pos, genome_version=config.genome_version)
                            fasta_id = f"{chr}{strand}_{start_pos}_{stop_pos}"
                            result_data = scanRBP_data.get((criteria_name, data["comparison"]), {})
                            result_data.setdefault(fasta_id, {})
                            result_data[fasta_id]["seq"] = donor_seq
                            result_data[fasta_id].setdefault("result_list", set()).add(data["result_id"])
                            scanRBP_data[(criteria_name, data["comparison"])] = result_data
                        if seq_field=="acceptor":
                            acceptor_site = stop if strand=="+" and feature_type=="junctions" else start
                            start_pos = acceptor_site - scanRBP_area
                            stop_pos = acceptor_site + scanRBP_area
                            acceptor_seq = pybio.core.genomes.seq_direct(config.species, chr, strand, start_pos, stop_pos, genome_version=config.genome_version)
                            fasta_id = f"{chr}{strand}_{start_pos}_{stop_pos}"
                            result_data = scanRBP_data.get((criteria_name, data["comparison"]), {})
                            result_data[fasta_id] = {}
                            result_data[fasta_id]["seq"] = acceptor_seq
                            result_data[fasta_id].setdefault("result_list", set()).add(data["result_id"])
                            scanRBP_data[(criteria_name, data["comparison"])] = result_data
            r = f.readline()
        f.close()

        for (criteria_name, comparison), result_data in scanRBP_data.items():
            if criteria_name not in scanRBP_to_process:
                continue
            print(f"[splicekit.core.motifs.make_scanRBP] {comparison}: {criteria_name} {feature_type}")
            fasta_fname = f"results/motifs/scanRBP/fasta/{comparison}_{criteria_name}_scanRBP.fasta"
            f = open(fasta_fname, "wt")
            for fasta_id, data in result_data.items():
                result_list = " ".join(list(data["result_list"]))
                seq = data["seq"]
                f.write(f">{fasta_id} {result_list}\n{seq}\n")
            f.close()
            os.system(f"scanRBP {fasta_fname} --output_folder results/motifs/scanRBP/data  -nonzero -simple_names")

def smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth

def bootstrap_logfc(matrix_signal, matrix_control, smoothing=6, iterations=1000):
    def matrix_vector(matrix):
        vector = matrix[0]
        for row in matrix[1:]:
            vector = [e1+e2 for e1,e2 in zip(vector, row)]
        vector = [e1/len(matrix) for e1 in vector] # normalize
        vector = smooth(vector, smoothing)[25:-25]
        return vector
    def boot():
        matrix_signal, matrix_control = [], []
        for index, value in enumerate(indices):
            if value==1:
                matrix_signal.append(master_matrix[index])
            else:
                matrix_control.append(master_matrix[index])
        vector_signal = matrix_vector(matrix_signal)
        vector_control = matrix_vector(matrix_control)
        if sum(vector_control)==0 or sum(vector_signal)==0:
            logfc = 0
        else:
            logfc = float(sum(vector_signal))/sum(vector_control)
            logfc = math.log(logfc, 2)
        boot_values.append(logfc)
    boot_values = []
    master_matrix = matrix_signal + matrix_control
    indices = [1]*len(matrix_signal) + [0]*len(matrix_control)
    boot()
    for i in range(iterations):
        random.shuffle(indices)
        boot()  
    return boot_values

def plot_scanRBP():
    print("[splicekit.core.motifs.plot_scanRBP] start")

    smoothing = 6
    protein = config.protein
    protein_label = config.protein_label

    def read_matrix_vector(fname, dtype="PWM"):
        matrix = []
        vector = []
        rows = []
        f = pybio.data.Fasta(fname)
        while f.read():
            id_original = f.id.split(" ")[0]
            id = f.id.replace("+", "plus").replace("-", "minus").split(" ")[0]
            if dtype=="PWM":
                scanRBP_fname = f"results/motifs/scanRBP/data/{id}_PWM.tab"
            elif dtype=="CLIP":
                scanRBP_fname = f"results/motifs/scanRBP/data/{id}_CLIP.tab"
            df = pd.read_csv(scanRBP_fname, sep='\t', header=0, index_col=0)
            data = list(df.loc[protein])
            #data = [x if x>=0 else 0 for x in data]
            data = [1 if x>0 else 0 for x in data]
            rows.append(id_original)
            matrix.append(data)
            if vector==[]:
                vector=data
            else:
                vector = [e1+e2 for e1,e2 in zip(vector, data)]
        return matrix, vector, rows

    for cdata in splicekit.core.annotation.comparisons:
        comparison = cdata[0]
        #for dtype in ["CLIP", "PWM"]:
        for dtype in ["PWM"]:
            for (cname, signal_up, signal_down, control_up, control_down) in scanRBP_pairs:
                up_fasta = f"results/motifs/scanRBP/fasta/{comparison}_{signal_up}_scanRBP.fasta"
                down_fasta = f"results/motifs/scanRBP/fasta/{comparison}_{signal_down}_scanRBP.fasta"
                upcontrol_fasta = f"results/motifs/scanRBP/fasta/{comparison}_{control_up}_scanRBP.fasta"
                downcontrol_fasta = f"results/motifs/scanRBP/fasta/{comparison}_{control_down}_scanRBP.fasta"
                if not os.path.exists(up_fasta) or not os.path.exists(down_fasta) or not os.path.exists(upcontrol_fasta) or not os.path.exists(downcontrol_fasta):
                    continue
                print(f"[splicekit.core.motifs.plot_scanRBP] processing {cname} {comparison} {dtype}")
                matrix_up, vector_up, rows_up = read_matrix_vector(up_fasta, dtype=dtype)
                matrix_down, vector_down, rows_down = read_matrix_vector(down_fasta, dtype=dtype)
                matrix_upcontrol, vector_upcontrol, rows_upcontrol = read_matrix_vector(upcontrol_fasta, dtype=dtype)
                matrix_downcontrol, vector_downcontrol, rows_downcontrol = read_matrix_vector(downcontrol_fasta, dtype=dtype)

                bootup = bootstrap_logfc(matrix_up, matrix_upcontrol, smoothing=smoothing, iterations=100000)
                value_up = bootup[0]
                bootup.sort(reverse=True)
                index_up = bootup.index(value_up)
                print("[splicekit.core.motifs.plot_scanRBP] up logFC=", value_up)
                print("[splicekit.core.motifs.plot_scanRBP] p-value up = ", index_up/float(len(bootup)))
                up_file = open(f"results/motifs/scanRBP/{comparison}_{signal_up}_bootstrap.tab", "wt")
                up_file.write(f"logFC\t{value_up}\np_value\t{index_up/float(len(bootup))}")
                up_file.close()

                bootdown = bootstrap_logfc(matrix_down, matrix_downcontrol, smoothing=smoothing, iterations=100000)
                value_down = bootdown[0]
                bootdown.sort(reverse=True)
                index_down = bootdown.index(value_down)
                print("[splicekit.core.motifs.plot_scanRBP] down logFC=", value_down)
                print("[splicekit.core.motifs.plot_scanRBP] p-value down = ", index_down/float(len(bootdown)))
                down_file = open(f"results/motifs/scanRBP/{comparison}_{signal_down}_bootstrap.tab", "wt")
                down_file.write(f"logFC\t{value_down}\np_value\t{index_down/float(len(bootdown))}")
                down_file.close()

                # heatmap
                """
                for htype, matrix_data, rows_data in [("up", matrix_up, rows_up), ("down", matrix_down, rows_down)]:
                    data = pd.DataFrame(matrix_data, index=rows_data)
                    plt.figure()
                    sns.set(font="Arial")
                    sns.set(font_scale=0.4)
                    sns.set_style("dark")
                    sns.set_style("ticks")
                    if htype=="up":
                        rdgn = sns.diverging_palette(h_neg=130, h_pos=10, s=99, l=55, sep=3, as_cmap=True)
                    if htype=="down":
                        rdgn = sns.diverging_palette(h_neg=250, h_pos=250, s=99, l=55, sep=3, as_cmap=True)
                    optional_params = {"linewidths":0.0}
                    optional_params["col_cluster"] = False
                    optional_params["annot"] = False
                    optional_params["center"] = 0
                    optional_params["fmt"] = ".0f"
                    optional_params["figsize"] = (10,3)
                    optional_params["yticklabels"] = True
                    optional_params["xticklabels"] = False
                    optional_params["cmap"] = rdgn
                    optional_params["cbar_pos"] = None
                    fig = sns.clustermap(data, **optional_params)
                    fig.ax_col_dendrogram.set_visible(False)
                    plt.tight_layout() 
                    plt.savefig(f"results/motifs/scanRBP/heatmap_{protein_label}_{comparison}_{dtype}_{htype}_{cname}.png", dpi=300)
                    plt.close()
                """

                plt.figure(figsize=(10,3))
                sns.set(font="Arial")
                sns.set(font_scale=0.6)
                sns.set_style("dark")
                sns.set_style("ticks")

                vector_up = [e1/len(matrix_up) for e1 in vector_up] # normalize
                vector_up = smooth(vector_up, smoothing)[25:-25]
                x = list(range(len(vector_up)))
                middle = int(len(x)/2)
                x = [el-middle for el in x]
                fig = sns.lineplot(data={"x":x, "y":vector_up}, x="x", y="y", linewidth=2, color='r')
                fig.fill_between(x, vector_up, color='#fb6767')
                fig.set(xlabel=f"distance from site [nt]", ylabel='fraction of sites (RBP binding)')
                fig.spines['left'].set_linewidth(0.5)
                fig.spines['left'].set_color('#333333')
                fig.spines['bottom'].set_linewidth(0.5)
                fig.spines['bottom'].set_color('#333333')
                fig.spines['top'].set_linewidth(0.5)
                fig.spines['top'].set_color('#333333')
                fig.spines['right'].set_linewidth(0.5)
                fig.spines['right'].set_color('#333333')
                fig.tick_params(axis='x', colors='#333333', width=0.5)
                fig.tick_params(axis='y', colors='#333333', width=0.5)

                vector_down = [e1/len(matrix_down) for e1 in vector_down] # normalize
                vector_down = smooth(vector_down, smoothing)[25:-25]
                max_val = max(vector_up)
                max_val = max(max_val, max(vector_down))
                max_val += max_val*0.01
                vector_down = [-e1 for e1 in vector_down]
                fig = sns.lineplot(data={"x":x, "y":vector_down}, x="x", y="y", linewidth=2, color='b')
                fig.fill_between(x, vector_down, color='#6767fb')

                # plot the controls
                vector_upcontrol = [e1/len(matrix_upcontrol) for e1 in vector_upcontrol] # normalize
                vector_upcontrol = smooth(vector_upcontrol, smoothing)[25:-25]
                fig = sns.lineplot(data={"x":x, "y":vector_upcontrol}, x="x", y="y", linewidth=2, color='#ffab40')
                #fig.fill_between(x, vector_upcontrol, color='#fb6767')

                vector_downcontrol = [e1/len(matrix_downcontrol) for e1 in vector_downcontrol] # normalize
                vector_downcontrol = smooth(vector_downcontrol, smoothing)[25:-25]
                vector_downcontrol = [-e1 for e1 in vector_downcontrol]
                fig = sns.lineplot(data={"x":x, "y":vector_downcontrol}, x="x", y="y", linewidth=2, color='#ffab40')
                #fig.fill_between(x, vector_downcontrol, color='#6767fb')

                plt.plot([-80, 80], [0, 0], color='#999999', linestyle='--', linewidth=0.3, alpha=0.5)
                plt.plot([0, 0], [-max_val, max_val], color='#999999', linestyle='--', linewidth=0.3, alpha=0.5)
                fig.set(ylim=(-max_val, max_val))
                fig.set(xlim=(-80, 80))
                plt.title(f"{dtype} {protein_label}, FDR<0.05, up={len(matrix_up)}, down={len(matrix_down)}, up_control={len(matrix_upcontrol)}, down_control={len(matrix_downcontrol)}, smoothing={smoothing}")
                plt.tight_layout() 
                plt.savefig(f"results/motifs/scanRBP/{protein_label}_{comparison}_{dtype}_{cname}.png", dpi=300)
                plt.close()

                """
                logfc = float(sum(vector_up))/sum(vector_upcontrol)
                logfc = math.log(logfc, 2)
                print(logfc)

                logfc = float(sum(vector_down))/sum(vector_downcontrol)
                logfc = math.log(logfc, 2)
                print(logfc)
                """


feature_added = {} # do not add sequences for the same donor/acceptor site multiple times (e.g. several junctions can represent the same donor site)
compound_seq_bytype = {}
def make_logos():
    print("[splicekit.core.motifs.make_logos] start")

    # evaluetes criteria from motif_criteria (see top of file) and distributes sequences accordingly
    # write fasta files
    # produces report in images and html

    for feature_type, mcriteria in motif_criteria.items():
        h_donors = {} # handle to fasta donor sequence files
        feature_file = f"results/results_edgeR_{feature_type}_all.tab"
        f = open(feature_file, "rt")
        header = f.readline().replace("\r", "").replace("\n", "").split("\t")
        r = f.readline()
        while r:
            r = r.replace("\r", "").replace("\n", "").split("\t")
            data = dict(zip(header, r))
            fdr = float(data["fdr"])
            coords = data["feature_id"].split('_')
            start = int(coords[-2])
            stop = int(coords[-1])
            strand = coords[-3][-1]
            chr = '_'.join(coords[:-2])[:-1]
            if feature_type in ["junctions"]:
                if strand=="+":
                    donor_site, acceptor_site = start, stop
                else:
                    donor_site, acceptor_site = stop, start
            if feature_type in ["exons", "donor_anchors"]:
                if strand=="+":
                    donor_site, acceptor_site = stop, start
                else:
                    donor_site, acceptor_site = start, stop
            donor_seq = pybio.core.genomes.seq(config.species, chr, strand, donor_site, splice_sites_area[0], splice_sites_area[1], genome_version=config.genome_version)
            if donor_seq.find("N")!=-1: # ignore sequences containing N
                r = f.readline()
                continue
            acceptor_seq = pybio.core.genomes.seq(config.species, chr, strand, acceptor_site, splice_sites_area[0], splice_sites_area[1], genome_version=config.genome_version)
            donor_id = f"{chr}{strand}_{donor_site}"
            acceptor_id = f"{chr}{strand}_{acceptor_site}"
            rank = data["rank"]
            if rank=="":
                rank = 0
            else:
                rank = int(rank)
            compound = str(rank)+"_"+data["compound"]
            ck = (rank, data["compound"]) # compound key (rank and compound name)

            # add control sequences to control fasta file
            if donor_id not in feature_added.get(f"{feature_type}_donor_controls", set()):
                fasta_id = "junction="+data["feature_id"]
                f_donor = h_donors.get(f"{feature_type}_donor_controls", open(f"results/motifs/fasta/{feature_type}_donor_controls.fasta", "wt"))
                f_donor.write(">{fasta_id}\n{donor_seq}\n".format(fasta_id=fasta_id, donor_seq=donor_seq))
                h_donors[f"{feature_type}_donor_controls"] = f_donor       

                feature_set = feature_added.get(f"{feature_type}_donor_controls", set())
                feature_set.add(donor_id)
                feature_added[f"{feature_type}_donor_controls"] = feature_set

            # add sequences to each motif_criteria, if the criteria is met (eval)
            for cryteria_name, seq_field, criteria in mcriteria:
                try: # because some anchors are not reported (too low expression or other reasons), data["fdr_anchor"] is sometimes empty for a junction
                    criteria_met = eval(criteria)
                except:
                    criteria_met = False

                #if criteria_met:
                if criteria_met and donor_id not in feature_added.get((cryteria_name, seq_field, data["comparison"]), set()):
                    temp = compound_seq_bytype.get(cryteria_name, {})
                    if seq_field=="donor":
                        temp_seq = donor_seq
                    if seq_field=="acceptor":
                        temp_seq = acceptor_seq
                    temp.setdefault(ck, []).append(temp_seq)
                    compound_seq_bytype[cryteria_name] = temp
                    # write sequence to fasta file
                    fasta_id = "junction="+data["feature_id"]
                    f_donor = h_donors.get((cryteria_name, data["comparison"]), open("results/motifs/fasta/"+data["comparison"]+"_" + cryteria_name + ".fasta", "wt"))
                    f_donor.write(">{fasta_id}\n{donor_seq}\n".format(fasta_id=fasta_id, donor_seq=temp_seq))
                    h_donors[(cryteria_name, data["comparison"])] = f_donor

                    feature_set = feature_added.get((cryteria_name, seq_field, data["comparison"]), set())
                    if seq_field=="donor":
                        feature_set.add(donor_id)
                    if seq_field=="acceptor":
                        feature_set.add(acceptor_id)
                    feature_added[(cryteria_name, seq_field, data["comparison"])] = feature_set


            temp = compound_seq_bytype.get(f"{feature_type}_all", {})
            temp.setdefault(ck, []).append(donor_seq)
            temp.setdefault((0, "overall"), []).append(donor_seq)
            compound_seq_bytype[f"{feature_type}_all"] = temp
            r = f.readline()
        f.close()
        for f_name, f_handle in h_donors.items():
            f_handle.close()

        print(f"[splicekit.core.motifs.make_logos] fasta files saved: {mcriteria}")
        to_remove = []
        for feature_type, compound_data in compound_seq_bytype.items():
            for ck, seqs in compound_data.items():
                if len(seqs)<3:
                    to_remove.append((feature_type, ck))
        for (feature_type, ck) in to_remove:
            del compound_seq_bytype[feature_type][ck]
        f.close()
    # end of reading data in from result files, now process and plot

    len_background, background = get_background()
    background_info = lm.transform_matrix(background, from_type='probability', to_type='information')
    for feature_type, compound_seq in compound_seq_bytype.items():
        m = len(compound_seq.keys())
        compound_seq = list(compound_seq.items())
        compound_seq.sort()
        logo = lm.Logo(background_info)
        logo.fig.savefig("results/motifs/images/background.png")
        plt.close()
        image_blocks = []
        image_blocks_md = []

        # write avg. levensthein dist per compound as rows to file per junction-type
        dist_fname = f"results/motifs/donor_pattern_per_compound_distance_{feature_type}.tab"
        dist_file = open(dist_fname, 'w')
        
        for (rank, compound_name), seqs in compound_seq:
            # compute per-compound and junction-type pairwise levensthein dist of patterns (so we can judge the accuracy of the splicing modifier) 
            seqs_mean_lev_pariwise_dist = compute_mean_levdist(seq_list=seqs)
            dist_file_row = f'{compound_name}\t{str(seqs_mean_lev_pariwise_dist)}\n'
            dist_file.write(dist_file_row)
            
            image_fname = "results/motifs/images/seqlogo_{feature_type}_rank{compound}.png".format(feature_type=feature_type, compound=f"{rank}_{compound_name}")
            image2_fname = "results/motifs/images/seqlogo2_{feature_type}_rank{compound}.png".format(feature_type=feature_type, compound=f"{rank}_{compound_name}")
            image_url = "images/seqlogo_{feature_type}_rank{compound}.png".format(feature_type=feature_type, compound=f"{rank}_{compound_name}")
            image2_url = "images/seqlogo2_{feature_type}_rank{compound}.png".format(feature_type=feature_type, compound=f"{rank}_{compound_name}")

            counts_mat = lm.alignment_to_matrix(seqs)
            counts_mat = fix_matrix(counts_mat)
            prob_mat = lm.transform_matrix(counts_mat, from_type='counts', to_type='probability')
            info_mat = lm.transform_matrix(counts_mat, from_type='counts', to_type='information')
            weight_mat = lm.transform_matrix(counts_mat, background=background, from_type='counts', to_type='weight')
            info_mat_back = lm.transform_matrix(weight_mat, from_type='weight', to_type='information')
            weight_mat_centered = lm.transform_matrix(weight_mat, center_values=True)
            weight_mat_info = lm.transform_matrix(weight_mat, from_type='weight', to_type='information')

            logo = lm.Logo(info_mat_back)
            logo.fig.savefig(image_fname)
            plt.close()

            logo = lm.Logo(info_mat)
            logo.fig.savefig(image2_fname)
            plt.close()

            seqs_mean_lev_pariwise_dist = "%.3f" % seqs_mean_lev_pariwise_dist
            image_blocks.append((rank, "<tr><td style='border-bottom: 1px dashed #aaaaaa; border-right: 1px dashed #cccccc;'>{compound} ({len_seq} sequences) (minus background)<br><img src='{image_fname}' style='width: 300px'></td><td style='border-bottom: 1px dashed #aaaaaa;'>{compound} ({len_seq} sequences, lev={seqs_mean_lev_pariwise_dist})<br><img src='{image2_fname}' style='width: 300px'></td></tr>".format(seqs_mean_lev_pariwise_dist=seqs_mean_lev_pariwise_dist, len_seq=len(seqs), image2_fname=image2_url, image_fname=image_url, compound=compound_name)))
            print("[splicekit.core.motifs.make_logos] seqlogo: feature_type={feature_type}, compound={compound_name}, sequences={len_seq}".format(feature_type=feature_type, compound_name=compound_name, len_seq=len(seqs)))

        image_blocks = [el[1] for el in image_blocks]
        f = open("results/motifs/index_{feature_type}.html".format(feature_type=feature_type), "wt")
        f.write(html_code.format(background="images/background.png", num_background_sites=round(len_background/1000), image_blocks="\n".join(image_blocks)))
        f.close()

        dist_file.close()

def dreme():
    print("[splicekit.core.motifs.dreme] start")
    for comparison in splicekit.core.annotation.comparisons:
        comp_id = comparison[0]
        for feature_type, mcriteria in motif_criteria.items():
            for (criteria_name, _, _) in mcriteria:
                if criteria_name.endswith("_controls"):
                    continue
                if os.path.exists(f"results/motifs/fasta/{feature_type}_donor_controls.fasta"):
                    if os.path.exists(f"results/motifs/fasta/{comp_id}_{criteria_name}.fasta"):
                        command = f"{splicekit.config.container} dreme -png -k 7 -norc -m 5 -p results/motifs/fasta/{comp_id}_{criteria_name}.fasta -n results/motifs/fasta/{feature_type}_donor_controls.fasta -oc results/motifs/dreme/{comp_id}_{criteria_name}"
                        os.system(command)

def make_distance():
    print("[splicekit.core.motifs.make_distance] start")
    len_background, background = get_background()
    for feature_type in motif_criteria.keys():
        compound_seq = compound_seq_bytype[f"{feature_type}_all"] # use all to do the clustering and dendrogram
        f = open(f"results/motifs/{feature_type}_donor_pattern_dm.tab", "wt")
        f.write("\t".join(["compound1", "compound2", "donnor_pattern_distance"]) + "\n")
        m = len(compound_seq.keys()) # 179 for ps180
        for index1, ((rank1, compound_name1), seqs1) in enumerate(compound_seq.items()):
            for index2, ((rank2, compound_name2), seqs2) in enumerate(compound_seq.items()):
                if not (index1<index2<m):
                    continue
                print("[splicekit.core.motifs.make_distance] distance for:", compound_name1, compound_name2)
                # https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.distance.pdist.html#scipy.spatial.distance.pdist
                # Returns a condensed distance matrix Y. For each i and j (where i<j<m), where m is the number of original observations. The metric dist(u=X[i], v=X[j]) is computed and stored in entry m * i + j - ((i + 2) * (i + 1)) // 2
                mij = m * index1 + index2 - ((index1 + 2) * (index1 + 1)) // 2
                counts_mat1 = lm.alignment_to_matrix(seqs1)
                counts_mat1 = fix_matrix(counts_mat1)
                prob_mat1 = lm.transform_matrix(counts_mat1, from_type='counts', to_type='probability')
                weight_mat1 = lm.transform_matrix(counts_mat1, background=background, from_type='counts', to_type='weight')
                weight_mat1_info = lm.transform_matrix(weight_mat1, from_type='weight', to_type='information')
                counts_mat2 = lm.alignment_to_matrix(seqs2)
                counts_mat2 = fix_matrix(counts_mat2)
                prob_mat2 = lm.transform_matrix(counts_mat2, from_type='counts', to_type='probability')
                weight_mat2 = lm.transform_matrix(counts_mat2, background=background, from_type='counts', to_type='weight')
                weight_mat2_info = lm.transform_matrix(weight_mat2, from_type='weight', to_type='information')
                dm = mat_distance(prob_mat1, prob_mat2)
                row = [compound_name1, compound_name2, LA.norm(dm, 'fro')]
                f.write("\t".join([str(el) for el in row]) + "\n")
        f.close()

def cluster(cutoff=9):
    print("[splicekit.core.motifs.cluster] start")
    from scipy.spatial.distance import pdist
    from scipy.cluster.hierarchy import dendrogram, linkage, cut_tree
    from matplotlib import pyplot as plt
    for feature_type in motif_criteria.keys():
        compounds = set()
        for ck, seqs in compound_seq_bytype[f"{feature_type}_all"].items():
            compounds.add(ck)

        compounds = list(compounds)
        compounds.sort()
        print("[splicekit.core.motifs.cluster] all compounds = ", len(compounds))
        print("[splicekit.core.motifs.cluster] cluster cutoff = ", cutoff)

        f = open(f"results/motifs/{feature_type}_donor_pattern_dm.tab", "rt")
        header = f.readline()
        r = f.readline()
        cdm = []
        while r:
            r = r.replace("\r", "").replace("\n", "").split("\t")
            cdm.append(float(r[-1]))
            r = f.readline()
        f.close()

        Z = linkage(cdm, 'ward')
        cutree = cut_tree(Z, n_clusters=cutoff) # cut dendrogram at point where we have args.cutoff clusters
        clusters = {}
        for compound_index, cluster_nr in enumerate(cutree):
            cluster = clusters.get(cluster_nr[0], [])
            cluster.append(compounds[compound_index])
            clusters[cluster_nr[0]] = cluster

        fig = plt.figure(figsize=(10, 25))
        dn = dendrogram(Z, color_threshold=0.9, labels=compounds, truncate_mode="lastp", p=cutoff)
        plt.tight_layout()
        plt.savefig(f"results/motifs/{feature_type}_dendrogram_truncated.png", dpi=300)
        plt.close()

        fig = plt.figure(figsize=(10, 25))
        dn = dendrogram(Z, color_threshold=0.9, labels=compounds, orientation="left")
        plt.tight_layout()
        plt.savefig(f"results/motifs/{feature_type}_dendrogram.png", dpi=300)
        plt.close()

        # print out clusters
        html_code = """
        <html>
            <title>splicekit: motifs</title>
            <link rel="preconnect" href="https://fonts.googleapis.com">
            <link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
            <link href="https://fonts.googleapis.com/css2?family=Open+Sans:wght@300;400;700&display=swap" rel="stylesheet">
            <link href="fontawesome-free-6.1.2-web/css/all.css" rel="stylesheet">
        <body>

        <style>
            * {{
                font-family: "Open Sans", sans-serif;
                font-size: 11px;
                text-color: #f1f1f1;
            }}
        </style>

        <center>

        <table>
        <tr>
        <td style='vertical-align: top;'><a href="dendrogram.png"><img src="dendrogram.png" style="width:350px;"></a></td>
        <td style='vertical-align: top;'>
            <table style='border-left: 1px dashed #555555;'>
            {cluster_blocks}
            </table>
        </td>
        </tr>
        </table>

        </body>
        </html>
        """

        cluster_blocks = []
        cluster_row = []
        for cluster_id, cluster_data in clusters.items():
            if len(cluster_data)<5:
                continue
            print("[splicekit.core.motifs.cluster] cluster id = ", cluster_id, "size = ", len(cluster_data))
            temp = ["<td style='vertical-align: top; border-bottom: 1px dashed #aaaaaa; padding-bottom: 10px; border-right: 1px dashed #cccccc;'>"]
            temp.append("<b>cluster</b> id {cluster}, <b>size = {size}</b> (displaying first 3 compounds)<br>".format(size=len(cluster_data), cluster=cluster_id))
            for (rank, compound) in cluster_data[:3]:
                image_fname = f"images/seqlogo2_{feature_type}_rank{rank}_{compound}.png"
                temp.append(f"{compound} (rank {rank})<br><a href='{image_fname}'><img src='{image_fname}' style='width: 300px'></a><br>")
            temp.append("</td>")
            cluster_row.append("\n".join(temp)+"\n\n")
            if len(cluster_row)==3:
                cluster_blocks.append("\n".join(["<tr>"] + cluster_row + ["</tr>"])+"\n\n")
                cluster_row = []
        cluster_blocks.append("\n".join(["<tr>"] + cluster_row + ["</tr>"])+"\n\n")    

        f = open("results/motifs/clusters.html", "wt")
        f.write(html_code.format(cluster_blocks="\n".join(cluster_blocks)))
        f.close()

def process():
    if config.scanRBP:
        make_scanRBP()
        plot_scanRBP()
    make_logos()
    dreme()
    make_distance()
    cluster()
