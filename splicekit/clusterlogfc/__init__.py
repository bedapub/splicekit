import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import logomaker as lm
from scipy.spatial import distance_matrix
from numpy import linalg as LA
import os
import sys
import itertools
import splicekit
import seaborn as sns
import copy
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import dendrogram, linkage, cut_tree
from matplotlib import pyplot as plt
import matplotlib.font_manager
import seaborn as sns

def compute_distance(compound_name1, compound_name2, fname):
    if not os.path.exists(fname):
        return 0
    f_temp = open(fname)
    r = f_temp.readline()
    r = f_temp.readline()
    r = f_temp.readline()
    vector_x = []
    vector_y = []
    header = f_temp.readline().replace("\r", "").replace("\n", "").split("\t")
    r = f_temp.readline()
    while r:
        r = r.replace("\r", "").replace("\n", "").split("\t")
        data = dict(zip(header, r))
        vector_x.append(float(data[f"logFC_{compound_name1}"]))
        vector_y.append(float(data[f"logFC_{compound_name2}"]))
        r = f_temp.readline()
    f_temp.close()
    compound_distance = np.corrcoef(vector_x, vector_y)
    return compound_distance[0][1]

def make_distance(fdr_thr=0.05):
    fdr_string = str(fdr_thr).replace(".", "_")
    if not os.path.exists(f"results/clusterlogfc/data_fdr{fdr_string}"):
        os.makedirs(f"results/clusterlogfc/data_fdr{fdr_string}")
    comparisons = copy.copy(splicekit.core.annotation.comparisons)
    comparisons = [comp_id for (comp_id, _, _, _, _) in comparisons]
    
    dm = np.zeros([len(comparisons), len(comparisons)])
    distance_files = {}
    for feature_type in ["junctions", "exons", "donor_anchors", "acceptor_anchors", "genes"]:
        distance_files[feature_type] = open(f"results/clusterlogfc/distance_{feature_type}_fdr{fdr_string}.tab", "wt")
        distance_files[feature_type].write("\t".join(["compound1", "compound2", "distance"]) + "\n")
    m = len(comparisons)
    for index1, comparison_id1 in enumerate(comparisons):
        for index2, comparison_id2 in enumerate(comparisons):
            if not (index1<index2<m):
                continue
            print("[clusterlogfc] distance for:", comparison_id1, comparison_id2)
            if os.path.exists(f"data/{comparison_id1}_{comparison_id2}_junctions.tab"):
                continue
            # https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.distance.pdist.html#scipy.spatial.distance.pdist
            # Returns a condensed distance matrix Y. For each i and j (where i<j<m), where m is the number of original observations. The metric dist(u=X[i], v=X[j]) is computed and stored in entry m * i + j - ((i + 2) * (i + 1)) // 2
            mij = m * index1 + index2 - ((index1 + 2) * (index1 + 1)) // 2

            # no intersect
            os.system(f"splicecompare . {comparison_id1} {comparison_id1}_altsplice.tab . {comparison_id2} {comparison_id2}_altsplice.tab results/clusterlogfc/data_fdr{fdr_string}/{comparison_id1}_{comparison_id2} -fdr {fdr_thr}")
            os.system(f"splicecompare -features genes . {comparison_id1} {comparison_id1}_difffeature.tab . {comparison_id2} {comparison_id2}_difffeature.tab results/clusterlogfc/data_fdr{fdr_string}/{comparison_id1}_{comparison_id2} -fdr {fdr_thr}")
            for feature_type in ["junctions", "exons", "donor_anchors", "acceptor_anchors", "genes"]:
                compound_distance = compute_distance(comparison_id1, comparison_id2, f"results/clusterlogfc/data_fdr{fdr_string}/{comparison_id1}_{comparison_id2}_{feature_type}.tab")
                row = [comparison_id1, comparison_id2, compound_distance]
                distance_files[feature_type].write("\t".join([str(el) for el in row]) + "\n")
    for feature_type in ["junctions", "exons", "donor_anchors", "acceptor_anchors", "genes"]:
        distance_files[feature_type].close()

def make_cluster(fdr_thr=0.05):
    fdr_string = str(fdr_thr).replace(".", "_")
    sns.set(font_scale=2)
    plt.figure()
    comparisons = copy.copy(splicekit.core.annotation.comparisons)
    comparisons = [comp_id for (comp_id, _, _, _, _) in comparisons]

    for feature_type in ["junctions", "exons", "donor_anchors", "acceptor_anchors", "genes"]:
        f = open(f"results/clusterlogfc/distance_{feature_type}_fdr{fdr_string}.tab", "rt")
        header = f.readline()
        r = f.readline()
        cdm = []
        while r:
            r = r.replace("\r", "").replace("\n", "").split("\t")
            if r[-1]=="nan":
                cdm.append(0) # similarity = 0, non-similar, we don't know enough to say anything about this pair
            else:
                cdm.append(1-abs(float(r[-1]))) # similarity = 1-distance
            r = f.readline()
        f.close()

        Z = linkage(cdm, 'ward')
        cutree = cut_tree(Z)
        clusters = {}
        for compound_index, cluster_nr in enumerate(cutree):
            cluster = clusters.get(cluster_nr[0], [])
            cluster.append(comparisons[compound_index])
            clusters[cluster_nr[0]] = cluster

        fig, ax = plt.subplots(figsize=(15, 15))
        dn = dendrogram(Z, color_threshold=0.5, labels=comparisons, orientation="left", leaf_font_size=12, ax=ax)
        plt.setp(ax.collections, linewidth=3)

        plt.tight_layout()
        plt.savefig(f"results/clusterlogfc/dendrogram_{feature_type}_fdr{fdr_string}.png", dpi=300)
        plt.savefig(f"results/clusterlogfc/dendrogram_{feature_type}_fdr{fdr_string}.pdf")
        plt.close()

def process():
    make_distance()
    make_cluster()