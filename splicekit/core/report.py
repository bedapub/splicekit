"""
# Description
Reads edgeR results files and generates tab reports.
"""

import os
import sys
import gzip
import pybio
import glob
import splicekit.config as config
import splicekit.core.annotation as annotation
import splicekit.core.features as features
import numpy as np

def toint(temp):
    try:
        temp = str(int(temp))
        return temp
    except:
        return temp

def edgeR_feature(feature_name):
    # feature_name = "genes", "exons", "junctions", "donor_anchors", "acceptor_anchors"
    print("Generating edgeR results file={fname}, considering all results with FDR<={edgeR_FDR_thr}".format(fname=f"results/results_edgeR_{feature_name}.tab", edgeR_FDR_thr=config.edgeR_FDR_thr))
    print()
    comparisons = {}

    database = {}
    # read junctions reference
    if feature_name=="junctions":
        features.load_genes() # for loading promoters
        f = open("reference/junctions.tab", "rt")
        header = f.readline().replace("\r", "").replace("\n", "").split("\t")
        r = f.readline()
        while r:
            r = r.replace("\r", "").replace("\n", "").split("\t")
            data = dict(zip(header, r))
            database[data["junction_id"]] = data
            r = f.readline()
        f.close()        

    # since we construct {sample_id}_{treatment_id} and {treatment_id} can still contain _,
    # we go back to the original samples table to extract {sample_id}
    # since the column can also be named otherwise (not only sample_id), we use the correct name from config.sample_column
    samples = []
    f = open("samples.tab", "rt")
    header = f.readline().replace("\r", "").replace("\n", "").split("\t")
    r = f.readline()
    while r:
        r = r.replace("\r", "").replace("\n", "").split("\t")
        data = dict(zip(header, r))
        samples.append(data[config.sample_column])
        r = f.readline()
    f.close()

    # read comparisons, to determine which tracks (samples) to show in JBrowse
    f = open("annotation/comparisons.tab", "rt")
    r = f.readline()
    r = f.readline()
    while r:
        compound_samples = []
        dmso_samples = []
        r = r.replace("\r", "").replace("\n", "").split("\t")
        comparison = r[0]
        temp = r[1].split(",")
        for el in temp:
            my_sample = ""
            for sample_id in samples:
                if el.startswith(sample_id+"_") and len(my_sample)<len(sample_id):
                    my_sample = sample_id
            if my_sample not in compound_samples and my_sample!="":
                compound_samples.append(my_sample)
        temp = r[2].split(",")
        for el in temp:
            my_sample = ""
            for sample_id in samples:
                if el.startswith(sample_id+"_") and len(my_sample)<len(sample_id):
                    my_sample = sample_id
            if my_sample not in dmso_samples and my_sample!="":
                dmso_samples.append(my_sample)
        comparisons[comparison] = (compound_samples, dmso_samples)
        r = f.readline()
    f.close()

    if feature_name=="genes":
        rfiles = glob.glob(f"results/results_edgeR_{feature_name}/*difffeature.tab")
    else:
        rfiles = glob.glob(f"results/results_edgeR_{feature_name}/*altsplice.tab")
    results = [] # only results with FDR<config.edgeR_FDR_thr, stored in results_edgeR_{feature}.tab
    results_all = [] # all results, stored in results_edgeR_{feature}_all.tab
    count = 0
    for fname in rfiles:
        bname = os.path.basename(fname)
        if feature_name=="genes":
            comparison = "_".join(bname.split("_difffeature.tab")[:-1])
        else:
            comparison = "_".join(bname.split("_altsplice.tab")[:-1])
        compound = bname.split("_"+config.control_name)[0]
        f = open(fname, "rt")
        header = f.readline().replace("\n", "").replace("\r", "").split("\t")
        r = f.readline()
        while r:
            r = r.replace("\n", "").replace("\r", "").split("\t")
            data = dict(zip(header, r))
            tracks_fixed = ["annotation_track"] # New JBrowse2
            tracks_samples, tracks_dmso = comparisons[comparison] 
            tracks_samples = [track+'_bw' for track in tracks_samples] # New JBrowse2 --> bigwig files are called by id_bw
            tracks_dmso = [track+'_bw' for track in tracks_dmso] # New JBrowse2 --> bigwig files are called by id_bw
            tracks = ",".join(tracks_samples[:3]+tracks_dmso[:3]+tracks_fixed)
            f_from=int(data["feature_start"])
            f_to=int(data["feature_stop"])
            delta = round(abs(f_from-f_to)*0.15) # 15% surrounding
            loc_from=int(data["feature_start"])-delta
            loc_to=int(data["feature_stop"])+delta
            row = [comparison, compound, "", data["feature_id"]]
            row += [data["chr"], data["strand"], data["feature_start"], data["feature_stop"], data["length"]]
            junction_start = f_from if data["strand"]=="+" else f_to
            first_exon_data = annotation.first_exons.get((data["chr"], data["strand"], junction_start), None)
            if first_exon_data!=None:
                transcript_id, gene_id, first_exon_start, first_exon_stop = first_exon_data
                junction_first_exon = f"first_exon_start_{first_exon_start}"
            else:
                first_exon_start = None
                junction_first_exon = ""
            if feature_name=="junctions":
                row.append(junction_first_exon)
            row += [data["gene_id"], data["gene_name"], "{chr}:{f_from}..{f_to}".format(chr=data["chr"], f_from=data["feature_start"], f_to=data["feature_stop"])]
            row.append("{jbrowse_url}&assembly=GenomeSequence&loc={chr}:{loc_from}..{loc_to}&tracks={tracks}".format(jbrowse_url=config.jbrowse2_url.format(compound=compound), compound=compound, chr=data["chr"], loc_from=loc_from, loc_to=loc_to, tracks=tracks))
            row += [data["sum_feature_test"], data["sum_feature_control"]]
            row += [data["test_pfi"], data["control_pfi"], data["delta_pfi"]]
            if feature_name=="junctions":
                row.append(database[data["feature_id"]]["annotated"]) 
                row.append(database[data["feature_id"]]["donor_anchor_id"])
                row.append(database[data["feature_id"]]["acceptor_anchor_id"])
                row.append(data["donor_DAI"])
                row.append(data["acceptor_DAI"])
                row.append(data["delta_DAI"])
                row.append(data["delta_DAI_pvalue"])
            if feature_name=="exons":
                row.append(data["test_PSI"])
                row.append(data["control_PSI"])
                row.append(data["delta_PSI"])
            if feature_name=="genes":
                row += [float(data["logFC"]), data["F"], float(data["PValue"]), float(data["FDR"])]
                row += [float(data["logFC"])* -np.log10(float(data["PValue"]))] # add pi value

            else:
                row += [float(data["logFC"]), data["exon.F"], float(data["P.Value"]), float(data["FDR"])]
                row += [float(data["logFC"])* -np.log10(float(data["P.Value"]))] # add pi value
            results_all.append(row)
            if float(data["FDR"])<=config.edgeR_FDR_thr:
                results.append(row)
            r = f.readline()
        f.close()
        count += 1
        print("Processed", fname, "{a}/{b} ({pdone}% done)".format(a=count, b=len(rfiles), pdone="%.2f" % (count/len(rfiles)*100)))

    headers = {}
    headers["genes"] = ["result_id", "comparison", "compound", "rank", "feature_id", "chr", "strand", "feature_start", "feature_stop", "feature_len", "gene_id", "gene_name", "jbrowse_loc", "jbrowse_url", "sum_feature_test", "sum_feature_control", "test_pfi", "control_pfi", "delta_pfi", "logFC", "exon.F", "p_value", "fdr","pi_value"]
    headers["junctions"] = ["result_id", "comparison", "compound", "rank", "feature_id", "chr", "strand", "feature_start", "feature_stop", "feature_len", "UTR", "gene_id", "gene_name", "jbrowse_loc", "jbrowse_url", "sum_feature_test", "sum_feature_control", "test_pfi", "control_pfi", "delta_pfi", "annotated", "donor_anchor_id", "acceptor_anchor_id", "donor_DAI", "acceptor_DAI", "delta_DAI", "delta_DAI_pvalue", "logFC", "exon.F", "p_value", "fdr","pi_value"]
    headers["exons"] = ["result_id", "comparison", "compound", "rank", "feature_id", "chr", "strand", "feature_start", "feature_stop", "feature_len", "gene_id", "gene_name", "jbrowse_loc", "jbrowse_url", "sum_feature_test", "sum_feature_control", "test_pfi", "control_pfi", "delta_pfi", "test_PSI", "control_PSI", "delta_PSI", "logFC", "exon.F", "p_value", "fdr","pi_value"]
    headers["donor_anchors"] = ["result_id", "comparison", "compound", "rank", "feature_id", "chr", "strand", "feature_start", "feature_stop", "feature_len", "gene_id", "gene_name", "jbrowse_loc", "jbrowse_url", "sum_feature_test", "sum_feature_control", "test_pfi", "control_pfi", "delta_pfi", "logFC", "exon.F", "p_value", "fdr","pi_value"]
    headers["acceptor_anchors"] = ["result_id", "comparison", "compound", "rank", "feature_id", "chr", "strand", "feature_start", "feature_stop", "feature_len", "gene_id", "gene_name", "jbrowse_loc", "jbrowse_url", "sum_feature_test", "sum_feature_control", "test_pfi", "control_pfi", "delta_pfi", "logFC", "exon.F", "p_value", "fdr","pi_value"]

    # sort by FDR
    FDR_index = headers[feature_name].index("fdr")-1 # -1 because of result_id
    results_all.sort(key = lambda el: float(el[FDR_index]))
    results.sort(key = lambda el: float(el[FDR_index]))

    f = open(f"results/results_edgeR_{feature_name}.tab", "wt")
    f.write("\t".join(headers[feature_name]) + "\n")
    f_all = open(f"results/results_edgeR_{feature_name}_all.tab", "wt")
    f_all.write("\t".join(headers[feature_name]) + "\n")
    if len(results_all)>0:
        assert(len(headers[feature_name])==len(results_all[0])+1) # header and data columns must match, +1 for result_id
    result_id = 1
    for row in results:
        f.write("\t".join(["r{result_id}".format(result_id=result_id)] + [str(el) for el in row]) + "\n")
        result_id += 1
    f.close()
    result_id = 1
    for row in results_all:
        f_all.write("\t".join(["ra{result_id}".format(result_id=result_id)] + [str(el) for el in row]) + "\n")
        result_id += 1
    f_all.close()
