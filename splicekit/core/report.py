# reads edgeR results files and generates TAB reports

import os
import sys
import gzip
import pybio
import glob
import splicekit.config as config
if config.jbrowse2_url==None:
    config.jbrowse2_config()
from splicekit.core import smart_number_format
import splicekit.core.annotation as annotation
import splicekit.core.features as features
import numpy as np

module_name = "splicekit | report |"

def edgeR_feature(feature_name, version=""):
    # feature_name = "genes", "exons", "junctions", "donor_anchors", "acceptor_anchors"
    fname = f"results/results_edgeR{version}_{feature_name}.tab.gz"
    print(f"{module_name} generating edgeR results file={fname}, considering all results with FDR<={config.edgeR_FDR_thr}")
    print()
    comparisons = {}

    database = {}
    # read junctions reference
    if feature_name=="junctions":
        features.load_genes() # for loading promoters
        f = gzip.open("reference/junctions.tab.gz", "rt")
        header = f.readline().replace("\r", "").replace("\n", "").split("\t")
        r = f.readline()
        while r:
            r = r.replace("\r", "").replace("\n", "").split("\t")
            data = dict(zip(header, r))
            database[data["junction_id"]] = data
            r = f.readline()
        f.close()        

    samples = []
    f = open("samples.tab", "rt")
    r = f.readline()
    while r.startswith("#"):
        r = f.readline()
    header = r.replace("\r", "").replace("\n", "").split("\t")
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
        control_samples = []
        test_samples = []
        r = r.replace("\r", "").replace("\n", "").split("\t")
        comparison = r[0]
        for el in r[1].split(","):
            control_samples.append(el)
        for el in r[2].split(","):
            test_samples.append(el)
        comparisons[comparison] = (control_samples, test_samples)
        r = f.readline()
    f.close()

    if feature_name=="genes":
        rfiles = glob.glob(f"results/edgeR/{feature_name}/*difffeature.tab.gz")
    else:
        rfiles = glob.glob(f"results/edgeR/{feature_name}/*altsplice.tab.gz")
    results = [] # only results with FDR<config.edgeR_FDR_thr, stored in results_edgeR_{feature}.tab
    results_all = [] # all results, stored in results_edgeR_{feature}_all.tab
    count = 0
    for fname in rfiles:
        bname = os.path.basename(fname)
        if feature_name=="genes":
            comparison = "_".join(bname.split("_difffeature.tab.gz")[:-1])
        else:
            comparison = "_".join(bname.split("_altsplice.tab.gz")[:-1])
        treatment = bname.split("_"+config.control_name)[0]
        f = gzip.open(fname, "rt")
        header = f.readline().replace("\n", "").replace("\r", "").split("\t")
        r = f.readline()
        while r:
            r = r.replace("\n", "").replace("\r", "").split("\t")
            data = dict(zip(header, r))
            tracks_fixed = ["annotation_track"] # New JBrowse2
            assembly = "GenomeSequence"
            if config.jbrowse2_url.startswith("https://genomebrowser"):
                assembly = {"homo_sapiens":"hg38", "mus_musculus":"mm39"}[config.species]
                tracks_fixed = ["gene-refseq", "transcript-refseq", "transcript-ensembl"]
            tracks_control, tracks_test = comparisons[comparison] 

            if config.jbrowse2_url.startswith("https://genomebrowser"):
                tracks_control = [track.lstrip("0")+'_bw' for track in tracks_control]
                tracks_test = [track.lstrip("0")+'_bw' for track in tracks_test]
            else:
                tracks_control = [track+'_bw' for track in tracks_control] # New JBrowse2 --> bigwig files are called by id_bw
                tracks_test = [track+'_bw' for track in tracks_test] # New JBrowse2 --> bigwig files are called by id_bw

            tracks = ",".join(tracks_test[:3]+tracks_control[:3]+tracks_fixed)
            try:
                for r1, r2 in config.track_name_replace:
                    tracks = tracks.replace(r1, r2)
            except:
                pass
            f_from=int(data["feature_start"])
            f_to=int(data["feature_stop"])
            delta = round(abs(f_from-f_to)*0.15) # 15% surrounding
            loc_from=int(data["feature_start"])-delta
            loc_to=int(data["feature_stop"])+delta
            row = [comparison, treatment, data["feature_id"]]
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
            row.append("{jbrowse_url}&assembly={assembly}&loc={chr}:{loc_from}..{loc_to}&tracks={tracks}&highlight={chr}:{hloc_from}..{hloc_to}".format(assembly=assembly, jbrowse_url=config.jbrowse2_url.format(treatment=treatment), treatment=treatment, chr=data["chr"], loc_from=loc_from, loc_to=loc_to, tracks=tracks, hloc_from=f_from, hloc_to=f_to))
            if feature_name=="junctions":
                row.append(database[data["feature_id"]]["annotated"]) 
                row.append(database[data["feature_id"]]["donor_anchor_id"])
                row.append(database[data["feature_id"]]["acceptor_anchor_id"])
            if feature_name=="genes":
                row += [smart_number_format(float(data["logFC"])), smart_number_format(float(data["F"])), smart_number_format(float(data["PValue"])), smart_number_format(float(data["FDR"]))]
                row += [smart_number_format(float(data["logFC"])* -np.log10(float(data["PValue"])))] # add pi value
            else:
                row += [smart_number_format(float(data["logFC"])), smart_number_format(float(data["exon.F"])), smart_number_format(float(data["P.Value"])), smart_number_format(float(data["FDR"]))]
                row += [smart_number_format(float(data["logFC"])* -np.log10(float(data["P.Value"])))] # add pi value
            results_all.append(row)
            if float(data["FDR"])<=config.edgeR_FDR_thr:
                results.append(row)
            r = f.readline()
        f.close()
        count += 1
        print(f"{module_name} {fname} {count}/{len(rfiles)}")

    headers = {}
    default_header = ["result_id", "comparison", "treatment", "feature_id", "chr", "strand", "feature_start", "feature_stop", "feature_len"]
    headers["genes"] = default_header + ["gene_id", "gene_name", "jbrowse_loc", "jbrowse_url", "logFC", "exon.F", "p_value", "fdr", "pi_value"]
    headers["junctions"] = default_header + ["junction_first_exon", "gene_id", "gene_name", "jbrowse_loc", "jbrowse_url", "annotated", "donor_anchor_id", "acceptor_anchor_id", "logFC", "exon.F", "p_value", "fdr", "pi_value"]
    headers["exons"] = default_header + ["gene_id", "gene_name", "jbrowse_loc", "jbrowse_url", "logFC", "exon.F", "p_value", "fdr", "pi_value"]
    headers["donor_anchors"] = default_header + ["gene_id", "gene_name", "jbrowse_loc", "jbrowse_url", "logFC", "exon.F", "p_value", "fdr", "pi_value"]
    headers["acceptor_anchors"] = default_header + ["gene_id", "gene_name", "jbrowse_loc", "jbrowse_url", "logFC", "exon.F", "p_value", "fdr", "pi_value"]

    # sort by FDR
    FDR_index = headers[feature_name].index("fdr")-1 # -1 because of result_id
    results_all.sort(key = lambda el: float(el[FDR_index]))
    results.sort(key = lambda el: float(el[FDR_index]))

    f = gzip.open(f"results/edgeR/{feature_name}_results_fdr005.tab.gz", "wt")
    f.write("\t".join(headers[feature_name]) + "\n")
    f_all = gzip.open(f"results/edgeR/{feature_name}_results_complete.tab.gz", "wt")
    f_all.write("\t".join(headers[feature_name]) + "\n")
    if len(results_all)>0:
        assert(len(headers[feature_name])==len(results_all[0])+1) # header and data columns must match, +1 for result_id
    result_id = 1
    for row in results:
        f.write("\t".join([f"r{result_id}"] + [str(el) for el in row]) + "\n")
        result_id += 1
    f.close()
    result_id = 1
    for row in results_all:
        f_all.write("\t".join([f"ra{result_id}"] + [str(el) for el in row]) + "\n")
        result_id += 1
    f_all.close()

def dexseq_feature(feature_name, version=""):
    # feature_name = "genes", "exons", "junctions", "donor_anchors", "acceptor_anchors"
    fname = f"results/results_dexseq{version}_{feature_name}.tab.gz"
    print(f"{module_name} generating dexseq results file={fname}, considering all results with FDR<={config.dexseq_FDR_thr}")
    print()
    comparisons = {}

    database = {}
    # read junctions reference
    if feature_name=="junctions":
        features.load_genes() # for loading promoters
        f = gzip.open("reference/junctions.tab.gz", "rt")
        header = f.readline().replace("\r", "").replace("\n", "").split("\t")
        r = f.readline()
        while r:
            r = r.replace("\r", "").replace("\n", "").split("\t")
            data = dict(zip(header, r))
            database[data["junction_id"]] = data
            r = f.readline()
        f.close()        

    samples = []
    f = open("samples.tab", "rt")
    r = f.readline()
    while r.startswith("#"):
        r = f.readline()
    header = r.replace("\r", "").replace("\n", "").split("\t")
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
        control_samples = []
        test_samples = []
        r = r.replace("\r", "").replace("\n", "").split("\t")
        comparison = r[0]
        for el in r[1].split(","):
            control_samples.append(el)
        for el in r[2].split(","):
            test_samples.append(el)
        comparisons[comparison] = (control_samples, test_samples)
        r = f.readline()
    f.close()

    if feature_name=="genes":
        rfiles = glob.glob(f"results/dexseq/{feature_name}/*difffeature.tab.gz")
    else:
        rfiles = glob.glob(f"results/dexseq/{feature_name}/*altsplice.tab.gz")
    results = [] # only results with FDR<config.dexseq_FDR_thr, stored in results_dexseq_{feature}.tab
    results_all = [] # all results, stored in results_dexseq_{feature}_all.tab
    count = 0
    for fname in rfiles:
        bname = os.path.basename(fname)
        if feature_name=="genes":
            comparison = "_".join(bname.split("_difffeature.tab.gz")[:-1])
        else:
            comparison = "_".join(bname.split("_altsplice.tab.gz")[:-1])
        treatment = bname.split("_"+config.control_name)[0]
        f = gzip.open(fname, "rt")
        header = f.readline().replace("\n", "").replace("\r", "").split("\t")
        r = f.readline()
        while r:
            r = r.replace("\n", "").replace("\r", "").split("\t")
            data = dict(zip(header, r))
            tracks_fixed = ["annotation_track"] # New JBrowse2
            assembly = "GenomeSequence"
            if config.jbrowse2_url.startswith("https://genomebrowser"):
                assembly = {"homo_sapiens":"hg38", "mus_musculus":"mm39"}[config.species]
                tracks_fixed = ["gene-refseq", "transcript-refseq", "transcript-ensembl"]
            tracks_control, tracks_test = comparisons[comparison] 

            if config.jbrowse2_url.startswith("https://genomebrowser"):
                tracks_control = [track.lstrip("0")+'_bw' for track in tracks_control]
                tracks_test = [track.lstrip("0")+'_bw' for track in tracks_test]
            else:
                tracks_control = [track+'_bw' for track in tracks_control] # New JBrowse2 --> bigwig files are called by id_bw
                tracks_test = [track+'_bw' for track in tracks_test] # New JBrowse2 --> bigwig files are called by id_bw

            tracks = ",".join(tracks_test[:3]+tracks_control[:3]+tracks_fixed)
            try:
                for r1, r2 in config.track_name_replace:
                    tracks = tracks.replace(r1, r2)
            except:
                pass
            f_from=int(data["feature_start"])
            f_to=int(data["feature_stop"])
            delta = round(abs(f_from-f_to)*0.15) # 15% surrounding
            loc_from=int(data["feature_start"])-delta
            loc_to=int(data["feature_stop"])+delta
            row = [comparison, treatment, data["feature_id"]]
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
            row.append("{jbrowse_url}&assembly={assembly}&loc={chr}:{loc_from}..{loc_to}&tracks={tracks}&highlight={chr}:{hloc_from}..{hloc_to}".format(assembly=assembly, jbrowse_url=config.jbrowse2_url.format(treatment=treatment), treatment=treatment, chr=data["chr"], loc_from=loc_from, loc_to=loc_to, tracks=tracks, hloc_from=f_from, hloc_to=f_to))
            if feature_name=="junctions":
                row.append(database[data["feature_id"]]["annotated"]) 
                row.append(database[data["feature_id"]]["donor_anchor_id"])
                row.append(database[data["feature_id"]]["acceptor_anchor_id"])
            if feature_name=="genes":
                row += [smart_number_format(float(data["logFC"])), smart_number_format(float(data["F"])), smart_number_format(float(data["PValue"])), smart_number_format(float(data["FDR"]))]
                row += [smart_number_format(float(data["logFC"])* -np.log10(float(data["PValue"])))] # add pi value
            else:
                row += [smart_number_format(float(data["logFC"])), smart_number_format(float(data["exon.F"])), smart_number_format(float(data["P.Value"])), smart_number_format(float(data["FDR"]))]
                row += [smart_number_format(float(data["logFC"])* -np.log10(float(data["P.Value"])))] # add pi value
            results_all.append(row)
            if float(data["FDR"])<=config.dexseq_FDR_thr:
                results.append(row)
            r = f.readline()
        f.close()
        count += 1
        print(f"{module_name} {fname} {count}/{len(rfiles)}")

    headers = {}
    default_header = ["result_id", "comparison", "treatment", "feature_id", "chr", "strand", "feature_start", "feature_stop", "feature_len"]
    headers["genes"] = default_header + ["gene_id", "gene_name", "jbrowse_loc", "jbrowse_url", "logFC", "exon.F", "p_value", "fdr", "pi_value"]
    headers["junctions"] = default_header + ["junction_first_exon", "gene_id", "gene_name", "jbrowse_loc", "jbrowse_url", "annotated", "donor_anchor_id", "acceptor_anchor_id", "logFC", "exon.F", "p_value", "fdr", "pi_value"]
    headers["exons"] = default_header + ["gene_id", "gene_name", "jbrowse_loc", "jbrowse_url", "logFC", "exon.F", "p_value", "fdr", "pi_value"]
    headers["donor_anchors"] = default_header + ["gene_id", "gene_name", "jbrowse_loc", "jbrowse_url", "logFC", "exon.F", "p_value", "fdr", "pi_value"]
    headers["acceptor_anchors"] = default_header + ["gene_id", "gene_name", "jbrowse_loc", "jbrowse_url", "logFC", "exon.F", "p_value", "fdr", "pi_value"]

    # sort by FDR
    FDR_index = headers[feature_name].index("fdr")-1 # -1 because of result_id
    results_all.sort(key = lambda el: float(el[FDR_index]))
    results.sort(key = lambda el: float(el[FDR_index]))

    f = gzip.open(f"results/dexseq/{feature_name}_results_fdr005.tab.gz", "wt")
    f.write("\t".join(headers[feature_name]) + "\n")
    f_all = gzip.open(f"results/dexseq/{feature_name}_results_complete.tab.gz", "wt")
    f_all.write("\t".join(headers[feature_name]) + "\n")
    if len(results_all)>0:
        assert(len(headers[feature_name])==len(results_all[0])+1) # header and data columns must match, +1 for result_id
    result_id = 1
    for row in results:
        f.write("\t".join([f"r{result_id}"] + [str(el) for el in row]) + "\n")
        result_id += 1
    f.close()
    result_id = 1
    for row in results_all:
        f_all.write("\t".join([f"ra{result_id}"] + [str(el) for el in row]) + "\n")
        result_id += 1
    f_all.close()
