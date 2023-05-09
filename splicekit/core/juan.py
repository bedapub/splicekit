import os
import splicekit
import glob
import pandas as pd

# append anchor edgeR results (donor+acceptor) to results_edgeR_junctions
def append_results():
    def read_database(anchor_type):
        temp = {}
        for comparison in splicekit.core.annotation.comparisons:
            comp_id = comparison[0]
            print("[juan] reading:", f"results/results_edgeR_{anchor_type}/{comp_id}_altsplice.tab")
            if not os.path.exists(f"results/results_edgeR_{anchor_type}/{comp_id}_altsplice.tab"):
                print("[juan] warning, do such file:", f"results/results_edgeR_{anchor_type}/{comp_id}_altsplice.tab")
                continue
            f = open(f"results/results_edgeR_{anchor_type}/{comp_id}_altsplice.tab", "rt")
            header = f.readline().replace("\r", "").replace("\n", "").split("\t")
            r = f.readline()
            while r:
                r = r.replace("\r", "").replace("\n", "").split("\t")
                data = dict(zip(header, r))
                key = comp_id + "_" + data["feature_id"]
                temp[key] = data
                r = f.readline()
            f.close()
        return temp
    
    database_donor_anchors = read_database("donor_anchors")
    database_acceptor_anchors = read_database("acceptor_anchors")

    for fname in ["results/results_edgeR_junctions", "results/results_edgeR_junctions_all"]:
        f = open(f"{fname}.tab", "rt")
        fout = open(f"{fname}2.tab", "wt")
        header = f.readline().replace("\r", "").replace("\n", "").split("\t")
        # already appended anchor info to the results file?
        if "donor_anchor_fdr" in header:
            continue
        header_new = header
        header_new.append("donor_anchor_fdr")
        header_new.append("donor_anchor_logfc")
        header_new.append("acceptor_anchor_fdr")
        header_new.append("acceptor_anchor_logfc")
        fout.write("\t".join(header_new) + "\n")
        r = f.readline()
        while r:
            r = r.replace("\r", "").replace("\n", "").split("\t")
            data = dict(zip(header, r))
            comp_id = data["comparison"]
            donor_anchor_data = database_donor_anchors.get(comp_id+"_"+data["donor_anchor_id"], None)
            acceptor_anchor_data = database_acceptor_anchors.get(comp_id+"_"+data["acceptor_anchor_id"], None)
            if donor_anchor_data!=None:
                r.append(donor_anchor_data["FDR"])
                r.append(donor_anchor_data["logFC"])
            else:
                r.append("")
                r.append("")
            if acceptor_anchor_data!=None:
                r.append(acceptor_anchor_data["FDR"])
                r.append(acceptor_anchor_data["logFC"])
            else:
                r.append("")
                r.append("")
            fout.write("\t".join(r) + "\n")
            r = f.readline()
        f.close()
        fout.close()
        os.system(f"mv {fname}2.tab {fname}.tab")
