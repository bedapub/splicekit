import os
import splicekit
import glob
import gzip
import pandas as pd

module_name = "splicekit | juan |"

# append anchor edgeR results (donor+acceptor) to results_edgeR_junctions
def append_results():
    def read_database(anchor_type):
        temp = {}
        for comparison in splicekit.core.annotation.comparisons:
            comp_id = comparison[0]
            altsplice_fname = f"results/edgeR/{anchor_type}/{comp_id}_altsplice.tab.gz"
            print(f"{module_name} reading {altsplice_fname}")
            if not os.path.exists(altsplice_fname):
                print(f"{module_name} warning, missing file {altsplice_fname}")
                continue
            f = gzip.open(altsplice_fname, "rt")
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

    for fname in ["results/edgeR/junctions_results_complete.tab.gz", "results/edgeR/junctions_results_fdr005.tab.gz"]:
        f = gzip.open(fname, "rt")
        header = f.readline().replace("\r", "").replace("\n", "").split("\t")
        header_new = header.copy()
        for el in ["donor_anchor_fdr", "donor_anchor_logfc", "acceptor_anchor_fdr", "acceptor_anchor_logfc"]:
            if el not in header:
                header_new.append(el)
        fout = gzip.open(f"{fname}.temp", "wt")
        fout.write("\t".join(header_new) + "\n")
        r = f.readline()
        while r:
            r = r.replace("\r", "").replace("\n", "").split("\t")
            data_new = dict(zip(header_new, r))
            comp_id = data_new["comparison"]
            data_new["donor_anchor_fdr"] = database_donor_anchors.get(comp_id+"_"+data_new["donor_anchor_id"], {}).get("FDR", "")
            data_new["donor_anchor_logfc"] = database_donor_anchors.get(comp_id+"_"+data_new["donor_anchor_id"], {}).get("logFC", "")
            data_new["acceptor_anchor_fdr"] = database_acceptor_anchors.get(comp_id+"_"+data_new["acceptor_anchor_id"], {}).get("FDR", "")
            data_new["acceptor_anchor_logfc"] = database_acceptor_anchors.get(comp_id+"_"+data_new["acceptor_anchor_id"], {}).get("logFC", "")
            fout.write("\t".join([str(data_new[x]) for x in header_new]) + "\n")
            r = f.readline()
        f.close()
        fout.close()
        os.system(f"mv {fname}.temp {fname}")
