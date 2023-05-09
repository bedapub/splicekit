import os
import sys
import splicekit.core.annotation as annotation
import splicekit.config as config

def toint(temp):
    try:
        temp = str(int(temp))
        return temp
    except:
        return temp

def same_junction(data1, data2):
    if data1["gene_id"]==data2["gene_id"] and data1["feature_start"]==data2["feature_start"] and data1["feature_stop"]==data2["feature_stop"]:
        return True
    return False

def process():
    annotation.read_comparisons()
    results = {}
    # populate result list
    f = open("results/results_edgeR_junctions.tab", "rt")
    header = f.readline().replace("\r", "").replace("\n", "").split("\t")
    r = f.readline()
    while r:
        r = r.replace("\r", "").replace("\n", "").split("\t")
        data = dict(zip(header, r))
        results[data["result_id"]] = data
        r = f.readline()
    f.close()

    # update results table
    count = 1
    for id, data in results.items():
        if count%100==0:
            print("[promisc] processed {a}/{b} ({c}% done)".format(a=count, b=len(results), c="%.2f" % (count/len(results)*100)))
        count += 1
        for rid, rdata in results.items():
            if same_junction(data, rdata):
                results[id][rdata["comparison"]+"_junction"] = 1

    reported_junction = {}
    fout_junction = open("results/results_edgeR_junctions_promisc.tab", "wt")
    header = ["gene_name", "gene_id", "chr", "strand", "feature_start", "feature_stop", "row_sum"]
    for (comparison, _, _, _, _) in annotation.comparisons:
        header.append(comparison)
    header.append("jbrowse_url")
    fout_junction.write("\t".join(header) + "\n")
    temp_junction = []
    for id, data in results.items():
        row_junction = []
        for hid in header[:6]:
            row_junction.append(data.get(hid, ""))
        for hid in header[7:]:
            row_junction.append(data.get(hid+"_junction", ""))
        row_sum_junction = sum([1 for el in row_junction if el==1])
        row_junction.insert(6, row_sum_junction)
        if reported_junction.get(tuple(row_junction[:6]), None)==None:
            temp_junction.append((row_sum_junction, row_junction))
            reported_junction[tuple(row_junction[:6])] = 1
    # junction level
    temp_junction.sort(reverse=True)
    for row_sum, row in temp_junction:
        fout_junction.write("\t".join([str(el) for el in row]) + "\n")
    fout_junction.close()
