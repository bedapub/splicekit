import os
import sys
import gzip
import itertools
import splicekit

module_desc = "splicekit | june |"

comparisons = {}
jbrowse_urls = {}

skipped_criteria = []
skipped_criteria.append('j1["feature_start"]==j2["feature_start"]')
skipped_criteria.append('j2["feature_stop"]<j1["feature_stop"]')
skipped_criteria.append('j3["feature_stop"]==j1["feature_stop"]')
skipped_criteria.append('j2["feature_stop"]<j3["feature_start"]')
skipped_criteria.append('float(j2["logFC"])*float(j3["logFC"])>0 and float(j1["logFC"])*float(j2["logFC"])<0')

mutual_criteria = []
mutual_criteria.append('j1["feature_start"]==j3["feature_start"]')
mutual_criteria.append('j2["feature_stop"]==j4["feature_stop"]')
mutual_criteria.append('j1["feature_stop"]<j2["feature_start"]')
mutual_criteria.append('j3["feature_stop"]<j4["feature_start"]')
mutual_criteria.append('j2["feature_stop"]==j4["feature_stop"]')
mutual_criteria.append('j3["feature_stop"]>j2["feature_start"]')
mutual_criteria.append('(float(j1["logFC"])>0 and float(j2["logFC"])>0 and float(j3["logFC"])<0 and float(j4["logFC"])<0) or (float(j1["logFC"])<0 and float(j2["logFC"])<0 and float(j3["logFC"])>0 and float(j4["logFC"])>0)')

def load_junctions():
    print(f"{module_desc} loading junctions")
    f = gzip.open("results/edgeR/junctions_results_fdr005.tab.gz", "rt")
    header = f.readline().replace("\r", "").replace("\n", "").split("\t")
    r = f.readline()
    while r:
        r = r.replace("\r", "").replace("\n", "").split("\t")
        data = dict(zip(header, r))
        junctions = comparisons.get(data["comparison"], {})
        temp = junctions.get(data["gene_id"], [])
        temp2 = data["feature_id"].split("_")
        chr, start, stop = "_".join(temp2[:-2]), temp2[-2], temp2[-1]
        chr, strand = chr[:-1], chr[-1]
        start, stop = int(start), int(stop)
        jbrowse_urls[data["feature_id"]] = data["jbrowse_url"]
        # lets expand data a little bit to make things easier
        data["feature_start"] = start
        data["feature_stop"] = stop
        temp.append(data)
        junctions[data["gene_id"]] = temp
        comparisons[data["comparison"]] = junctions
        r = f.readline()
    f.close()

def compute_results():
    print(f"{module_desc} computing results")
    cryptic_exons = {}
    for comp_id, junctions in comparisons.items():
        for gene_id, data in junctions.items():

            # exon skipping
            if len(data)>=3:
                combs = list(itertools.permutations(data, 3))
                for j1, j2, j3 in combs:
                    assert(j1["gene_id"]==j2["gene_id"]==j3["gene_id"])
                    add_event = True
                    for criteria in skipped_criteria:
                        if not eval(criteria):
                            add_event = False
                            break
                    if add_event:
                        record_id = f"{comp_id}_{j1['chr']}{j1['strand']}_{j2['feature_stop']}_{j3['feature_start']}"
                        exon_id = f"{j1['chr']}{j1['strand']}_{j2['feature_stop']}_{j3['feature_start']}"
                        if j1["strand"]=="+":
                            exon_annotation = f'{j2["annotated"][1]}{j3["annotated"][0]}'
                        else:
                            exon_annotation = f'{j3["annotated"][1]}{j2["annotated"][0]}'
                        tracks = j1["jbrowse_url"].split("tracks=")[1].split("&highlight=")[0]
                        exon_jbrowse = f"{splicekit.config.jbrowse2_url}&assembly=hg38&loc={j1['chr']}:{j1['feature_start']}..{j1['feature_stop']}&tracks={tracks}&highlight={j1['chr']}:{j2['feature_stop']}..{j3['feature_start']}"
                        exon_data = {}
                        exon_data["june_type"] = "skipping"
                        exon_data["comp_id"] = comp_id
                        exon_data["exon_id"] = exon_id
                        exon_data["gene_id"] = j1["gene_id"]
                        exon_data["annotation"] = exon_annotation
                        exon_data["j1"] = j1["feature_id"]
                        exon_data["j2"] = j2["feature_id"]
                        exon_data["j3"] = j3["feature_id"]
                        exon_data["jbrowse"] = exon_jbrowse
                        exon_data["j1_fdr"] = float(j1['fdr'])
                        exon_data["j1_logFC"] = float(j1['logFC'])
                        delta_logFC = float(j1['logFC'])-(float(j2['logFC'])+float(j3['logFC']))
                        exon_data["delta_logFC"] = delta_logFC
                        cryptic_exons[record_id] = exon_data

            # mutually exclusive exons
            if len(data)>=4:
                combs = list(itertools.permutations(data, 4))
                for j1, j2, j3, j4 in combs:
                    assert(j1["gene_id"]==j2["gene_id"]==j3["gene_id"]==j4["gene_id"])
                    add_event = True
                    for criteria in mutual_criteria:
                        if not eval(criteria):
                            add_event = False
                            break
                    if add_event:
                        record_id = f"{comp_id}_{j1['chr']}{j1['strand']}_{j1['feature_start']}_{j2['feature_stop']}"
                        exon_id_1 = f"{j1['chr']}{j1['strand']}_{j1['feature_stop']}_{j2['feature_start']}"
                        exon_id_2 = f"{j1['chr']}{j1['strand']}_{j3['feature_stop']}_{j4['feature_start']}"
                        tracks = j1["jbrowse_url"].split("tracks=")[1].split("&highlight=")[0]
                        exon_jbrowse = f"{splicekit.config.jbrowse2_url}&assembly=hg38&loc={j1['chr']}:{j1['feature_start']}..{j2['feature_stop']}&tracks={tracks}&highlight={j1['chr']}:{j1['feature_stop']}..{j2['feature_start']}"
                        exon_data = {}
                        exon_data["june_type"] = "mutual"
                        exon_data["comp_id"] = comp_id
                        exon_data["exon_id"] = f"{exon_id_1},{exon_id_2}"
                        exon_data["gene_id"] = j1["gene_id"]
                        if j1["strand"]=="+":
                            exon_annotation = f'{j1["annotated"][1]}{j2["annotated"][0]}|{j3["annotated"][1]}{j4["annotated"][0]}'
                        else:
                            exon_annotation = f'{j4["annotated"][1]}{j3["annotated"][0]}|{j2["annotated"][1]}{j1["annotated"][0]}'
                        exon_data["annotation"] = exon_annotation
                        exon_data["j1"] = j1["feature_id"]
                        exon_data["j2"] = j2["feature_id"]
                        exon_data["j3"] = j3["feature_id"]
                        exon_data["jbrowse"] = exon_jbrowse
                        exon_data["j1_fdr"] = float(j1['fdr'])
                        exon_data["j1_logFC"] = float(j1['logFC'])
                        delta_logFC = (float(j1['logFC'])+float(j2['logFC']))-(float(j3['logFC'])+float(j4['logFC']))
                        exon_data["delta_logFC"] = delta_logFC
                        cryptic_exons[record_id] = exon_data
    results = []
    for exon_id, exon_data in cryptic_exons.items():
        results.append((exon_data["delta_logFC"], exon_data))
    results = sorted(results, key=lambda x: x[0])
    fout = gzip.open("results/edgeR/june.tab.gz", "wt")
    header = ["comparison", "june_type", "exon_id", "delta_logFC", "gene_id", "j1", "exon_annotation", "exon_jbrowse"]
    fout.write("\t".join(header)+"\n")
    for exon_id, exon_data in results:
        fout.write(f"{exon_data['comp_id']}\t{exon_data['june_type']}\t{exon_data['exon_id']}\t{exon_data['delta_logFC']}\t{exon_data['gene_id']}\t{exon_data['j1']}\t{exon_data['annotation']}\t{exon_data['jbrowse']}\n")
    fout.close()

    # GPT ready results
    template_description = "This file contains splicing events across genes for diverse comparisons. Each block of lines describes a splicing events with its parameters and values. Each block starts with [splicing_event]."
    template="""
[splicing_event]
comparison={comp_id}
gene_id={gene_id}
exon_id={exon_id}
exon_annotation={exon_annotation}
event_type={june_type}
logFC={delta_logFC}
jbrowse={exon_jbrowse}
description=The comparison {comp_id} has a splicing event in the gene {gene_id} of type {june_type}, with logFC of {delta_logFC}. The affected exon(s) are {exon_id}. The annotation of the exon is {exon_annotation}.
    """
    fout = open("results/edgeR/june.gpt.txt", "wt")
    fout.write(template_description+"\n")
    for exon_id, exon_data in results:
        text = template.format(gene_id=exon_data['gene_id'], comp_id=exon_data['comp_id'], june_type=exon_data['june_type'], exon_id=exon_data['exon_id'], delta_logFC=exon_data['delta_logFC'], exon_annotation=exon_data['annotation'], exon_jbrowse=exon_data['jbrowse'])
        fout.write(text)
    fout.close()

def process():
    load_junctions()
    compute_results()