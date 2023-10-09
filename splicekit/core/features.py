"""
# Description
Generates juncion, exon, genes, anchor count tables.
"""

import os
import sys
import math
import gzip
import pybio
import splicekit.config as config
import splicekit.core.annotation as annotation
import scipy
import scipy.stats

# load genes with exons, junctions and anchors
def load_genes():
    annotation.junctions_genes = {}
    annotation.donor_anchors_genes = {}
    annotation.acceptor_anchors_genes = {}
    annotation.exons_genes = {}
    transcript_exons = {}
    print("[features] Reading annotation:", config.gtf_path)
    f = gzip.open(config.gtf_path, "rt")
    r = f.readline()
    while r:
        if r.startswith("#"):
            r = f.readline()
            continue
        r = r.replace("\r", "").replace("\n", "").split("\t")
        if r[2]!="exon":
            r = f.readline()
            continue
        chr, pos_start, pos_stop, strand, temp_atts = r[0], r[3], r[4], r[6], r[8]
        atts = {}
        pos_start, pos_stop = int(pos_start)-1, int(pos_stop)-1
        temp_atts = temp_atts.split("; ")
        for el in temp_atts:
            el = el.split(" ")
            if len(el)!=2:
                continue
            atts[el[0]] = el[1].replace("\"", "")
        gene_id = atts["gene_id"]
        gene = annotation.genes.get(gene_id, {})
        gene["chr"] = chr
        gene["strand"] = strand
        gene["gene_name"] = atts.get("gene_name", "")

        # because we have a feature called "genes" (and exons, donor_anchors, acceptor_anchors and junctions)
        # this stores gene counts
        gene_info = gene.get("genes", {})
        gene_info["start"] = min(gene_info.get("start", math.inf), pos_start)
        gene_info["stop"] = max(gene_info.get("start", 0), pos_stop)
        gene["genes"] = gene_info

        # exons
        exon_id = gene["chr"]+gene["strand"]+"_" + str(pos_start) + "_" + str(pos_stop)
        annotation.exons_genes[exon_id] = gene_id
        exons = gene.get("exons", {})
        exons[exon_id] = {"start":pos_start, "stop":pos_stop} # sample counts dict
        gene["exons"] = exons

        exon_first_pos = pos_start if strand=="+" else pos_stop
        exon_last_pos = pos_stop if strand=="+" else pos_start
        transcript_exons.setdefault((atts["transcript_id"], chr, strand, gene_id), []).append((exon_first_pos, exon_last_pos))

        annotation.genes[gene_id] = gene
        r = f.readline()
    f.close()

    def identify_exons(transcript_exons, annotation_data, exon_index):
        for (transcript_id, chr, strand, gene_id), transcript_list in transcript_exons.items():
            # sort exons by start position
            transcript_list = sorted(transcript_list, key = lambda x: (x[0])) if strand=="+" else sorted(transcript_list, key = lambda x: (-x[0]))
            # keep only first exons, even if there are >1 (same start position, but diff stop position)
            transcript_list = [(x[0], x[1]) for x in transcript_list if x[0]==transcript_list[exon_index][0]]
            for (exon_start, exon_stop) in transcript_list:
                store_data = True
                # store exon last pos, since we find this exon by first nucleotide of junction (last nucleotide of first exon)
                current_data = annotation_data.get((chr, strand, exon_stop), None)
                if current_data!=None: # already an exon stored at this "last" position
                    current_transcript_id, current_gene_id, current_exon_start, current_exon_stop = current_data
                    current_len = abs(current_exon_start - current_exon_stop + 1)
                    new_len = abs(exon_start - exon_stop + 1)
                    if new_len>current_len:
                        store_data = False
                if store_data:
                    annotation_data[chr, strand, exon_stop] = (transcript_id, gene_id, exon_start, exon_stop)

    # save last nucleotide of first/second exon of each transcript)
    print("splicekit.features: saving last nucleotide of first/second exon of each transcript")
    annotation.first_exons = {}
    identify_exons(transcript_exons, annotation.first_exons, 0)
    #annotation.second_exons = {}
    #identify_exons(transcript_exons, annotation.second_exons, 1)
    
    print("splicekit.features: reading junctions and anchors annotation from: reference/junctions.tab")
    f = open("reference/junctions.tab", "rt")
    header = f.readline().replace("\r", "").replace("\n", "").split("\t")
    r = f.readline()
    while r:
        r = r.replace("\r", "").replace("\n", "").split("\t")
        data = dict(zip(header, r))
        gene_id = data["gene_id"]
        gene = annotation.genes.get(gene_id, {})
        gene["chr"] = data["chr"]
        gene["strand"] = data["strand"]
        gene["gene_name"] = data["gene_name"]
        # junctions
        annotation.junctions_genes[data["junction_id"]] = gene_id
        junctions = gene.get("junctions", {})
        junction_start, junction_stop = int(data["junction_id"].split("_")[-2]), int(data["junction_id"].split("_")[-1])
        junctions[data["junction_id"]] = {"start":junction_start, "stop":junction_stop} # + sample counts later
        gene["junctions"] = junctions
        # anchors
        for anchor_type in ["donor", "acceptor"]:
            if anchor_type=="donor":
                annotation.donor_anchors_genes[data[f"{anchor_type}_anchor_id"]] = gene_id
            if anchor_type=="acceptor":
                annotation.acceptor_anchors_genes[data[f"{anchor_type}_anchor_id"]] = gene_id
            anchors = gene.get(f"{anchor_type}_anchors", {})
            anchor_start, anchor_stop = int(data[f"{anchor_type}_anchor_id"].split("_")[-2]), int(data[f"{anchor_type}_anchor_id"].split("_")[-1])
            anchors[data[f"{anchor_type}_anchor_id"]] = {"start":anchor_start, "stop":anchor_stop} # + sample counts later
            gene[f"{anchor_type}_anchors"] = anchors
        annotation.genes[gene_id] = gene
        r = f.readline()
    f.close()

# simply reads in the data from the sample_junctions_data files
def read_junctions():
    novel_ok = 0
    novel_fail = 0
    count = 0
    count_all = len(annotation.samples)
    for sample_id in annotation.samples:
        f = open(f"data/sample_junctions_data/sample_{sample_id}.tab", "rt")
        r = f.readline() # header
        r = f.readline()
        while r:
            r = r.replace("\r", "").replace("\n", "").split("\t")
            junction_id = r[0]
            coords = junction_id.split('_')
            start = int(coords[-2])
            stop = int(coords[-1])
            strand = coords[-3][-1]
            chr = '_'.join(coords[:-2])[:-1]            
            raw_count = int(r[-1])
            gene_id = annotation.junctions_genes[junction_id]
            gene = annotation.genes.get(gene_id, {})
            junctions = gene.get("junctions", {})
            junctions.setdefault(junction_id, {})[sample_id] = raw_count # sample counts dict
            gene["junctions"] = junctions
            annotation.genes[gene_id] = gene
            r = f.readline()
        count += 1
        print("[features] Junction read OK, sample={sample_id}, from=data/sample_junctions_data/sample_{sample_id}.tab, ({a}/{b}, {c}% done)".format(a=count, b=count_all, c="%.2f" % (count/count_all*100), sample_id=sample_id))
        f.close()

def read_anchors(anchor_type):
    """
    # Description
    Read anchors sample raw counts and populate the annotation.genes ["anchors"] field.
    The input is read from the output of juan scripts (data/sample_anchors_data) producing anchor read counts.
    """
    count = 0
    count_all = len(annotation.samples)
    for sample_id in annotation.samples:
        f = open(f"data/sample_{anchor_type}_anchors_data/sample_{sample_id}.tab", "rt")
        r = f.readline() # header
        r = f.readline()
        while r:
            r = r.replace("\r", "").replace("\n", "").split("\t")
            anchor_id = r[0]
            coords = anchor_id.split('_')
            start = int(coords[-2])
            stop = int(coords[-1])
            strand = coords[-3][-1]
            chr = '_'.join(coords[:-2])[:-1]            
            raw_count = int(r[-1])
            if anchor_type=="donor":
                gene_id = annotation.donor_anchors_genes[anchor_id]
            if anchor_type=="acceptor":
                gene_id = annotation.acceptor_anchors_genes[anchor_id]
            gene = annotation.genes.get(gene_id, {})
            anchors = gene.get(f"{anchor_type}_anchors", {})
            anchors.setdefault(anchor_id, {})[sample_id] = raw_count # sample counts dict
            gene[f"{anchor_type}_anchors"] = anchors
            annotation.genes[gene_id] = gene
            r = f.readline()
        count += 1
        print("[features] Anchors read OK, sample={sample_id}, from=data/sample_{anchor_type}_anchors_data/sample_{sample_id}.tab, ({a}/{b}, {c}% done)".format(a=count, b=count_all, c="%.2f" % (count/count_all*100), sample_id=sample_id, anchor_type=anchor_type))
        f.close()

def read_exons():
    """
    # Description
    Read exons sample raw counts and populate the annotation.genes ["exons"] field.
    The input is read from data/sample_exons_data
    """
    count = 0
    count_all = len(annotation.samples)
    for sample_id in annotation.samples:
        f = open(f"data/sample_exons_data/sample_{sample_id}.tab", "rt")
        r = f.readline() # header
        r = f.readline()
        while r:
            r = r.replace("\r", "").replace("\n", "").split("\t")
            exon_id = r[0]
            coords = exon_id.split('_')
            start = int(coords[-2])
            stop = int(coords[-1])
            strand = coords[-3][-1]
            chr = '_'.join(coords[:-2])[:-1]            
            raw_count = int(r[-1])
            gene_id = annotation.exons_genes[exon_id]
            gene = annotation.genes.get(gene_id, {})
            exons = gene.get("exons", {})
            exons.setdefault(exon_id, {})[sample_id] = raw_count # sample counts dict
            gene["exons"] = exons
            annotation.genes[gene_id] = gene
            r = f.readline()
        count += 1
        print("[features] Exons read OK, sample={sample_id}, from=data/sample_exons_data/sample_{sample_id}.tab, ({a}/{b}, {c}% done)".format(a=count, b=count_all, c="%.2f" % (count/count_all*100), sample_id=sample_id))
        f.close()

def read_genes():
    """
    # Description
    Read gene sample raw counts and populate the annotation.genes ["genes"] field.
    The input is read from data/sample_genes_data
    """
    count = 0
    count_all = len(annotation.samples)
    for sample_id in annotation.samples:
        f = open(f"data/sample_genes_data/sample_{sample_id}.tab", "rt")
        r = f.readline() # header
        r = f.readline()
        while r:
            r = r.replace("\r", "").replace("\n", "").split("\t")
            gene_id = r[0]
            raw_count = int(r[-1])
            # to make it compatible with other features (exons, anchors), we have a "gene_id" at the gene_counts here
            gene = annotation.genes.get(gene_id, {})
            gene_counts = gene.get("genes", {})
            gene_counts.setdefault(gene_id, {})[sample_id] = raw_count # sample counts dict
            gene["genes"] = gene_counts
            annotation.genes[gene_id] = gene
            r = f.readline()
        count += 1
        print("[features] Genes read OK, sample={sample_id}, from=data/sample_genes_data/sample_{sample_id}.tab, ({a}/{b}, {c}% done)".format(a=count, b=count_all, c="%.2f" % (count/count_all*100), sample_id=sample_id))
        f.close()

def save_comps_feature_data(feature_type):
    """
    # Description
    Take files from data/samples_{feature_type}_data/*.tab and read in the counts
    Write the data/comparison_{feature_type}_data files

    # Parameters
    feature_type = exons, junctions, donor_anchors, acceptor_anchors

    # Structure and debug
    comp1 = [(36, 'DMSO', 'DMSO_rep1'), (48, 'DMSO', 'DMSO_rep1'), (60, 'DMSO', 'DMSO_rep1'), (72, 'DMSO', 'DMSO_rep1'), (84, 'DMSO', 'DMSO_rep1'), (96, 'DMSO', 'DMSO_rep1'), (132, 'DMSO', 'DMSO_rep1'), (144, 'DMSO', 'DMSO_rep1'), (156, 'DMSO', 'DMSO_rep1'), (168, 'DMSO', 'DMSO_rep1'), (180, 'DMSO', 'DMSO_rep1'), (192, 'DMSO', 'DMSO_rep1'), (228, 'DMSO', 'DMSO_rep2'), (240, 'DMSO', 'DMSO_rep2'), (252, 'DMSO', 'DMSO_rep2'), (264, 'DMSO', 'DMSO_rep2'), (276, 'DMSO', 'DMSO_rep2'), (288, 'DMSO', 'DMSO_rep2'), (324, 'DMSO', 'DMSO_rep2'), (336, 'DMSO', 'DMSO_rep2'), (348, 'DMSO', 'DMSO_rep2'), (360, 'DMSO', 'DMSO_rep2'), (372, 'DMSO', 'DMSO_rep2'), (384, 'DMSO', 'DMSO_rep2'), (420, 'DMSO', 'DMSO_rep3'), (432, 'DMSO', 'DMSO_rep3'), (444, 'DMSO', 'DMSO_rep3'), (456, 'DMSO', 'DMSO_rep3'), (468, 'DMSO', 'DMSO_rep3'), (480, 'DMSO', 'DMSO_rep3'), (516, 'DMSO', 'DMSO_rep3'), (528, 'DMSO', 'DMSO_rep3'), (540, 'DMSO', 'DMSO_rep3'), (552, 'DMSO', 'DMSO_rep3'), (564, 'DMSO', 'DMSO_rep3'), (576, 'DMSO', 'DMSO_rep3')]
    comp2 = [(1, 'RO7282701', 'RO7282701-000-001_rep1'), (193, 'RO7282701', 'RO7282701-000-001_rep2'), (385, 'RO7282701', 'RO7282701-000-001_rep3')]
    """
    count = 0
    count_all = len(annotation.comparisons)
    for (comp_name, comp1, comp2, _, _) in annotation.comparisons:
        fout = open("data/comparison_{feature}_data/{comp_name}.tab".format(feature=feature_type, comp_name=comp_name), "wt")
        header = ["gene_id", "gene_name", "chr", "strand", "feature_start", "feature_stop", "length", "feature_id"]
        for (sample_id, compound, rep, _) in comp2:
            header.append("{sample}_{compound}".format(sample=sample_id, compound=compound))
        for (sample_id, compound, rep, _) in comp1:
            header.append("{sample}_{compound}".format(sample=sample_id, compound=compound))
        header.append("sum_gene_test")
        header.append("sum_gene_control")
        header.append("sum_feature_test")
        header.append("sum_feature_control")
        header.append("test_pfi")
        header.append("control_pfi")
        header.append("delta_pfi")
        fout.write("\t".join(header)+"\n")
        for gene_id, gene_data in annotation.genes.items():
            if gene_data.get(feature_type, None)==None:
                continue
            sum_gene_control = 0
            sum_gene_test = 0
            for feature_id, feature_data in gene_data[feature_type].items():
                # propagate gene start stop from gene object to count object
                if feature_type=="genes" and feature_id in ["start", "stop"]:
                    continue
                if feature_type=="genes":
                    feature_data["start"] = gene_data["genes"]["start"]
                    feature_data["stop"] = gene_data["genes"]["stop"]
                for (sample_id, compound, rep, _) in comp2:
                    sum_gene_test += feature_data.get(sample_id, 0)
                for (sample_id, compound, rep, _) in comp1:
                    sum_gene_control += feature_data.get(sample_id, 0)
            for feature_id, feature_data in gene_data[feature_type].items():
                # propagate gene start stop from gene object to count object
                if feature_type=="genes" and feature_id in ["start", "stop"]:
                    continue
                if feature_type=="genes":
                    feature_data["start"] = gene_data["genes"]["start"]
                    feature_data["stop"] = gene_data["genes"]["stop"]
                row = [gene_id, "{symbol}\t{chr}\t{strand}\t{feature_start}\t{feature_stop}\t{length}".format(chr=gene_data["chr"], strand=gene_data["strand"], symbol=gene_data["gene_name"], feature_start=feature_data["start"], feature_stop=feature_data["stop"], length=feature_data["stop"]-feature_data["start"]+1)]
                row.append(feature_id)
                sum_feature_control = 0
                sum_feature_test = 0
                for (sample_id, compound, rep, _) in comp2:
                    row.append(feature_data.get(sample_id, 0))
                    sum_feature_test += feature_data.get(sample_id, 0)
                for (sample_id, compound, rep, _) in comp1:
                    row.append(feature_data.get(sample_id, 0))
                    sum_feature_control += feature_data.get(sample_id, 0)
                row.append(sum_gene_test)
                row.append(sum_gene_control)
                row.append(sum_feature_test)
                row.append(sum_feature_control)
                test_pfi = 100 * (sum_feature_test/(sum_gene_test+1))
                control_pfi = 100 * (sum_feature_control/(sum_gene_control+1))
                delta_pfi = test_pfi - control_pfi
                row.append("%.2f" % test_pfi)
                row.append("%.2f" % control_pfi)
                row.append("%.2f" % delta_pfi)
                fout.write("\t".join(str(el) for el in row)+"\n")
        fout.close()
        count += 1
        print("[features] Saved comparison_{feature}_data/{comp_name}.tab ({a}/{b}, {c}% done)".format(comp_name=comp_name, a=count, b=count_all, c="%.2f" % (count/count_all*100), feature=feature_type))
    return True

def add_psi_cluster():
    commands = []
    commands.append("export BSUB_QUIET=Y")
    for (comp_name, comp1, comp2, _, _) in annotation.comparisons:
        commands.append("bsub -q short -o /dev/null -e /dev/null -K python -c 'import splicekit; splicekit.core.features.add_psi(\"" + comp_name + "\")'")
    commands = "export BSUB_QUIET=Y; " + " & ".join(commands) + "; wait"
    os.system(commands)

def add_psi(comp_name):
    read_length = 150
    """
    # Description
    Adds delta_PSI to exon results tables
    Reference: https://currentprotocols.onlinelibrary.wiley.com/doi/full/10.1002/0471142905.hg1116s87

    IR = IR_exon / (exon_length + read_length - 1) # reads covering the exon
    ER = ER_exon / (read_length - 1) # reads skipping the exon
    PSI = IR / (IR+ER) * 100
    """

    def find_junctions(data, chr, strand, feature_start, feature_stop):
        sum_test = 0
        sum_control = 0
        candidates = data.get(chr, {}).get(strand, [])
        for (junction_start, junction_stop, junction_test, junction_control) in candidates:
            if junction_start < feature_start < feature_stop < junction_stop:
                sum_test += junction_test
                sum_control += junction_control
        return sum_test, sum_control

    fin_exons = open("data/comparison_exons_data/{comp_name}.tab".format(comp_name=comp_name), "rt")
    fin_junctions = open("data/comparison_junctions_data/{comp_name}.tab".format(comp_name=comp_name), "rt")
    data_junctions = {} # chr / strand + list of junctions
    # first read in the junctions
    print("[features] reading junctions from: data/comparison_junctions_data/{comp_name}.tab".format(comp_name=comp_name))
    header = fin_junctions.readline().replace("\r", "").replace("\n", "").split("\t")
    r = fin_junctions.readline()
    while r:
        r = r.replace("\r", "").replace("\n", "").split("\t")
        data = dict(zip(header, r))
        chr = data["chr"]
        strand = data["strand"]
        junction_start = int(data["feature_start"])
        junction_stop = int(data["feature_stop"])
        sum_feature_test = int(data["sum_feature_test"])
        sum_feature_control = int(data["sum_feature_control"])
        data_junctions.setdefault(chr, {}).setdefault(strand, []).append((junction_start, junction_stop, sum_feature_test, sum_feature_control))
        r = fin_junctions.readline()
    fin_junctions.close()

    fout_exons = open("data/comparison_exons_data/{comp_name}_new.tab".format(comp_name=comp_name), "wt")
    print("[features] processing exons from: data/comparison_exons_data/{comp_name}.tab".format(comp_name=comp_name))
    header = fin_exons.readline().replace("\r", "").replace("\n", "").split("\t")
    header_out = header.copy()
    if "test_PSI" not in header_out:
        header_out.append("test_PSI")
    if "control_PSI" not in header_out:
        header_out.append("control_PSI")
    if "delta_PSI" not in header_out:
        header_out.append("delta_PSI")
    fout_exons.write("\t".join(header_out) + "\n")
    r = fin_exons.readline()
    while r:
        r = r.replace("\r", "").replace("\n", "").split("\t")
        data = dict(zip(header, r))
        exon_length = int(data["length"])
        exon_start = int(data["feature_start"])
        exon_stop = int(data["feature_stop"])
        chr = data["chr"]
        strand = data["strand"]
        sum_exon_test = int(data["sum_feature_test"])
        sum_exon_control = int(data["sum_feature_control"])
        sum_junction_test, sum_junction_control = find_junctions(data_junctions, chr, strand, exon_start, exon_stop)
        test_IR = sum_exon_test/float(exon_length + read_length - 1)
        test_ER = sum_junction_test/float(read_length - 1)
        control_IR = sum_exon_control/float(exon_length + read_length - 1)
        control_ER = sum_junction_control/float(read_length - 1)
        try:
            test_PSI = test_IR/(test_IR + test_ER) * 100
        except:
            test_PSI = 0
            #print("ZERO", chr, strand, exon_length, sum_exon_test, sum_exon_control, sum_junction_test, sum_junction_control, test_IR, test_ER)
        try:
            control_PSI = control_IR/(control_IR + control_ER) * 100
        except:
            control_PSI = 0
            #print("ZERO", chr, strand, exon_length, sum_exon_test, sum_exon_control, sum_junction_test, sum_junction_control, test_IR, test_ER)
        delta_PSI = test_PSI - control_PSI
        #if abs(delta_PSI>10):
        #    print(chr, strand, exon_length, sum_exon_test, sum_exon_control, sum_junction_test, sum_junction_control, test_PSI, control_PSI, delta_PSI)
        data_out = dict(zip(header_out, r))
        data_out["test_PSI"] = "%.2f" % test_PSI
        data_out["control_PSI"] = "%.2f" % control_PSI
        data_out["delta_PSI"] = "%.2f" % delta_PSI
        fout_exons.write("\t".join(str(data_out[h]) for h in header_out) + "\n")
        r = fin_exons.readline()
    fin_exons.close()
    fout_exons.close()
    os.system("mv {fout_exons_new} {fout_exons}".format(fout_exons_new="data/comparison_exons_data/{comp_name}_new.tab".format(comp_name=comp_name), fout_exons="data/comparison_exons_data/{comp_name}.tab".format(comp_name=comp_name)))
    print("[features] Added PSI for comparison_exons_data/{comp_name}.tab".format(comp_name=comp_name))
    return True    

def add_dai():
    """
    Description
    -----------
    Adds delta_DAI to junction results tables
    donor_DAI = control_donor / (control_donor + test_donor)
    acceptor_DAI = control_acceptor / (control_acceptor + test_acceptor)
    delta_DAI = (donor_DAI - acceptor_DAI) * 100
    """

    count = 0
    count_all = len(annotation.comparisons)
    # read junctions, donors, anchors mapping
    mapping_junction_donor = {}
    mapping_junction_acceptor = {}
    fin_database = open("reference/junctions.tab", "rt")
    header = fin_database.readline().replace("\r", "").replace("\n", "").split("\t")
    r = fin_database.readline()
    while r:
        r = r.replace("\r", "").replace("\n", "").split("\t")
        data = dict(zip(header, r))
        mapping_junction_donor[data["junction_id"]] = data["donor_anchor_id"]
        mapping_junction_acceptor[data["junction_id"]] = data["acceptor_anchor_id"]
        r = fin_database.readline()
    fin_database.close()

    for (comp_name, comp1, comp2, _, _) in annotation.comparisons:
        fin_junctions = open("data/comparison_junctions_data/{comp_name}.tab".format(comp_name=comp_name), "rt")
        fin_donor_anchors = open("data/comparison_donor_anchors_data/{comp_name}.tab".format(comp_name=comp_name), "rt")
        fin_acceptor_anchors = open("data/comparison_acceptor_anchors_data/{comp_name}.tab".format(comp_name=comp_name), "rt")

        # read in donors
        data_donor_anchors = {}
        print("[features] reading donor anchors from: data/comparison_donor_anchors_data/{comp_name}.tab".format(comp_name=comp_name))
        header = fin_donor_anchors.readline().replace("\r", "").replace("\n", "").split("\t")
        r = fin_donor_anchors.readline()
        while r:
            r = r.replace("\r", "").replace("\n", "").split("\t")
            data = dict(zip(header, r))
            chr = data["chr"]
            strand = data["strand"]
            feature_id = data["feature_id"]
            sum_feature_test = int(data["sum_feature_test"])
            sum_feature_control = int(data["sum_feature_control"])
            data_donor_anchors[feature_id] = (sum_feature_test, sum_feature_control)
            r = fin_donor_anchors.readline()
        fin_donor_anchors.close()

        # read in anchors
        data_acceptor_anchors = {}
        print("[features] reading acceptor anchors from: data/comparison_acceptor_anchors_data/{comp_name}.tab".format(comp_name=comp_name))
        header = fin_acceptor_anchors.readline().replace("\r", "").replace("\n", "").split("\t")
        r = fin_acceptor_anchors.readline()
        while r:
            r = r.replace("\r", "").replace("\n", "").split("\t")
            data = dict(zip(header, r))
            chr = data["chr"]
            strand = data["strand"]
            feature_id = data["feature_id"]
            sum_feature_test = int(data["sum_feature_test"])
            sum_feature_control = int(data["sum_feature_control"])
            data_acceptor_anchors[feature_id] = (sum_feature_test, sum_feature_control)
            r = fin_acceptor_anchors.readline()
        fin_acceptor_anchors.close()

        fout_junctions = open("data/comparison_junctions_data/{comp_name}_new.tab".format(comp_name=comp_name), "wt")
        print("[features] processing junctions from: data/comparison_junctions_data/{comp_name}.tab".format(comp_name=comp_name))
        header = fin_junctions.readline().replace("\r", "").replace("\n", "").split("\t")
        header_out = header.copy()
        if "donor_DAI" not in header_out:
            header_out.append("donor_DAI")
        if "acceptor_DAI" not in header_out:
            header_out.append("acceptor_DAI")
        if "delta_DAI" not in header_out:
            header_out.append("delta_DAI")
        if "delta_DAI_pvalue" not in header_out:
            header_out.append("delta_DAI_pvalue")
        fout_junctions.write("\t".join(header_out) + "\n")

        r = fin_junctions.readline()
        while r:
            r = r.replace("\r", "").replace("\n", "").split("\t")
            data = dict(zip(header, r))
            exon_length = int(data["length"])
            exon_start = int(data["feature_start"])
            exon_stop = int(data["feature_stop"])
            chr = data["chr"]
            strand = data["strand"]
            junction_id = data["feature_id"]
            donor_id = mapping_junction_donor[junction_id]
            acceptor_id = mapping_junction_acceptor[junction_id]
            sum_donor_test, sum_donor_control = data_donor_anchors[donor_id]
            sum_acceptor_test, sum_acceptor_control = data_acceptor_anchors[acceptor_id]
            if (sum_donor_test+sum_donor_control)>0:
                donor_DAI = sum_donor_test / (sum_donor_test+sum_donor_control)
            else:
                donor_DAI = 0
            if (sum_acceptor_test+sum_acceptor_control)>0:
                acceptor_DAI = sum_acceptor_test / (sum_acceptor_test+sum_acceptor_control)
            else:
                acceptor_DAI = 0
            delta_DAI = (donor_DAI-acceptor_DAI)
            delta_DAI_pvalue = "nc" # not-computed
            # time costly to compute, only for >10% delta_DAI
            #if abs(delta_DAI)>0.5:
            #    _, delta_DAI_pvalue = scipy.stats.fisher_exact([[sum_donor_test, sum_donor_control], [sum_acceptor_test, sum_acceptor_control]])
            data_out = dict(zip(header_out, r))
            data_out["donor_DAI"] = "%.2f" % donor_DAI
            data_out["acceptor_DAI"] = "%.2f" % acceptor_DAI
            data_out["delta_DAI"] = "%.2f" % delta_DAI
            data_out["delta_DAI_pvalue"] = delta_DAI_pvalue
            fout_junctions.write("\t".join(str(data_out[h]) for h in header_out) + "\n")
            r = fin_junctions.readline()
        fin_junctions.close()
        fout_junctions.close()
        count += 1
        os.system("mv data/comparison_junctions_data/{comp_name}_new.tab data/comparison_junctions_data/{comp_name}.tab".format(comp_name=comp_name))
        print("[features] Added DAI for comparison_junctions_data/{comp_name}.tab ({a}/{b}, {c}% done)".format(comp_name=comp_name, a=count, b=count_all, c="%.2f" % (count/count_all*100)))
    return True    

def save_feature_data(feature_type, filter=None):
    count = 0
    count_all = len(annotation.samples)
    for sample_id in annotation.samples:
        if filter==None:
            fname = "data/sample_{feature}_data/sample_{sample_id}.tab".format(feature=feature_type, sample_id=sample_id)
        else:
            fname = "data/sample_{feature}_data/sample_{sample_id}_filtered.tab".format(feature=feature_type, sample_id=sample_id)
        fout = open(fname, "wt")
        for gene_id, gene_data in annotation.genes.items():
            if gene_data.get(feature_type, None)==None:
                continue
            for (feature_start, feature_stop), feature_data in gene_data[feature_type].items():
                row = ["{gene_id}_{gene_name}:{feature_start}-{feature_stop}".format(gene_id=gene_id, gene_name=gene_data["gene_name"], feature_start=feature_start, feature_stop=feature_stop), feature_data.get(sample_id, 0)]
                filter_key = "{chr}:{feature_start}-{feature_stop}".format(chr=gene_data["chr"], feature_start=feature_start, feature_stop=feature_stop)
                if filter==None or filter_key not in filter:
                    fout.write("\t".join(str(el) for el in row)+"\n")
        fout.close()
        count += 1
        print("[features] Saved {feature_type} data, sample={sample_id}, file={fname} ({a}/{b}, {c}% done)".format(a=count, b=count_all, c="%.2f" % (count/count_all*100), feature_type=feature_type, sample_id=sample_id, fname=fname))
