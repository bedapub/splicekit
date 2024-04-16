# module features
# generates count tables: juncion, exon, genes, anchor

import os
import sys
import math
import gzip
import pybio
import splicekit.config as config
import splicekit.core.annotation as annotation
import scipy
import scipy.stats

module_name = "splicekit | features |"

# load genes with exons, junctions and anchors
def load_genes():
    annotation.junctions_genes = {}
    annotation.donor_anchors_genes = {}
    annotation.acceptor_anchors_genes = {}
    annotation.exons_genes = {}
    transcript_exons = {}
    print(f"{module_name} reading annotation: {config.gtf_path}")
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
        if gene["gene_name"]=="":
            gene["gene_name"] = atts.get("gene", "")
            if gene["gene_name"]=="":
                gene["gene_name"] = gene_id
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

    # save last nucleotide of first exon of each transcript
    print(f"{module_name} saving last nucleotide of first/second exon of each transcript")
    annotation.first_exons = {}
    identify_exons(transcript_exons, annotation.first_exons, 0)
    
    print(f"{module_name} reading junctions and anchors annotation from: reference/junctions.tab.gz")
    f = gzip.open("reference/junctions.tab.gz", "rt")
    header = f.readline().replace("\r", "").replace("\n", "").split("\t")
    r = f.readline()
    while r:
        r = r.replace("\r", "").replace("\n", "").split("\t")
        data = dict(zip(header, r))
        gene_id = data["gene_id"]
        gene = annotation.genes.get(gene_id, {})
        gene["chr"], gene["strand"] = data["chr"], data["strand"]
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
    count = 0
    count_all = len(annotation.samples)
    for sample_id in annotation.samples:
        f = gzip.open(f"data/sample_junctions_data/sample_{sample_id}.tab.gz", "rt")
        r = f.readline() # header
        r = f.readline()
        while r:
            r = r.replace("\r", "").replace("\n", "").split("\t")
            junction_id = r[0]
            coords = junction_id.split('_')
            start, stop = int(coords[-2]), int(coords[-1])
            chr, strand = '_'.join(coords[:-2])[:-1], coords[-3][-1]
            raw_count = int(r[-1])
            gene_id = annotation.junctions_genes[junction_id]
            gene = annotation.genes.get(gene_id, {})
            junctions = gene.get("junctions", {})
            junctions.setdefault(junction_id, {})[sample_id] = raw_count # sample counts dict
            gene["junctions"] = junctions
            annotation.genes[gene_id] = gene
            r = f.readline()
        count += 1
        print(f"{module_name} junction read OK, sample={sample_id}, from=data/sample_junctions_data/sample_{sample_id}.tab {count}/{count_all}")
        f.close()

def read_anchors(anchor_type):
    """
    * read anchors sample raw counts and populate the annotation.genes ["anchors"] field
    * input is read from the output of juan scripts (data/sample_anchors_data) producing anchor read counts
    """

    count = 0
    count_all = len(annotation.samples)
    for sample_id in annotation.samples:
        f = gzip.open(f"data/sample_{anchor_type}_anchors_data/sample_{sample_id}.tab.gz", "rt")
        r = f.readline() # header
        r = f.readline()
        while r:
            r = r.replace("\r", "").replace("\n", "").split("\t")
            anchor_id = r[0]
            coords = anchor_id.split('_')
            start, stop = int(coords[-2])-1, int(coords[-1])-1 # reports by featureCounts are on 1-indexed gtf regions, convert back
            chr, strand = '_'.join(coords[:-2])[:-1], coords[-3][-1]
            anchor_id = f"{chr}{strand}_{start}_{stop}" # reconstruct anchor_id with 0-indexed coords
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
        print(f"{module_name} anchors read OK, sample={sample_id}, from=data/sample_{anchor_type}_anchors_data/sample_{sample_id}.tab.gz {count}/{count_all}")
        f.close()

def read_exons():
    """
    * read exons sample raw counts and populate the annotation.genes ["exons"] field
    * input is read from data/sample_exons_data
    """
    count = 0
    count_all = len(annotation.samples)
    for sample_id in annotation.samples:
        f = gzip.open(f"data/sample_exons_data/sample_{sample_id}.tab.gz", "rt")
        r = f.readline() # header
        r = f.readline()
        while r:
            r = r.replace("\r", "").replace("\n", "").split("\t")
            exon_id = r[0]
            coords = exon_id.split('_')
            start, stop = int(coords[-2])-1, int(coords[-1])-1 # reports by featureCounts are on 1-indexed gtf regions, convert back
            chr, strand = '_'.join(coords[:-2])[:-1], coords[-3][-1]
            exon_id = f"{chr}{strand}_{start}_{stop}" # reconstruct exon_id with 0-indexed coords
            raw_count = int(r[-1])
            gene_id = annotation.exons_genes[exon_id]
            gene = annotation.genes.get(gene_id, {})
            exons = gene.get("exons", {})
            exons.setdefault(exon_id, {})[sample_id] = raw_count # sample counts dict
            gene["exons"] = exons
            annotation.genes[gene_id] = gene
            r = f.readline()
        count += 1
        print(f"{module_name} exons read OK, sample={sample_id}, from=data/sample_exons_data/sample_{sample_id}.tab.gz, ({count}/{count_all})")
        f.close()

def read_genes():
    """
    * read gene sample raw counts and populate the annotation.genes ["genes"] field
    * input is read from data/sample_genes_data
    """
    count = 0
    count_all = len(annotation.samples)
    for sample_id in annotation.samples:
        f = gzip.open(f"data/sample_genes_data/sample_{sample_id}.tab.gz", "rt")
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
        print(f"{module_name} genes read OK, sample={sample_id}, from=data/sample_genes_data/sample_{sample_id}.tab.gz, ({count}/{count_all}")
        f.close()

def save_comps_feature_data(feature_type):
    """
    * take files from data/samples_{feature_type}_data/*.tab and read in the counts
    * write the data/comparison_{feature_type}_data files
    * feature_type = exons, junctions, donor_anchors, acceptor_anchors
    """

    count = 0
    count_all = len(annotation.comparisons)
    for (comp_name, comp1, comp2, _, _) in annotation.comparisons:
        fout = gzip.open(f"data/comparison_{feature_type}_data/{comp_name}.tab.gz", "wt")
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
        print(f"{module_name} saved comparison_{feature_type}_data/{comp_name}.tab.gz ({count}/{count_all})")
    return True

def make_counts_table(feature_type):
    # general matrix (all samples), with a slightly different structure (no control/test sum etc.)
    # all sample ids are in: annotation.samples
    print(f"{module_name} make counts table: data/samples_{feature_type}_counts.tab.gz")
    fout = gzip.open(f"data/samples_{feature_type}_counts.tab.gz", "wt")
    header = ["gene_id", "gene_name", "chr", "strand", "feature_start", "feature_stop", "length", "feature_id"]
    for sample_id in annotation.samples:
        header.append(sample_id)
    fout.write("\t".join(header)+"\n")
    for gene_id, gene_data in annotation.genes.items():
        if gene_data.get(feature_type, None)==None:
            continue
        for feature_id, feature_data in gene_data[feature_type].items():
            # propagate gene start stop from gene object to count object
            if feature_type=="genes" and feature_id in ["start", "stop"]:
                continue
            if feature_type=="genes":
                feature_data["start"] = gene_data["genes"]["start"]
                feature_data["stop"] = gene_data["genes"]["stop"]
        for feature_id, feature_data in gene_data[feature_type].items():
            # propagate gene start stop from gene object to count object
            if feature_type=="genes" and feature_id in ["start", "stop"]:
                continue
            if feature_type=="genes":
                feature_data["start"] = gene_data["genes"]["start"]
                feature_data["stop"] = gene_data["genes"]["stop"]
            feature_len = feature_data["stop"]-feature_data["start"]+1
            row = [gene_id, f"{gene_data['gene_name']}\t{gene_data['chr']}\t{gene_data['strand']}\t{feature_data['start']}\t{feature_data['stop']}\t{feature_len}"]
            row.append(feature_id)
            for sample_id in annotation.samples:
                row.append(feature_data.get(sample_id, 0))
            fout.write("\t".join(str(el) for el in row)+"\n")
    fout.close()

def add_psi_cluster():
    if config.platform=="cluster":
        commands = []
        commands.append("export BSUB_QUIET=Y")
        for (comp_name, comp1, comp2, _, _) in annotation.comparisons:
            commands.append("bsub -q short -o /dev/null -e /dev/null -K python -c 'import splicekit; splicekit.core.features.add_psi(\"" + comp_name + "\")'")
        commands = "export BSUB_QUIET=Y; " + " & ".join(commands) + "; wait"
        os.system(commands)
    if config.platform=="desktop":
        for (comp_name, comp1, comp2, _, _) in annotation.comparisons:
            command = "python -c 'import splicekit; splicekit.core.features.add_psi(\"" + comp_name + "\")'"
            os.system(command)

def add_psi(comp_name):
    read_length = 150
    """
    * adds delta_PSI to exon results tables
    * reference: https://currentprotocols.onlinelibrary.wiley.com/doi/full/10.1002/0471142905.hg1116s87

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

    fin_exons = gzip.open(f"data/comparison_exons_data/{comp_name}.tab.gz", "rt")
    fin_junctions = gzip.open(f"data/comparison_junctions_data/{comp_name}.tab.gz", "rt")
    data_junctions = {} # chr / strand + list of junctions
    # first read in the junctions
    print(f"{module_name} reading junctions from: data/comparison_junctions_data/{comp_name}.tab.gz")
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

    fout_exons = gzip.open(f"data/comparison_exons_data/{comp_name}.tab.gz.temp", "wt")
    print(f"{module_name} processing exons from: data/comparison_exons_data/{comp_name}.tab.gz")
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
        try:
            control_PSI = control_IR/(control_IR + control_ER) * 100
        except:
            control_PSI = 0
        delta_PSI = test_PSI - control_PSI
        data_out = dict(zip(header_out, r))
        data_out["test_PSI"] = "%.2f" % test_PSI
        data_out["control_PSI"] = "%.2f" % control_PSI
        data_out["delta_PSI"] = "%.2f" % delta_PSI
        fout_exons.write("\t".join(str(data_out[h]) for h in header_out) + "\n")
        r = fin_exons.readline()
    fin_exons.close()
    fout_exons.close()
    os.system(f"mv data/comparison_exons_data/{comp_name}.tab.gz.temp data/comparison_exons_data/{comp_name}.tab.gz")
    print(f"{module_name} added PSI for comparison_exons_data/{comp_name}.tab.gz")
    return True    

def add_dai():
    """
    * adds delta_DAI to junction results tables
    * donor_DAI = control_donor / (control_donor + test_donor)
    * acceptor_DAI = control_acceptor / (control_acceptor + test_acceptor)
    * delta_DAI = (donor_DAI - acceptor_DAI) * 100
    """

    count = 0
    count_all = len(annotation.comparisons)
    # read junctions, donors, anchors mapping
    mapping_junction_donor = {}
    mapping_junction_acceptor = {}
    fin_database = gzip.open("reference/junctions.tab.gz", "rt")
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
        fin_junctions = gzip.open(f"data/comparison_junctions_data/{comp_name}.tab.gz", "rt")
        fin_donor_anchors = gzip.open("data/comparison_donor_anchors_data/{comp_name}.tab.gz", "rt")
        fin_acceptor_anchors = gzip.open("data/comparison_acceptor_anchors_data/{comp_name}.tab.gz", "rt")

        # read in donors
        data_donor_anchors = {}
        print(f"{module_name} reading donor anchors from: data/comparison_donor_anchors_data/{comp_name}.tab.gz")
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
        print(f"{module_name} reading acceptor anchors from: data/comparison_acceptor_anchors_data/{comp_name}.tab.gz")
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

        fout_junctions = gzip.open("data/comparison_junctions_data/{comp_name}.tab.gz.temp", "wt")
        print(f"{module_name} processing junctions from: data/comparison_junctions_data/{comp_name}.tab.gz")
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
        os.system(f"mv data/comparison_junctions_data/{comp_name}.tab.gz.temp data/comparison_junctions_data/{comp_name}.tab.gz")
        print(f"{module_name} Added DAI for comparison_junctions_data/{comp_name}.tab.gz ({count}/{count_all}")
    return True