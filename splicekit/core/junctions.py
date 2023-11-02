# script to produce junction database + counts from bam files
# usage: python junctions.py filename.bam result
# outputs result.tab, result.bed

# https://www.biostars.org/p/51504/
# pysam coordinates are 0-based
# IGV help file is misleading the file is not 1 based, what they show is one based. Similarly when you open a BED file in IGV (BED is also zero based) you will see that it is drawn as a 1 based file. But the file is of course still zero based.

# GFF and GTF are 1-indexed
# http://www.ensembl.org/info/website/upload/gff.html

import os
import sys
import seaborn as sns
import matplotlib.pyplot as plt
import pysam
import argparse
import splicekit
import gzip

# per sample per raw junction count
min_count = 0

# across ALL samples in the project, include junction in master table?
try:
    min_master_junction = splicekit.config.min_master_junction 
except:
    min_master_junction = 30

cigar_dist = {}
junctions = {}
junctions_stats = []
read_junctions_dist = {}
exons_gtf = {}
genes_gtf = {}
gene_name_gid = {}
gene_id_strand = {}

def read_exons():
    print(f"[junctions] reading {splicekit.config.gtf_path} to get gene and exon coordinates")
    f = gzip.open(splicekit.config.gtf_path, "rt")
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
        pos_start, pos_stop = int(pos_start), int(pos_stop)
        pos_start -= 1
        pos_stop -= 1
        temp_atts = temp_atts.split("; ")
        for el in temp_atts:
            el = el.split(" ")
            if len(el)!=2:
                continue            
            atts[el[0]] = el[1].replace("\"", "")
        gene_id = atts["gene_id"]
        gene_name = atts.get("gene_name", "")
        gene_name_gid[gene_name] = gene_id
        gene_id_strand[gene_id] = strand
        genes_chrstrand = genes_gtf.get(f"{chr}{strand}", {})
        gene_start, gene_stop, _, _ = genes_chrstrand.get(gene_id, (float("inf"), 0, gene_id, gene_name))
        gene_start = min(gene_start, pos_start)
        gene_stop = max(gene_stop, pos_stop)
        genes_chrstrand[gene_id] = (gene_start, gene_stop, gene_id, gene_name)
        genes_gtf[f"{chr}{strand}"] = genes_chrstrand
        if r[2]=="exon":
            chr = r[0]
            exons_gtf[f"{chr}_{strand}_{pos_start}"] = (gene_id, gene_name)
            exons_gtf[f"{chr}_{strand}_{pos_stop}"] = (gene_id, gene_name)
        r = f.readline()
    f.close()

def find_genes(chr, strand, start, stop):
    detected_genes = []
    genes = genes_gtf.get(f"{chr}{strand}", {})
    for gene_id, (gene_start, gene_stop, gene_id, gene_name) in genes.items():
        if gene_start<=start<=stop<=gene_stop:
            annotation = ["N", "N"]
            if exons_gtf.get(f"{chr}_{strand}_{start}", None)!=None:
                if strand=="+":
                    annotation[0] = "A"
                else:
                    annotation[1] = "A"
            if exons_gtf.get(f"{chr}_{strand}_{stop}", None)!=None:
                if strand=="+":
                    annotation[1] = "A"
                else:
                    annotation[0] = "A"
            # this sorting classifies genes based on annotation (NN = 0..1, AN and NA = 1..2, AA = 2..3
            # shorter genes get higher values and are at the top of the list
            annotation = "".join(annotation)
            factor = {"NN":0, "NA":1, "AN":1, "AA":2}[annotation]
            gene_len = gene_stop-gene_start+1
            detected_genes.append((factor + 1.0/gene_len, gene_id, gene_name, strand, annotation))
    detected_genes.sort(reverse=True)
    return detected_genes

def make_jobs():
    job_junctions="""
#!/bin/bash
#BSUB -J {job_name}                               # Job name
#BSUB -n 4                                        # number of tasks
#BSUB -R "span[hosts=1]"                          # Allocate all tasks in 1 host
#BSUB -M 16GB                                     # Memory
#BSUB -q short                                    # Select queue
#BSUB -o logs/logs_junctions/{sample_id}.out # Output file
#BSUB -e logs/logs_junctions/{sample_id}.err # Error file

python {core_path}/junctions.py {bam_fname} data/sample_junctions_data/sample_{sample_id}
    """

    job_sh_junctions="""
python {core_path}/junctions.py {bam_fname} data/sample_junctions_data/sample_{sample_id}
    """

    fsh = open("jobs/jobs_junctions/process.sh", "wt")
    for sample_id in splicekit.core.annotation.samples:
        core_path=os.path.dirname(splicekit.core.__file__)
        bam_fname = f"{splicekit.config.bam_path}/{sample_id}.bam"
        f = open("jobs/jobs_junctions/sample_{sample_id}.job".format(sample_id=sample_id), "wt")
        f.write(job_junctions.format(sample_id=sample_id, core_path=core_path, bam_fname=bam_fname, job_name="detect_junctions_{sample_id}".format(sample_id=sample_id)))
        f.close()
        fsh.write(job_sh_junctions.format(sample_id=sample_id, core_path=core_path, bam_fname=bam_fname))
    fsh.close()
    
def detect_junctions(positions, cigar, read, pair=None):
    chr = read.reference_name
    strand = "-" if read.is_reverse else "+"
    if pair==1 and splicekit.config.library_strand=="SECOND_READ_TRANSCRIPTION_STRAND": # paired-end, second read is in transcription direction, reverse first read
        strand = {"+":"-", "-":"+"}[strand]
    if pair==2 and splicekit.config.library_strand=="FIRST_READ_TRANSCRIPTION_STRAND": # paired-end, first read is in transciption direction, reverse second read
        strand = {"+":"-", "-":"+"}[strand]
    if pair==None and splicekit.config.library_strand=="SINGLE_REVERSE": # single-end reverse direction, reverse read strand
        strand = {"+":"-", "-":"+"}[strand]
    cpos = 0
    num_junctions = 0
    for index, (t, v) in enumerate(cigar):
        if t==3:
            num_junctions += 1
            cigar_left = cigar[index-1]
            cigar_right = cigar[index+1]
            assert(cigar_left[0]==0)
            assert(cigar_right[0]==0)
            coverage_left = cigar_left[1]
            coverage_right = cigar_right[1]
            junction_size = v
            junctions_stats.append((junction_size, coverage_left, coverage_right))
            junction_start = positions[cpos-1] # last nucleotide of exon, 5'ss
            junction_stop = junction_start + v + 1  # first nucleotide of exon, 3'ss
            junction_key = (chr, strand, junction_start, junction_stop)
            junctions[junction_key] = junctions.get(junction_key, 0) + 1
        if t in [0, 1, 4, 5, 6, 7, 8, 9]:
            cpos += v
    return num_junctions

def plot_hist(data, label):
    print("plotting", label)
    plt.figure()
    sns.set(font_scale=0.5)
    sns.set_style("dark")
    sns.set_style("ticks")
    fig = sns.histplot(data=data, x=label, element="bars", discrete=True, kde=True)
    fig.spines['left'].set_linewidth(0.5)
    fig.spines['left'].set_color('#333333')
    fig.spines['bottom'].set_linewidth(0.5)
    fig.spines['bottom'].set_color('#333333')
    fig.tick_params(axis='x', colors='#333333', width=0.5)
    fig.tick_params(axis='y', colors='#333333', width=0.5)
    sns.despine()
    plt.savefig(f"{label}.png", dpi=300)
    plt.close()

def parse_sam(sam_fname, out_fname, output_bed=True):
    print(f"[junctions] parsing sam {sam_fname} to {out_fname}")
    cigar_types = {0:"M", 1:"I", 2:"D", 3:"N", 4:"S", 5:"H", 6:"P", 7:"=", 8:"X"}
    sam_fname = sys.argv[1]
    out_fname = sys.argv[2]
    count = 0
    samfile = pysam.AlignmentFile(sam_fname, "rb")
    for read in samfile.fetch():
        count += 1
        if count%1e5==0:
            print(f"[junctions] processed %.1f M alignments from {sam_fname}" % (count/1e6))
        cigar = read.cigar
        pair = None
        if read.is_paired:
            if read.is_read1:
                pair = 1
            else:
                pair = 2
        cigar_types_present = [t for (t, v) in cigar]
        cigar_length = [v for (t, v) in cigar if t==0]
        if 3 not in cigar_types_present:
            continue
        cigar_dist[len(cigar)] = cigar_dist.get(len(cigar), 0) + 1
        read_junctions = detect_junctions(read.get_reference_positions(full_length=True), cigar, read, pair=pair)
        read_junctions_dist[read_junctions] = read_junctions_dist.get(read_junctions, 0) + 1

    junction_results = []
    for (chr, strand, start, stop), count in junctions.items():
        row = [chr, strand, start, stop, count]
        junction_results.append(row)

    junction_results = sorted(junction_results, key = lambda x: (x[0], -x[-1]))    
    f = open(f"{out_fname}_raw.tab", "wt")
    header = ["chr", "strand", "start", "stop", "count"]
    f.write("\t".join(header) + "\n")
    for (chr, strand, start, stop, count) in junction_results:
        if count>=min_count:
            row = [chr, strand, start, stop, count]
            f.write("\t".join([str(el) for el in row]) + "\n")
    f.close()

    if output_bed:
        f = open(f"{out_fname}_raw.bed", "wt")
        for (chr, strand, start, stop, count) in junction_results:
            if count>=min_count:
                bed_row = [chr, start, stop, count]
                f.write("\t".join([str(el) for el in bed_row]) + "\n")
        f.close()

# given junction_id <chr><str>_<start>_<stop> and anchor_type, return corresponding anchor_id
def get_anchor_id(junction_id: str, anchor_type, anchor_width: int = 15):
    assert(anchor_type in ["donor", "acceptor"])
    coords = junction_id.split('_')
    start = int(coords[-2])
    stop = int(coords[-1])
    strand = coords[-3][-1]
    chr = '_'.join(coords[:-2])[:-1]
    if anchor_type=="donor":
        start_stop_tuple = (stop, stop + (anchor_width-1)) if (strand == '-') else (start - (anchor_width-1),start)
        anchor_id = f"{chr}{strand}_{start_stop_tuple[0]}_{start_stop_tuple[1]}"
    if anchor_type=="acceptor":
        start_stop_tuple = (stop, stop + (anchor_width-1)) if (strand == '+') else (start - (anchor_width-1),start)
        anchor_id = f"{chr}{strand}_{start_stop_tuple[0]}_{start_stop_tuple[1]}"
    return anchor_id

# alpha numeric sort
def jsort(item):
    try:
        chr = int(item[0].split("chr")[1])
    except:
        chr = float('inf')
    return (chr, item[1], int(item[2]))

def sort_junctions(data):
    return sorted(data, key=jsort)

def read_raw():
    junctions = {}
    for sample_id in splicekit.core.annotation.samples:
        raw_fname = f"data/sample_junctions_data/sample_{sample_id}_raw.tab"
        print("processing", raw_fname)
        f = open(raw_fname, "rt")
        header = f.readline().replace("\r", "").replace("\n", "").split("\t")
        r = f.readline()
        while r:
            r = r.replace("\r", "").replace("\n", "").split("\t")
            data = dict(zip(header, r))
            chr, strand, start, stop = data["chr"], data["strand"], int(data["start"]), int(data["stop"])
            count = int(data["count"])
            jkey = (chr, strand, start, stop)
            junctions[jkey] = junctions.get(jkey, 0) + count
            r = f.readline()
        f.close()
    return junctions

# create master junctions table
def make_master():
    read_exons()
    junctions = read_raw()
    junctions_annotated = []

    # stranded data protocol
    if splicekit.config.library_strand!="NONE": # stranded
        stats = {"unresolved":0, "resolved":0, "annotated_NN":0, "annotated_AN":0, "annotated_NA":0, "annotated_AA":0}
        for (chr, strand, start, stop), count in junctions.items():
            found_gene_id = ""
            found_gene_name = ""
            found_annotation = ""
            # find gene(s) for this junction
            detected_genes = find_genes(chr, strand, start, stop)
            if len(detected_genes)>0:
                found_gene_id = detected_genes[0][1]
                found_gene_name = detected_genes[0][2]
                found_strand = detected_genes[0][3]
                assert(found_strand==strand)
                found_annotation = detected_genes[0][4]
                junctions_annotated.append((chr, start, stop, strand, found_gene_id, found_gene_name, found_annotation, count))
                stats["resolved"] = stats["resolved"] + 1
                stats[f"annotated_{found_annotation}"] = stats[f"annotated_{found_annotation}"] + 1
            else:
                stats["unresolved"] = stats["unresolved"] + 1

    # unstranded data protocol
    if splicekit.config.library_strand=="NONE":
        stats = {"unresolved":0, "resolved":0, "annotated_NN":0, "annotated_AN":0, "annotated_NA":0, "annotated_AA":0}
        junctions_unstranded = {}
        # make a list of junctions without strand
        for (chr, strand, start, stop), count in junctions.items():
            jkey = (chr, start, stop)
            junctions_unstranded[jkey] = junctions_unstranded.get(jkey, 0) + count
        for (chr, start, stop), count in junctions_unstranded.items():
            # search for genes for this junction on both strands
            detected_genes = find_genes(chr, "+", start, stop) + find_genes(chr, "-", start, stop)
            detected_genes.sort(reverse=True) # resort, since we searched on both strands
            if len(detected_genes)>0:
                found_gene_id = detected_genes[0][1]
                found_gene_name = detected_genes[0][2]
                found_strand = detected_genes[0][3]
                found_annotation = detected_genes[0][4]
                junctions_annotated.append((chr, start, stop, found_strand, found_gene_id, found_gene_name, found_annotation, count))
                stats["resolved"] = stats["resolved"] + 1
                stats[f"annotated_{found_annotation}"] = stats[f"annotated_{found_annotation}"] + 1
            else:
                stats["unresolved"] = stats["unresolved"] + 1

    junctions_annotated.sort(key=jsort)
    f = open("reference/junctions.tab", "wt")
    header = ["junction_id", "donor_anchor_id", "acceptor_anchor_id", "gene_id", "gene_name", "chr", "strand", "annotated", "count"]
    f.write("\t".join(header) + "\n")
    for (chr, start, stop, strand, gene_id, gene_name, annotated, count) in junctions_annotated:
        if count>=min_master_junction or annotated==1: # write annotated junctions in the table, no matter if lowly expressed
            junction_id = f"{chr}{strand}_{start}_{stop}"
            donor_anchor_id = get_anchor_id(junction_id, anchor_type="donor")
            acceptor_anchor_id = get_anchor_id(junction_id, anchor_type="acceptor")
            row = [junction_id, donor_anchor_id, acceptor_anchor_id, gene_id, gene_name, chr, strand, annotated, count]
            f.write("\t".join([str(el) for el in row]) + "\n")
    f.close()

# read in and return reference/junctions.tab (junctions master table)
def read_junctions():
    junctions = []
    f = open("reference/junctions.tab", "rt")
    header = f.readline().replace("\r", "").replace("\n", "").split("\t")
    r = f.readline()
    while r:
        r = r.replace("\r", "").replace("\n", "").split("\t")
        data = dict(zip(header, r))
        junctions.append(data)
        r = f.readline()
    f.close()
    return junctions
   
def junctions_per_sample():
    if splicekit.config.library_strand!="NONE":
        stranded = True
    if splicekit.config.library_strand=="NONE":
        stranded = False
    junctions = read_junctions()
    for sample_id in splicekit.core.annotation.samples:
        sample_counts = {}
        print(f"[core.junctions] reading raw junctions: data/sample_junctions_data/sample_{sample_id}_raw.tab")
        f = open(f"data/sample_junctions_data/sample_{sample_id}_raw.tab", "rt")
        header = f.readline().replace("\r", "").replace("\n", "").split("\t")
        r = f.readline()
        while r:
            r = r.replace("\r", "").replace("\n", "").split("\t")
            data = dict(zip(header, r))
            if not stranded:
                junction_id = data["chr"] + "+" + "_" + data["start"] + "_" + data["stop"]
                sample_counts[junction_id] = sample_counts.get(junction_id, 0) + int(data["count"])
                junction_id = data["chr"] + "-" + "_" + data["start"] + "_" + data["stop"]
                sample_counts[junction_id] = sample_counts.get(junction_id, 0) + int(data["count"])
            else:
                junction_id = data["chr"] + data["strand"] + "_" + data["start"] + "_" + data["stop"]
                sample_counts[junction_id] = int(data["count"])
            r = f.readline()
        f.close()

        print(f"[core.junctions] writting junctions counts file: data/sample_junctions_data/sample_{sample_id}.tab")
        fout = open(f"data/sample_junctions_data/sample_{sample_id}.tab", "wt")
        header = ["junction_id", "count"]
        fout.write("\t".join(header) + "\n")
        for data in junctions:
            junction_id = data["junction_id"]
            count = sample_counts.get(junction_id, 0)
            fout.write("\t".join([junction_id, str(count)]) + "\n")
        fout.close()

if __name__ == "__main__":
    parse_sam(sys.argv[1], sys.argv[2])
