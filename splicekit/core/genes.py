# splicekit genes module
# generate gene count data
# use featureCounts and provided gtf file + bam files

import os
import sys
import splicekit.core as core
import splicekit.config as config
import gzip
import re

def split_ignore_quoted(input_str):
    # Regular expression to split by semicolons not inside quotes
    pattern = r'(?:[^";]|"(?:\\.|[^"])*")+'
    parts = re.findall(pattern, input_str)
    return [part.strip() for part in parts]

def write_genes_gtf():

    def make_row(r):
        chr = r[0]
        start, stop = int(r[3]), int(r[4]) # ! GTF file to GTF file, no change of coordinates here
        strand = r[6]
        attributes = core.split_ignore_quoted(r[-1]) # some attributes are quotes (gene_name) and can contain ; inside quotes
        new_attributes = []
        for att in attributes:
            att = att.lstrip(" ")
            att = att.split(" ")
            name, val = att[0], " ".join(att[1:])
            name = name.lstrip(" ")
            name = name.rstrip(" ")
            if name in ["gene_id", "gene_name", "transcript_id"]:
                val = val.lstrip(" ").rstrip(" ") # [1:-1] # do not remove quotes
                new_attributes.append(f"{name}={val}")
        exon_id = f"{chr}{strand}_{start}_{stop}"
        new_attributes.append("exon_id="+exon_id)
        attributes_str = '; '.join(new_attributes)
        row = '\t'.join([chr, 'splicekit', "exon", str(start), str(stop), '.', strand, '0', attributes_str])+'\n'
        return row

    """
    create genes.gtf to parse to featureCount jobs
    we need 9 columns https://www.ensembl.org/info/website/upload/gff.html/
    <seqname> <source> <feature> <start> <end> <score> <strand> <frame> [attributes]
    """
    # iterate over original gtf file
    gtf_file = gzip.open(config.gtf_path, 'rt')
    r = gtf_file.readline()
    exons_file = gzip.open("reference/genes.gtf.gz", "wt")
    while r:
        if r.startswith("#"):
            r = gtf_file.readline()
            continue
        r = r.replace("\n", "").replace("\r", "").split("\t")
        if r[2]!="exon":
            r = gtf_file.readline()
            continue
        row = make_row(r)
        exons_file.write(row)
        r = gtf_file.readline()
    gtf_file.close()
    exons_file.close()

def write_jobs_featureCounts(library_type='single-end', library_strand='NONE'):
    """
    This function write a featureCounts job per comparison for every bamfile it can find
    It takes in two parameters:
        library_type: str
        library_strand: str
    Both are needed to call featureCounts correctly (see string formating with library_type_insert and library_strand_insert variables)
    """

    #translate to be used in featureCoutns command
    library_type_insert = {"single-end":"", "paired-end":"-p "}[library_type]
    library_strand_insert = {"FIRST_READ_TRANSCRIPTION_STRAND":1, "SINGLE_STRAND":1, "SINGLE_REVERSE":1, "SECOND_READ_TRANSCRIPTION_STRAND":2, "NONE":0}[library_strand]
    
    gtf_fname = f"reference/genes.gtf.gz"
    bam_dir = f"{config.bam_path}" # files inside end with <sample_id>.bam
    out_dir = f"data/sample_genes_data"
    jobs_dir = f"jobs/count_genes"
    logs_dir = f"logs/count_genes"

    if config.platform == 'SLURM':
        job_genes="""
#!/bin/bash
#SBATCH --job-name=count_genes_{sample_id}            # Job name
#SBATCH --ntasks=12                                   # Number of tasks
#SBATCH --nodes=1                                     # All tasks on one node
#SBATCH --partition=short                             # Select queue
#SBATCH --output={logs_dir}/genes_{sample_id}.out     # Output file
#SBATCH --error={logs_dir}/genes_{sample_id}.err      # Error file    
            
{container} featureCounts {library_type_insert}-s {library_strand_insert} -M -O -T 12 -F GTF -f -t exon -g gene_id -a {gtf_fname} -o {out_fname} {sam_fname} 
# featureCount outputs command as first line of file, get rid of this first line and replace header for further parsing
# next, we are only interested in the 1st and 7th column (gene_id and count)
cp {out_fname} {out_fname}_temp
# make header line of file and overwrite out_fname as new file
echo "{header_line}" >| {out_fname}
# sum exon counts to gene counts with awk
tail -n +3 {out_fname}_temp| cut -f1,7 | awk '{{sum[$1]+=$2}} END {{OFS="\\t"; for (i in sum) print i, sum[i]}}' | sort -n >> {out_fname}
rm {out_fname}_temp
# move summary from featureCount to logs
mv {out_fname}.summary {logs_dir}/
gzip -f {out_fname}
        """

    else:
        job_genes="""
#!/bin/bash
#BSUB -J count_genes_{sample_id}            # Job name
#BSUB -n 12                                 # number of tasks
#BSUB -R "span[hosts=1]"                    # Allocate all tasks in 1 host
#BSUB -q short                              # Select queue
#BSUB -o {logs_dir}/genes_{sample_id}.out   # Output file
#BSUB -e {logs_dir}/genes_{sample_id}.err   # Error file    
    
{container} featureCounts {library_type_insert}-s {library_strand_insert} -M -O -T 12 -F GTF -f -t exon -g gene_id -a {gtf_fname} -o {out_fname} {sam_fname} 
# featureCount outputs command as first line of file, get rid of this first line and replace header for further parsing
# next, we are only interested in the 1st and 7th column (gene_id and count)
cp {out_fname} {out_fname}_temp
# make header line of file and overwrite out_fname as new file
echo "{header_line}" >| {out_fname}
# sum exon counts to gene counts with awk
tail -n +3 {out_fname}_temp| cut -f1,7 | awk '{{sum[$1]+=$2}} END {{OFS="\\t"; for (i in sum) print i, sum[i]}}' | sort -n >> {out_fname}
rm {out_fname}_temp
# move summary from featureCount to logs
mv {out_fname}.summary {logs_dir}/
gzip -f {out_fname}
        """

    job_sh_genes="""
    {container} featureCounts {library_type_insert}-s {library_strand_insert} -M -O -T 12 -F GTF -f -t exon -g gene_id -a {gtf_fname} -o {out_fname} {sam_fname} 
    cp {out_fname} {out_fname}_temp
    echo "{header_line}" >| {out_fname}
    # sum exon counts to gene counts with awk
    tail -n +3 {out_fname}_temp| cut -f1,7 | awk '{{sum[$1]+=$2}} END {{OFS="\\t"; for (i in sum) print i, sum[i]}}' | sort -n >> {out_fname}
    rm {out_fname}_temp
    mv {out_fname}.summary {logs_dir}/
    gzip -f {out_fname}
    """

    bam_files = [fi for fi in os.listdir(bam_dir) if fi.endswith('.bam')]
    sample_ids = [fi.replace('.bam', '') for fi in bam_files]
    header_line = '\t'.join(['gene_id', 'count'])
    fsh = open(f"{jobs_dir}/process.sh", "wt")
    for sample in sample_ids:
        out_fname = f"{out_dir}/sample_{sample}.tab"
        job_fname = f'{jobs_dir}/genes_{sample}.job'
        sam_fname = f"{bam_dir}/{sample}.bam"
        # cluster job
        job_out = job_genes.format(container=config.container, library_type_insert=library_type_insert, library_strand_insert=library_strand_insert, gtf_fname=gtf_fname, sample_id=sample, sam_fname=sam_fname, out_fname=out_fname, logs_dir=logs_dir, header_line=header_line)
        job_file = open(job_fname, "w")
        job_file.write(job_out)
        job_file.close()
        # shell job
        job_out = job_sh_genes.format(container=config.container, library_type_insert=library_type_insert, library_strand_insert=library_strand_insert, gtf_fname=gtf_fname, sam_fname=sam_fname, out_fname=out_fname, logs_dir=logs_dir, header_line=header_line)
        fsh.write(job_out)
    fsh.close()
