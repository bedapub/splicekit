# splicekit anchors module
# reference/junction.tab = all detected junctions and their corresponding anchors in the columns
# junction_id = <chr><strand>_<start>_<stop>
# annotated = AA/AN/NA/NN
# create GTF file from ancher coordinates for featureCount

import os
import sys
import splicekit.config as config
import gzip

def write_anchor_gtf():

    # field = donor_anchor, acceptor_anchor
    def make_row(data, field):
        ldata = data.copy()
        att_keys = ["gene_id", "gene_name", "chr", "strand", "annotated", "count", f"{field}_id", "junction_id"]
        temp_id = data[field+"_id"]
        coords = temp_id.split('_')
        start = int(coords[-2]) + 1 # GTF file 1-index coordinates (junctions.tab.gz is 0-indexed)
        stop = int(coords[-1]) + 1 # GTF file 1-index coordinates (junctions.tab.gz is 0-indexed)
        strand = coords[-3][-1]
        chr = '_'.join(coords[:-2])[:-1]
        field_id = f"{chr}{strand}_{start}_{stop}"
        ldata[f"{field}_id"] = field_id
        for key, val in ldata.items():
            if key in ["gene_id", "gene_name"]:
                ldata[key] = "\"" + ldata[key] + "\""
        attributes_str = '; '.join([key + "=" + ldata[key] for key in att_keys])
        row = '\t'.join([chr, 'splicekit', "anchor", str(start), str(stop), '.', strand, '0', attributes_str])+'\n'
        return row

    """
    create donor/acceptor_anchor.gtf to parse to featureCount jobs.
    we need 9 columns https://www.ensembl.org/info/website/upload/gff.html/
    <seqname> <source> <feature> <start> <end> <score> <strand> <frame> [attributes]
    """
    # iterate over junctions.tab.gz file and for every junction, create the anchor_[donor/acceptor].gtf.gz
    junction_file = gzip.open("reference/junctions.tab.gz", "rt")
    header = junction_file.readline().replace("\r", "").replace("\n", "").split("\t")
    r = junction_file.readline()
    anchor_files = {}
    anchor_files["donor"] = gzip.open("reference/donor_anchors.gtf.gz", "wt")
    anchor_files["acceptor"] = gzip.open("reference/acceptor_anchors.gtf.gz", "wt")
    while r:
        r = r.replace("\n", "").replace("\r", "").split("\t")
        data = dict(zip(header, r))
        row = make_row(data, "donor_anchor")
        anchor_files["donor"].write(row)
        row = make_row(data, "acceptor_anchor")
        anchor_files["acceptor"].write(row)
        r = junction_file.readline()
    junction_file.close()
    anchor_files["donor"].close()
    anchor_files["acceptor"].close()

def write_jobs_featureCounts(library_type='single-end', library_strand='NONE'):
    '''
    This function write a featureCounts job per comparison for every bamfile it can find
    It takes in two parameters:
        library_type: str
        library_strand: str
    Both are needed to call featureCounts correctly (see string formating with library_type_insert and library_strand_insert variables)
    '''

    #translate to be used in featureCoutns command
    library_type_insert = {"single-end":"", "paired-end":"-p "}[library_type]
    library_strand_insert = {"FIRST_READ_TRANSCRIPTION_STRAND":1, "SINGLE_STRAND":1, "SINGLE_REVERSE":1, "SECOND_READ_TRANSCRIPTION_STRAND":2, "NONE":0}[library_strand]
    
    for anchor_type in ["donor", "acceptor"]:
        gtf_fname = f"reference/{anchor_type}_anchors.gtf.gz"
        bam_dir = f"{config.bam_path}" # files inside end with <sample_id>.bam
        out_dir = f'data/sample_{anchor_type}_anchors_data'
        jobs_dir = f'jobs/count_{anchor_type}_anchors'
        logs_dir = f'logs/count_{anchor_type}_anchors'

        if config.platform == 'SLURM':

            job_anchors="""
#!/bin/bash
#SBATCH --job-name={anchor_type}_anchors_{sample_id}  # Job name
#SBATCH --ntasks=12                                   # Number of tasks
#SBATCH --nodes=1                                     # All tasks on one node
#SBATCH --partition=short                             # Select queue
#SBATCH --output={logs_dir}/{anchor_type}_anchors_{sample_id}.out # Output file
#SBATCH --error={logs_dir}/{anchor_type}_anchors_{sample_id}.err  # Error file

{container} featureCounts {library_type_insert}-s {library_strand_insert} -M -O -T 12 -F GTF -f -t anchor -g {anchor_type}_anchor_id -a {gtf_fname} -o {out_fname} {sam_fname} 
# featureCount outputs command as first line of file, get rid of this first line and replace header for further parsing
# next, we are only interested in the 1st and 7th column (anchor_id and count)
cp {out_fname} {out_fname}_temp
# make header line of file and overwrite out_fname as new file
echo "{header_line}" >| {out_fname}
tail -n +3 {out_fname}_temp| cut -f1,7 >> {out_fname} 
rm {out_fname}_temp
# move summary from featureCount to logs
mv {out_fname}.summary {logs_dir}/
gzip -f {out_fname}
            """
        else:

            job_anchors="""
#!/bin/bash
#BSUB -J {anchor_type}_anchors_{sample_id}  # Job name
#BSUB -n 12                                 # number of tasks
#BSUB -R "span[hosts=1]"                    # Allocate all tasks in 1 host
#BSUB -q short                              # Select queue
#BSUB -o {logs_dir}/{anchor_type}_anchors_{sample_id}.out # Output file
#BSUB -e {logs_dir}/{anchor_type}_anchors_{sample_id}.err # Error file    
    
{container} featureCounts {library_type_insert}-s {library_strand_insert} -M -O -T 12 -F GTF -f -t anchor -g {anchor_type}_anchor_id -a {gtf_fname} -o {out_fname} {sam_fname} 
# featureCount outputs command as first line of file, get rid of this first line and replace header for further parsing
# next, we are only interested in the 1st and 7th column (anchor_id and count)
cp {out_fname} {out_fname}_temp
# make header line of file and overwrite out_fname as new file
echo "{header_line}" >| {out_fname}
tail -n +3 {out_fname}_temp| cut -f1,7 >> {out_fname} 
rm {out_fname}_temp
# move summary from featureCount to logs
mv {out_fname}.summary {logs_dir}/
gzip -f {out_fname}
            """

        job_sh_anchors="""
{container} featureCounts {library_type_insert}-s {library_strand_insert} -M -O -T 12 -F GTF -f -t anchor -g {anchor_type}_anchor_id -a {gtf_fname} -o {out_fname} {sam_fname} 
cp {out_fname} {out_fname}_temp
echo "{header_line}" >| {out_fname}
tail -n +3 {out_fname}_temp| cut -f1,7 >> {out_fname} 
rm {out_fname}_temp
mv {out_fname}.summary {logs_dir}/
gzip -f {out_fname}
        """

        bam_files = [fi for fi in os.listdir(bam_dir) if fi.endswith('.bam')]
        sample_ids = [fi.replace('.bam', '') for fi in bam_files]
        header_line = '\t'.join(['anchor_id', 'count'])
        fsh = open(f"{jobs_dir}/process.sh", "wt")
        for sample in sample_ids:
            out_fname = f"{out_dir}/sample_{sample}.tab"
            job_fname = f'{jobs_dir}/{anchor_type}_anchors_{sample}.job'
            sam_fname = f"{bam_dir}/{sample}.bam"
            # cluster job
            job_out = job_anchors.format(container=config.container, library_type_insert=library_type_insert, library_strand_insert=library_strand_insert, anchor_type=anchor_type, gtf_fname=gtf_fname, sample_id=sample, sam_fname=sam_fname, out_fname=out_fname, logs_dir=logs_dir, header_line=header_line)
            job_file = open(job_fname, "w")
            job_file.write(job_out)
            job_file.close()
            # shell job
            job_out = job_sh_anchors.format(container=config.container, library_type_insert=library_type_insert, library_strand_insert=library_strand_insert, anchor_type=anchor_type, gtf_fname=gtf_fname, sam_fname=sam_fname, out_fname=out_fname, logs_dir=logs_dir, header_line=header_line)
            fsh.write(job_out)
        fsh.close()


