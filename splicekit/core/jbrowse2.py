import os
import splicekit.config as config
import splicekit.core as core
import http.server
import socket
import socketserver
import RangeHTTPServer

module_desc = "splicekit | jbrowse2 |"

# server part

def start():
    setup()
    server()

def server():
    hostname=socket.gethostname()
    ip_addr=socket.gethostbyname(hostname)
    Handler = RangeHTTPServer.RangeRequestHandler
    socketserver.TCPServer.allow_reuse_address = True
    with socketserver.TCPServer(("", config.jbrowse2_port), Handler) as httpd:
        print(f"{module_desc} http://{ip_addr}:{config.jbrowse2_port}/jbrowse2/?config=splicekit_data/config.json")
        httpd.serve_forever()

def setup():
    if not os.path.exists("jbrowse2"):
        os.makedirs("jbrowse2")
    if not os.path.exists("jbrowse2/version.txt"):
        os.system("wget https://github.com/GMOD/jbrowse-components/releases/download/v2.4.2/jbrowse-web-v2.4.2.zip -O jbrowse2/jbrowse-web-v2.4.2.zip")
        os.system("unzip -qq -d jbrowse2 -o jbrowse2/jbrowse-web-v2.4.2.zip")
        os.system("rm jbrowse2/jbrowse-web-v2.4.2.zip")

# JBrowse2 part

# for every bam file create bigwig and cram file plus a bed file that shows the junctions
bam_dir = config.bam_path
container = config.container
try:
    bam_files = [fi for fi in os.listdir(bam_dir) if fi.endswith('.bam')]
except:
    print(f"{module_desc} no bam files found in provided bam folder")
genome_fa = config.fasta_path
gff_fname =  config.gff3_path
dirs_to_check = ['logs/logs_jbrowse/', 'jobs/jobs_jbrowse/','results/results_jbrowse/'] 

def check_genome():
    # jbrowse requires indexed genome assembly
    if not os.path.exists(genome_fa + ".fai"):
        print(f"{module_desc} no fasta index file found, creating it")
        print(f"{module_desc} start indexing genome: {genome_fa}")
        os.system(f"{container} samtools faidx {genome_fa}")
        print(f"{module_desc} done indexing genome: {genome_fa}")

def write_sample_jobs(force_samples):
    os.system("rm -r jobs/jobs_jbrowse/* >/dev/null 2>&1") # clean up previous jobs

    # create bigwig and then cram files
    if config.platform == 'SLURM':
        job_bw="""
#!/bin/bash
#SBATCH --job-name={sample}_jbrowse               # Job name
#SBATCH --ntasks=4                               # Number of tasks
#SBATCH --nodes=1                                # All tasks on one node
#SBATCH --partition=short                        # Select queue
#SBATCH --output=logs/logs_jbrowse/{sample}.out  # Output file
#SBATCH --error=logs/logs_jbrowse/{sample}.err   # Error file

{container} bamCoverage --ignoreDuplicates --binSize {bamCoverage_binSize}  -b {bam_fname} -o {bigwig_fname} -of bigwig
{container} samtools view -C -T {genome_fa} {bam_fname} -O CRAM -o {cram_fname}
{container} samtools index {cram_fname}
        """
    
    else:

        job_bw="""
#!/bin/bash
#BSUB -J {sample}_jbrowse               # job name
#BSUB -n 4                              # number of tasks
#BSUB -R "span[hosts=1]"                # allocate all tasks to 1 host
#BSUB -q short                          # select queue
#BSUB -o logs/logs_jbrowse/{sample}.out # output file
#BSUB -e logs/logs_jbrowse/{sample}.err # error file

{container} bamCoverage --ignoreDuplicates --binSize {bamCoverage_binSize}  -b {bam_fname} -o {bigwig_fname} -of bigwig
{container} samtools view -C -T {genome_fa} {bam_fname} -O CRAM -o {cram_fname}
{container} samtools index {cram_fname}
        """

    job_sh_bw="""
{container} bamCoverage --ignoreDuplicates --binSize {bamCoverage_binSize}  -b {bam_fname} -o {bigwig_fname} -of bigwig 2> logs/logs_jbrowse/{sample}.err
{container} samtools view -C -T {genome_fa} {bam_fname} -O CRAM -o {cram_fname} 2> logs/logs_jbrowse/{sample}.err
{container} samtools index {cram_fname} 2>> logs/logs_jbrowse/{sample}.err
"""
    
    # start job creation
    for dir in dirs_to_check:
        if not os.path.exists(dir):
            os.mkdir(dir)
    
    # create jobs
    bamCoverage_binSize= 3
    fsh = open(f"jobs/jobs_jbrowse/process.sh", "wt")
    for sample in bam_files:
        bam_fname = bam_dir +'/'+ sample
        cram_fname = dirs_to_check[2] + '/'+sample.replace('.bam', '.cram')
        bigwig_fname = dirs_to_check[2] +'/'+ sample.replace('.bam', '.bw')
        if (not (os.path.exists(cram_fname) and os.path.exists(bigwig_fname))) or (force_samples == True):
            os.system(f"rm {bigwig_fname} >/dev/null 2>&1")
            os.system(f"rm {cram_fname} >/dev/null 2>&1")
            print(f"{module_desc} create sample processing job for {sample}")
            fout = open("jobs/jobs_jbrowse/{sample}.job".format(sample=sample), "wt")
            job_bw_out = job_bw.format(container=container, bamCoverage_binSize=bamCoverage_binSize, sample=sample,
                                        bam_fname=bam_fname,cram_fname=cram_fname, bigwig_fname=bigwig_fname, genome_fa=genome_fa)
            fout.write(job_bw_out)
            fout.close()        
            job_sh_bw_out = job_sh_bw.format(container=container, sample=sample, bamCoverage_binSize=bamCoverage_binSize,
                                            bam_fname=bam_fname,cram_fname=cram_fname, bigwig_fname=bigwig_fname, genome_fa=genome_fa)
            fsh.write(job_sh_bw_out)
        else:
            print(f"{module_desc} sample {sample} already processed (.bw and .cram in results/results_jbrowse) --> use 'splicekit jbrowse2 process -force' to overwrite")
    fsh.close()

def create_jbrowse_config(force_annotation):
    # add reference genome, samples files and genome annotation to config.json
    if not os.path.exists("jbrowse2/splicekit_data"):
        os.system("mkdir jbrowse2/splicekit_data")

    # initialize config.json with adding fasta file
    if not os.path.exists("jbrowse2/splicekit_data/config.json") or force_annotation==True:
        os.system(f"rm -r jbrowse2/splicekit_data/* >/dev/null 2>&1") # clean up previous jobs
        for sample in bam_files:
            os.system(f"rm logs/logs_jbrowse/{sample}.out >/dev/null 2>&1")
            os.system(f"rm logs/logs_jbrowse/{sample}.err >/dev/null 2>&1")
        print(f'{module_desc} creating config.json in jbrowse2/splicekit_data')
        if os.path.isabs(genome_fa): # check if genome_fa path is absolute or relative
            os.system(f"cd jbrowse2/splicekit_data; {container} jbrowse add-assembly {genome_fa} --name GenomeSequence --load symlink") # initialize config.json with genome sequence
        else:
            os.system(f"cd jbrowse2/splicekit_data; {container} jbrowse add-assembly ../../{genome_fa} --name GenomeSequence --load symlink") # initialize config.json with genome sequence
        # add annotation
        if gff_fname.endswith('.gz'):
            new_gff_fname = gff_fname.replace(".gz", "")
            os.system(f"{container} gunzip -c {gff_fname} >| {new_gff_fname}")
        else:
            new_gff_fname = gff_fname
        if os.path.isabs(gff_fname):
            new_gff_fname = new_gff_fname
        else:
            new_gff_fname = "../../" + new_gff_fname
        print(f"{module_desc} start annotation parsing")
        os.system(f'cd jbrowse2/splicekit_data; grep -v "^#" {new_gff_fname} | sed "s/gene/AAA/; s/mRNA/BBB/" | sort -k1,1 -k4,4n | sed "s/AAA/gene/; s/BBB/mRNA/" >| {new_gff_fname}.sorted.gff'); print('[jbrowse2] annotation sorted') 
        os.system(f'cd jbrowse2/splicekit_data; {container} bgzip {new_gff_fname}.sorted.gff -f' ); print('[jbrowse2] annotation compressed')  
        os.system(f'cd jbrowse2/splicekit_data; {container} tabix {new_gff_fname}.sorted.gff.gz'); print('[jbrowse2] annotation indexed')  
        os.system(f'cd jbrowse2/splicekit_data; {container} jbrowse add-track {new_gff_fname}.sorted.gff.gz --trackId annotation_track --name "Gene Annotations" --category "Annotations" --load symlink' ) 
        # indexes all assemblies that it can find in the current directory's config.json
        print(f"{module_desc} adding indices for text searching")
        os.system(f"cd jbrowse2/splicekit_data; {container} jbrowse text-index --out .")
        print(f"{module_desc} annotation fully parsed")

        # now we add each sample's files in different track categories  --------------------------------------------------------------------------------------------------------------------------
        for sample in bam_files:
            sample_id = sample.replace(".bam", "")
            cram_fname = dirs_to_check[2] + sample.replace(".bam", ".cram")
            bigwig_fname = dirs_to_check[2] + sample.replace(".bam", ".bw")
            os.system(f'cd jbrowse2/splicekit_data; {container} jbrowse add-track ../../{bigwig_fname} --name {sample_id}_bw --trackId {sample_id}_bw  --category "Coverage" --load symlink')
            os.system(f'cd jbrowse2/splicekit_data; {container} jbrowse add-track ../../{cram_fname} --name {sample_id}_cram --trackId {sample_id}_cram --category "Reads" --load symlink')
    else:
        print(f"{module_desc} config.json already existing in jbrowse2/splicekit_data --> use 'splicekit jbrowse2 process -force' to overwrite")

def process(force_samples=False, force_annotation=False):
    check_genome()    
    os.system(f"rm -r jobs/jobs_jbrowse/* >/dev/null 2>&1") # clean up previous jobs
    write_sample_jobs(force_samples) # if force_samples is True, then jobs are created
    # run sample jobs if available
    if not os.path.exists("jbrowse2/splicekit_data/config.json") or force_samples:
        if config.platform=="cluster":
            os.system('export BSUB_QUIET=Y; jobs=( $(ls jobs/jobs_jbrowse/*.job) ); g=10; for((i=0; i < ${#jobs[@]}; i+=g)); do part=( "${jobs[@]:i:g}" ); for job_fname in ${part[*]}; do echo "[jbrowse] submitted $job_fname"; bsub -K < ${job_fname} & done; wait; echo "[jbrowse] processing next 10"; done; echo "[jbrowse] processing complete"')
        if config.platform=="desktop":
            core.mprocess("jobs/jobs_jbrowse/process.sh")
    create_jbrowse_config(force_annotation)
