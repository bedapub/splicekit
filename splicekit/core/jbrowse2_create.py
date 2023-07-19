import os
import splicekit

'''
create config.json pointing to 
genome annotation
ref sequence
per sample
    bigwig
    cram
    wihtin splicekit/jbrowse2/splicekit_data <-- this is where config.json will be 
    add this config.json to URL of Jbrwose instance that is served
    IMPORTANT: all data referenced within config.json, should be in (or simlinked from) splicekit/jbrowse2/splicekit_data
'''
# create for every bam file create bigwig and cram file plus a bed file that shows the junctions
bam_dir =   splicekit.config.bam_path
container = splicekit.config.container
bam_files = [fi for fi in os.listdir(bam_dir) if fi.endswith('.bam')]
genome_fa = splicekit.config.fasta_path
gff_fname =  splicekit.config.gff3_path
dirs_to_check = ['logs/logs_jbrowse/', 'jobs/jobs_jbrowse/','results/results_jbrowse/'] 

def check_genome():
    '''jbrowse requires index genome assembly'''
    if not os.path.exists(genome_fa + '.fai'):
        print('[jbrowse] No fasta index file found. Creating it...')
        print(f'[jbrowse/samtools] Indexing genome {genome_fa}')
        os.system(f'{container} samtools faidx {genome_fa}')
        print(f'[jbrowse/samtools] Indexing genome done')


def write_sample_jobs(force_samples):
    
    # clean up previous jobs
    os.system(f'rm -r jobs/jobs_jbrowse/*')

    '''We need to create bigwig and then cram files'''

    job_bw="""
#!/bin/bash
#BSUB -J {sample}_jbrowse                              # Job name
#BSUB -n 4                                       # number of tasks
#BSUB -R "span[hosts=1]"                         # Allocate all tasks in 1 host
#BSUB -q short                                   # Select queue
#BSUB -o logs/logs_jbrowse/{sample}.out # Output file
#BSUB -e logs/logs_jbrowse/{sample}.err # Error file

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
            os.system(f'rm {bigwig_fname} >/dev/null 2>&1')
            os.system(f'rm {cram_fname} >/dev/null 2>&1')
            print(f'[jbrowse] create sample processing job for {sample}')
            fout = open("jobs/jobs_jbrowse/{sample}.job".format(sample=sample), "wt")
            job_bw_out = job_bw.format(container=container, bamCoverage_binSize=bamCoverage_binSize, sample=sample,
                                        bam_fname=bam_fname,cram_fname=cram_fname, bigwig_fname=bigwig_fname, genome_fa=genome_fa)
            fout.write(job_bw_out)
            fout.close()        
            job_sh_bw_out = job_sh_bw.format(container=container, sample=sample, bamCoverage_binSize=bamCoverage_binSize,
                                            bam_fname=bam_fname,cram_fname=cram_fname, bigwig_fname=bigwig_fname, genome_fa=genome_fa)
            fsh.write(job_sh_bw_out)
        else:
            print(f'[jbrowse] sample {sample} already processed (.bw and .cram in results/results_jbrowse) --> use "splicekit jbrowse2_create samples -force" to overwrite')
    fsh.close()

def create_jbrowse_config(force_annotation):
    ''' We need to add reference genome, samples files and genome annotation to jbrowse's config.json and create it first'''
    if not os.path.exists('jbrowse2/splicekit_data'):
        os.system('mkdir jbrowse2/splicekit_data') # first create folder to put everything inside

    # initialize config.json with adding fasta file --------------------------------------------------------------------------------------------------------------------------------------
    if (not os.path.exists('jbrowse2/splicekit_data/config.json')) or (force_annotation==True):
        
        # clean up previous jobs
        os.system(f'rm -r jbrowse2/splicekit_data/*')
        for sample in bam_files:
            os.system(f'rm logs/logs_jbrowse/{sample}.out >/dev/null 2>&1')
            os.system(f'rm logs/logs_jbrowse/{sample}.err >/dev/null 2>&1')

        print('[jbrowse] creating config.json in jbrowse2/splicekit_data')

        if os.path.isabs(genome_fa): # check if genome_fa path is absolute or relative
            os.system(f'cd jbrowse2/splicekit_data; {container} jbrowse add-assembly {genome_fa} --name GenomeSequence --load symlink') # initialize config.json with genome sequence
        else:
            os.system(f'cd jbrowse2/splicekit_data; {container} jbrowse add-assembly ../../{genome_fa} --name GenomeSequence --load symlink') # initialize config.json with genome sequence

        # add the annotation ------------------------------------------------------------------------------------------------------------------------------------------------------------------
        if gff_fname.endswith('.gz'):
            new_gff_fname = gff_fname.replace(".gz","")
            os.system(f'{container} gunzip -c {gff_fname} >| {new_gff_fname}')
        else:
            new_gff_fname = gff_fname

        if os.path.isabs(gff_fname):
            new_gff_fname=new_gff_fname
        else:
            new_gff_fname='../../'+new_gff_fname

        print('[jbrowse2] start annotation parsing')
        os.system(f'cd jbrowse2/splicekit_data; grep -v "^#" {new_gff_fname} | sed "s/gene/AAA/; s/mRNA/BBB/" | sort -k1,1 -k4,4n | sed "s/AAA/gene/; s/BBB/mRNA/" >| {new_gff_fname}.sorted.gff'); print('[jbrowse2] annotation sorted') 
        os.system(f'cd jbrowse2/splicekit_data; {container} bgzip {new_gff_fname}.sorted.gff -f' ); print('[jbrowse2] annotation compressed')  
        os.system(f'cd jbrowse2/splicekit_data; {container} tabix {new_gff_fname}.sorted.gff.gz'); print('[jbrowse2] annotation indexed')  
        os.system(f'cd jbrowse2/splicekit_data; {container} jbrowse add-track {new_gff_fname}.sorted.gff.gz --trackId annotation_track --name "Gene Annotations" --category "Annotations" --load symlink' ) 
        # indexes all assemblies that it can find in the current directory's config.json
        os.system(f'cd jbrowse2/splicekit_data; {container} jbrowse text-index --out .'); print('adding indices for text searching')
        print('[jbrowse2] annotation fully parsed')

        # now we add each sample's files in different track categories  --------------------------------------------------------------------------------------------------------------------------
        for sample in bam_files:
            sample_id = sample.replace('.bam', '')
            cram_fname = dirs_to_check[2] + sample.replace('.bam', '.cram')
            bigwig_fname = dirs_to_check[2] + sample.replace('.bam', '.bw')
            os.system(f'cd jbrowse2/splicekit_data; {container} jbrowse add-track ../../{bigwig_fname} --name {sample_id}_bw --trackId {sample_id}_bw  --category "Coverage" --load symlink')
            os.system(f'cd jbrowse2/splicekit_data; {container} jbrowse add-track ../../{cram_fname} --name {sample_id}_cram --trackId {sample_id}_cram --category "Reads" --load symlink')
    else:
        print('[jbrowse] config.json already existing in jbrowse2/splicekit_data --> use "splicekit jbrowse2_create annotation -force" to overwrite')

def process(force_samples=False, force_annotation=False):
    #os.system(f'rm -r jbrowse2/splicekit_data/*') # clean up before starting
    check_genome()
    # clean up previous jobs
    os.system(f'rm -r jobs/jobs_jbrowse/*')
    write_sample_jobs(force_samples) # if force_samples is True, then jobs are created
    # run  sample jobs if available
    if (not os.path.exists('jbrowse2/splicekit_data/config.json')) or (force_samples):
        if splicekit.config.platform=="cluster":
            os.system('export BSUB_QUIET=Y; jobs=( $(ls jobs/jobs_jbrowse/*.job) ); g=10; for((i=0; i < ${#jobs[@]}; i+=g)); do part=( "${jobs[@]:i:g}" ); for job_fname in ${part[*]}; do echo "[jbrowse] submitted $job_fname"; bsub -K < ${job_fname} & done; wait; echo "[jbrowse] processing next 10"; done; echo "[jbrowse] processing complete"')
        if splicekit.config.platform=="desktop":
            os.system("source jobs/jobs_jbrowse/process.sh")
    create_jbrowse_config(force_annotation)
