import os
import sys
import argparse
import splicekit.core as core

# load version
splicekit_path = os.path.abspath(__file__)
splicekit_folder = os.path.dirname(splicekit_path)
version = open(os.path.join(splicekit_folder, "version"), "rt").readlines()[0].replace("\n", "").replace("\r", "")

if not os.path.exists("splicekit.config"):
    print("splicekit | please run splicekit in a folder with splicekit.config present")
    sys.exit(0)

parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter, add_help=False)
parser.add_argument('command', help="command(s) to run", nargs='*')
parser.add_argument("-help", "-h", "--help", action="store_true")
parser.add_argument("-version", "--version", help="Print version", action="store_true")
parser.add_argument("-force", "--force", help="Force recreation / overwrite of results", action="store_true", default=False)
parser.add_argument("-hostname", "--hostname", help="If you are submitting 'splicekit process' to the cluster, you can provide a custom hostname. This will be used to create JBrowse2 URLs.", default=None)
parser.add_argument("-verbose", "--verbose", help="Verbose mode", action="store_true", default=False)
args, unknown_args = parser.parse_known_args()

verbose = args.verbose

import splicekit.config as config
import splicekit.core as core
verbose and print("splicekit | loading splicekit.core.annotation")
import splicekit.core.annotation
verbose and print("splicekit | loading splicekit.core.features")
import splicekit.core.features
verbose and print("splicekit | loading splicekit.core.exons")
import splicekit.core.exons
verbose and print("splicekit | loading splicekit.core.genes")
import splicekit.core.genes
verbose and print("splicekit | loading splicekit.core.report")
import splicekit.core.report
verbose and print("splicekit | loading splicekit.core.patterns")
import splicekit.core.patterns
verbose and print("splicekit | loading splicekit.core.promisc")
import splicekit.core.promisc
verbose and print("splicekit | loading splicekit.core.anchors ")
import splicekit.core.anchors 
verbose and print("splicekit | loading splicekit.core.junctions")
import splicekit.core.junctions
verbose and print("splicekit | loading splicekit.core.jbrowse2")
import splicekit.core.jbrowse2
verbose and print("splicekit | loading splicekit.core.juan")
import splicekit.core.juan
verbose and print("splicekit | loading splicekit.judge")
import splicekit.judge
verbose and print("splicekit | loading splicekit.core.delta_dar")
import splicekit.core.delta_dar
verbose and print("splicekit | loading splicekit.clusterlogfc")
import splicekit.clusterlogfc
verbose and print("splicekit | loading splicekit.june")
import splicekit.june
verbose and print("splicekit | loading splicekit.report")
import splicekit.report

# initialization (triggered once on "import splicekit")
splicekit.core.annotation.compounds = {}
splicekit.core.annotation.samples = set()
splicekit.core.annotation.comparisons = []
splicekit.core.annotation.genes = {}
splicekit.core.annotation.make_comparisons() # load existing comparisons

def setup():
    core_path = os.path.dirname(splicekit.__file__)
    folder_fname = os.path.join(core_path, "folders.setup")
    folders = open(folder_fname).read().split("\n")
    # create folders
    for folder_name in folders:
        if folder_name=="":
            continue
        os.system("mkdir -p {folder_name}".format(folder_name=folder_name))
    print("splicekit | setup | successfully setup splicekit v{version} analysis in {folder}\n".format(folder=os.getcwd(), version=version))

def clean():
    os.system("rm -f results/edgeR/junctions/*.tab > /dev/null 2>&1")
    os.system("rm -f jobs/edgeR/junctions/*.job > /dev/null 2>&1")
    os.system("rm -f results/edgeR/exons/*.tab > /dev/null 2>&1")
    os.system("rm -f jobs/edgeR/exons/*.job > /dev/null 2>&1")

def annotation():
    folders = ["annotation", "jobs/edgeR/junctions", "jobs/edgeR/exons", "jobs/edgeR/donor_anchors", "jobs/edgeR/acceptor_anchors", "jobs/edgeR/genes"]
    for folder in folders:
        os.system(f"rm -f {folder}/* > /dev/null 2>&1")
    splicekit.core.annotation.make_comparisons()
    splicekit.core.annotation.write_comparisons()
    splicekit.core.annotation.make_design_contrast()

def junctions_master():
    splicekit.core.annotation.make_comparisons()
    splicekit.core.junctions.make_jobs()
    if splicekit.config.platform=="cluster":
        os.system('export BSUB_QUIET=Y; jobs=( $(ls jobs/count_junctions/*.job) ); g=10; for((i=0; i < ${#jobs[@]}; i+=g)); do part=( "${jobs[@]:i:g}" ); for job_fname in ${part[*]}; do echo "splicekit | features | junctions | submitted $job_fname"; bsub -M 8GB -K < ${job_fname} & done; wait; echo "splicekit | features | junctions | processing next 10"; done; echo "splicekit | features | junctions | processing complete"')
        os.system("export BSUB_QUIET=Y; bsub -q short -M 16GB -o /dev/null -e /dev/null -K python -c 'import splicekit; splicekit.core.junctions.make_master()'; wait;") # run make master as cluster job
    if splicekit.config.platform=="desktop":
        splicekit.core.mprocess("jobs/count_junctions/process.sh")
        splicekit.core.junctions.make_master()
    if splicekit.config.platform=="SLURM":
        os.system("jobs=( $(ls jobs/count_junctions/*.job) ); g=10; " "for((i=0; i < ${#jobs[@]}; i+=g)); do " "part=( \"${jobs[@]:i:g}\" ); " "for job_fname in ${part[*]}; do " "echo \"splicekit | features | junctions | submitted $job_fname\"; " "sbatch --mem=8G --open-mode=append ${job_fname} & " "done; wait; " "echo \"splicekit | features | junctions | processing next 10\"; " "done; " "echo \"splicekit | features | junctions | processing complete\"")
        os.system("sbatch --partition=short --mem=16G --output=/dev/null --error=/dev/null " "--wrap=\"python -c 'import splicekit; splicekit.core.junctions.make_master()'\"")

def junctions():
    splicekit.core.annotation.make_comparisons()
    splicekit.core.junctions.junctions_per_sample()
    splicekit.core.features.load_genes()
    splicekit.core.features.read_junctions()
    os.system("rm -f data/comparison_junctions_data/*.tab.gz > /dev/null 2>&1")
    #splicekit.core.features.save_comps_feature_data("junctions")
    splicekit.core.features.make_counts_table("junctions")
    # TODO
    #splicekit.core.features.add_dai() # add dai to junctions

def exons():
    os.system(f"rm -f reference/exons.gtf > /dev/null 2>&1")
    os.system(f"rm -f jobs/count_exons/*.job > /dev/null 2>&1")
    os.system(f"rm -f logs/count_exons/* > /dev/null 2>&1")
    splicekit.core.exons.write_exons_gtf()
    splicekit.core.exons.write_jobs_featureCounts(library_type=config.library_type, library_strand=config.library_strand)
    if splicekit.config.platform=="cluster":
        os.system('export BSUB_QUIET=Y; jobs=( $(ls jobs/count_exons/*.job) ); g=10; for((i=0; i < ${#jobs[@]}; i+=g)); do part=( "${jobs[@]:i:g}" ); for job_fname in ${part[*]}; do echo "splicekit | features | exons | submitted $job_fname"; bsub -M 8GB -K < ${job_fname} & done; wait; echo "splicekit | features | exons | processing next 10"; done; echo "[exons] processing complete"')
    if splicekit.config.platform=="desktop":
        os.system(". jobs/count_exons/process.sh")
    if splicekit.config.platform=="SLURM":
        os.system('jobs=( $(ls jobs/count_exons/*.job) ); g=10; for((i=0; i < ${#jobs[@]}; i+=g)); do part=( "${jobs[@]:i:g}" ); job_ids=(); for job_fname in "${part[@]}"; do echo "splicekit | features | exons | submitted $job_fname"; job_id=$(sbatch --mem=8G --parsable ${job_fname}); job_ids+=($job_id); done; for job_id in "${job_ids[@]}"; do scontrol show job $job_id | grep -q "JobState=COMPLETED" || scontrol wait job $job_id; done; echo "splicekit | features | exons | processing next 10"; done; echo "[exons] processing complete"')
    splicekit.core.annotation.make_comparisons()
    splicekit.core.features.load_genes()
    splicekit.core.features.read_exons()
    os.system("rm -f data/comparison_exons_data/*.tab > /dev/null 2>&1")
    #splicekit.core.features.save_comps_feature_data("exons")
    splicekit.core.features.make_counts_table("exons")
    # TODO
    #splicekit.core.features.add_psi_cluster() # cluster, each comparison one job

def genes():
    os.system(f"rm -f reference/genes.gtf > /dev/null 2>&1")
    os.system(f"rm -f jobs/count_genes/*.job > /dev/null 2>&1")
    os.system(f"rm -f logs/count_genes/* > /dev/null 2>&1")
    splicekit.core.genes.write_genes_gtf()
    splicekit.core.genes.write_jobs_featureCounts(library_type=config.library_type, library_strand=config.library_strand)
    if splicekit.config.platform=="cluster":
        os.system('export BSUB_QUIET=Y; jobs=( $(ls jobs/count_genes/*.job) ); g=10; for((i=0; i < ${#jobs[@]}; i+=g)); do part=( "${jobs[@]:i:g}" ); for job_fname in ${part[*]}; do echo "[genes] submitted $job_fname"; bsub -M 8GB -K < ${job_fname} & done; wait; echo "splicekit | features | genes | processing next 10"; done; echo "splicekit | features | genes | processing complete"')
    if splicekit.config.platform=="desktop":
        os.system(". jobs/count_genes/process.sh")
    if splicekit.config.platform=="SLURM":
        os.system('jobs=( $(ls jobs/count_genes/*.job) ); g=10; for((i=0; i < ${#jobs[@]}; i+=g)); do part=( "${jobs[@]:i:g}" ); job_ids=(); for job_fname in "${part[@]}"; do echo "[genes] submitted $job_fname"; job_id=$(sbatch --mem=8G --parsable ${job_fname}); job_ids+=($job_id); done; for job_id in "${job_ids[@]}"; do scontrol show job $job_id | grep -q "JobState=COMPLETED" || scontrol wait job $job_id; done; echo "splicekit | features | genes | processing next 10"; done; echo "splicekit | features | genes | processing complete"')
    splicekit.core.annotation.make_comparisons()
    splicekit.core.features.load_genes()
    splicekit.core.features.read_genes()
    os.system("rm -f data/comparison_genes_data/*.tab.gz > /dev/null 2>&1")
    splicekit.core.features.make_counts_table("genes")

def features():
    junctions_master()
    anchors()
    junctions()
    exons()
    genes()

def anchors():
    for anchor_type in ["donor", "acceptor"]:
        os.system(f"rm -f reference/{anchor_type}_anchors.gtf > /dev/null 2>&1")
        os.system(f"rm -f jobs/count_{anchor_type}_anchors/*.job > /dev/null 2>&1")
        os.system(f"rm -f logs/count_{anchor_type}_anchors/* > /dev/null 2>&1")
    splicekit.core.anchors.write_anchor_gtf()
    splicekit.core.anchors.write_jobs_featureCounts(library_type=config.library_type, library_strand=config.library_strand) # takes ~ 10 minutes
    if splicekit.config.platform=="cluster":
        os.system('export BSUB_QUIET=Y; jobs=( $(ls jobs/count_donor_anchors/*.job) ); g=10; for((i=0; i < ${#jobs[@]}; i+=g)); do part=( "${jobs[@]:i:g}" ); for job_fname in ${part[*]}; do echo "[features] submitted $job_fname"; bsub -M 8GB -K < ${job_fname} & done; wait; echo "splicekit | features | donor_anchors | processing next 10"; done; echo "splicekit | features | donor_anchors | processing complete"')
        os.system('export BSUB_QUIET=Y; jobs=( $(ls jobs/count_acceptor_anchors/*.job) ); g=10; for((i=0; i < ${#jobs[@]}; i+=g)); do part=( "${jobs[@]:i:g}" ); for job_fname in ${part[*]}; do echo "[features] submitted $job_fname"; bsub -M 8GB -K < ${job_fname} & done; wait; echo "splicekit | features | acceptor_anchors | processing next 10"; done; echo "splicekit | features | acceptor_anchors | processing complete"')
    if splicekit.config.platform=="desktop":
        os.system(". jobs/count_donor_anchors/process.sh")
        os.system(". jobs/count_acceptor_anchors/process.sh")
    if splicekit.config.platform=="SLURM":
        os.system('jobs=( $(ls jobs/count_donor_anchors/*.job) ); g=10; for((i=0; i < ${#jobs[@]}; i+=g)); do part=( "${jobs[@]:i:g}" ); job_ids=(); for job_fname in "${part[@]}"; do echo "[features] submitted $job_fname"; job_id=$(sbatch --mem=8G --parsable ${job_fname}); job_ids+=($job_id); done; for job_id in "${job_ids[@]}"; do scontrol show job $job_id | grep -q "JobState=COMPLETED" || scontrol wait job $job_id; done; echo "splicekit | features | donor_anchors | processing next 10"; done; echo "splicekit | features | donor_anchors | processing complete"')
        os.system('jobs=( $(ls jobs/count_acceptor_anchors/*.job) ); g=10; for((i=0; i < ${#jobs[@]}; i+=g)); do part=( "${jobs[@]:i:g}" ); job_ids=(); for job_fname in "${part[@]}"; do echo "[features] submitted $job_fname"; job_id=$(sbatch --mem=8G --parsable ${job_fname}); job_ids+=($job_id); done; for job_id in "${job_ids[@]}"; do scontrol show job $job_id | grep -q "JobState=COMPLETED" || scontrol wait job $job_id; done; echo "splicekit | features | acceptor_anchors | processing next 10"; done; echo "splicekit | features | acceptor_anchors | processing complete"')
    splicekit.core.annotation.make_comparisons()
    splicekit.core.features.load_genes()
    splicekit.core.features.read_anchors("donor")
    splicekit.core.features.read_anchors("acceptor")
    os.system("rm -f data/comparison_donor_anchors_data/*.tab.gz > /dev/null 2>&1")
    os.system("rm -f data/comparison_acceptor_anchors_data/*.tab.gz > /dev/null 2>&1")
    splicekit.core.features.make_counts_table("donor_anchors")
    splicekit.core.features.make_counts_table("acceptor_anchors")

def patterns():
    splicekit.core.patterns.process()

def edgeR(run=None):

    if run=="junctions" or run==None:
        splicekit.core.annotation.make_comparisons()
        os.system(f"rm -f results/edgeR/junctions/*.tab.gz > /dev/null 2>&1")
        if splicekit.config.platform=="cluster":
            os.system('export BSUB_QUIET=Y; jobs=( $(ls jobs/edgeR/junctions/*.job) ); g=10; for((i=0; i < ${#jobs[@]}; i+=g)); do part=( "${jobs[@]:i:g}" ); for job_fname in ${part[*]}; do echo "[edgeR.junctions] submitted $job_fname"; bsub -M ' + config.edgeR_memory + ' -K < ${job_fname} & done; wait; echo "[edgeR.junctions] processing next 10"; done; echo "[edgeR.junctions] processing complete"')
        if splicekit.config.platform=="desktop":
            os.system(f". jobs/edgeR/junctions/process.sh")
        if splicekit.config.platform=="SLURM":
            os.system('jobs=( $(ls jobs/edgeR/junctions/*.job) ); g=10; for((i=0; i < ${#jobs[@]}; i+=g)); do part=( "${jobs[@]:i:g}" ); job_ids=(); for job_fname in "${part[@]}"; do echo "[edgeR.junctions] submitted $job_fname"; job_id=$(sbatch --mem=' + config.edgeR_memory + ' --parsable ${job_fname}); job_ids+=($job_id); done; for job_id in "${job_ids[@]}"; do scontrol show job $job_id | grep -q "JobState=COMPLETED" || scontrol wait job $job_id; done; echo "[edgeR.junctions] processing next 10"; done; echo "[edgeR.junctions] processing complete"')
        splicekit.core.report.edgeR_feature('junctions')
        splicekit.core.patterns.process() # adds patterns

    if run=="exons" or run==None:
        os.system(f"rm -f results/edgeR/exons/*.tab.gz > /dev/null 2>&1")
        if splicekit.config.platform=="cluster":
            os.system('export BSUB_QUIET=Y; jobs=( $(ls jobs/edgeR/exons/*.job) ); g=10; for((i=0; i < ${#jobs[@]}; i+=g)); do part=( "${jobs[@]:i:g}" ); for job_fname in ${part[*]}; do echo "[edgeR.exons] submitted $job_fname"; bsub -M ' + config.edgeR_memory + ' -K < ${job_fname} & done; wait; echo "[edgeR.exons] processing next 10"; done; echo "[edgeR.exons] processing complete"')
        if splicekit.config.platform=="desktop":
            os.system(f". jobs/edgeR/exons/process.sh")
        if splicekit.config.platform=="SLURM":
            os.system('jobs=( $(ls jobs/edgeR/exons/*.job) ); g=10; for((i=0; i < ${#jobs[@]}; i+=g)); do part=( "${jobs[@]:i:g}" ); job_ids=(); for job_fname in "${part[@]}"; do echo "[edgeR.exons] submitted $job_fname"; job_id=$(sbatch --mem=' + config.edgeR_memory + ' --parsable ${job_fname}); job_ids+=($job_id); done; for job_id in "${job_ids[@]}"; do scontrol show job $job_id | grep -q "JobState=COMPLETED" || scontrol wait job $job_id; done; echo "[edgeR.exons] processing next 10"; done; echo "[edgeR.exons] processing complete"')
        splicekit.core.report.edgeR_feature('exons')

    if run=="genes" or run==None:
        os.system(f"rm -f results/edgeR/genes/*.tab.gz > /dev/null 2>&1")
        if splicekit.config.platform=="cluster":
            os.system('export BSUB_QUIET=Y; jobs=( $(ls jobs/edgeR/genes/*.job) ); g=10; for((i=0; i < ${#jobs[@]}; i+=g)); do part=( "${jobs[@]:i:g}" ); for job_fname in ${part[*]}; do echo "[edgeR.genes] submitted $job_fname"; bsub -M ' + config.edgeR_memory + ' -K < ${job_fname} & done; wait; echo "[edgeR.genes] processing next 10"; done; echo "[edgeR.exons] processing complete"')
        if splicekit.config.platform=="desktop":
            os.system(f". jobs/edgeR/genes/process.sh")
        if splicekit.config.platform=="SLURM":
            os.system('jobs=( $(ls jobs/edgeR/genes/*.job) ); g=10; for((i=0; i < ${#jobs[@]}; i+=g)); do part=( "${jobs[@]:i:g}" ); job_ids=(); for job_fname in "${part[@]}"; do echo "[edgeR.genes] submitted $job_fname"; job_id=$(sbatch --mem=${config.edgeR_memory} --parsable ${job_fname}); job_ids+=($job_id); done; for job_id in "${job_ids[@]}"; do scontrol show job $job_id | grep -q "JobState=COMPLETED" || scontrol wait job $job_id; done; echo "[edgeR.genes] processing next 10"; done; echo "[edgeR.genes] processing complete"')
        splicekit.core.report.edgeR_feature('genes')

    if run=="anchors" or run==None:
        os.system(f"rm -f results/edgeR/donor_anchors/*.tab.gz > /dev/null 2>&1")
        if splicekit.config.platform=="cluster":
            os.system('export BSUB_QUIET=Y; jobs=( $(ls jobs/edgeR/donor_anchors/*.job) ); g=10; for((i=0; i < ${#jobs[@]}; i+=g)); do part=( "${jobs[@]:i:g}" ); for job_fname in ${part[*]}; do echo "[edgeR.donor_anchors] submitted $job_fname"; bsub -M ' + config.edgeR_memory + ' -K < ${job_fname} & done; wait; echo "[edgeR] processing next 10"; done; echo "[edgeR.donor_anchors] processing complete"')
        if splicekit.config.platform=="desktop":
            os.system(f". jobs/edgeR/donor_anchors/process.sh")
        if splicekit.config.platform=="SLURM":
            os.system('jobs=( $(ls jobs/edgeR/donor_anchors/*.job) ); g=10; for((i=0; i < ${#jobs[@]}; i+=g)); do part=( "${jobs[@]:i:g}" ); job_ids=(); for job_fname in "${part[@]}"; do echo "[edgeR.donor_anchors] submitted $job_fname"; job_id=$(sbatch --mem=${config.edgeR_memory} --parsable ${job_fname}); job_ids+=($job_id); done; for job_id in "${job_ids[@]}"; do scontrol show job $job_id | grep -q "JobState=COMPLETED" || scontrol wait job $job_id; done; echo "[edgeR] processing next 10"; done; echo "[edgeR.donor_anchors] processing complete"')
        splicekit.core.report.edgeR_feature('donor_anchors')
        os.system(f"rm -f results/edgeR/acceptor_anchors/*.tab.gz > /dev/null 2>&1")
        if splicekit.config.platform=="cluster":
            os.system('export BSUB_QUIET=Y; jobs=( $(ls jobs/edgeR/acceptor_anchors/*.job) ); g=10; for((i=0; i < ${#jobs[@]}; i+=g)); do part=( "${jobs[@]:i:g}" ); for job_fname in ${part[*]}; do echo "[edgeR.acceptor_anchors] submitted $job_fname"; bsub -M ' + config.edgeR_memory + ' -K < ${job_fname} & done; wait; echo "[edgeR] processing next 10"; done; echo "[edgeR.acceptor_anchors] processing complete"')
        if splicekit.config.platform=="desktop":
            os.system(f". jobs/edgeR/acceptor_anchors/process.sh")
        if splicekit.config.platform=="SLURM":
            os.system('jobs=( $(ls jobs/edgeR/acceptor_anchors/*.job) ); g=10; for((i=0; i < ${#jobs[@]}; i+=g)); do part=( "${jobs[@]:i:g}" ); job_ids=(); for job_fname in "${part[@]}"; do echo "[edgeR.acceptor_anchors] submitted $job_fname"; job_id=$(sbatch --mem=${config.edgeR_memory} --parsable ${job_fname}); job_ids+=($job_id); done; for job_id in "${job_ids[@]}"; do scontrol show job $job_id | grep -q "JobState=COMPLETED" || scontrol wait job $job_id; done; echo "[edgeR] processing next 10"; done; echo "[edgeR.acceptor_anchors] processing complete"')
        splicekit.core.report.edgeR_feature('acceptor_anchors')

def dexseq(run=None):

    if run=="junctions" or run==None:
        splicekit.core.annotation.make_comparisons()
        os.system(f"rm -f results/dexseq/junctions/*.tab.gz > /dev/null 2>&1")
        if splicekit.config.platform=="cluster":
            os.system('export BSUB_QUIET=Y; jobs=( $(ls jobs/dexseq/junctions/*.job) ); g=10; for((i=0; i < ${#jobs[@]}; i+=g)); do part=( "${jobs[@]:i:g}" ); for job_fname in ${part[*]}; do echo "[dexseq.junctions] submitted $job_fname"; bsub -M ' + config.dexseq_memory + ' -K < ${job_fname} & done; wait; echo "[dexseq.junctions] processing next 10"; done; echo "[dexseq.junctions] processing complete"')
        if splicekit.config.platform=="desktop":
            os.system(f". jobs/dexseq/junctions/process.sh")
        if splicekit.config.platform=="SLURM":
            os.system('jobs=( $(ls jobs/dexseq/junctions/*.job) ); g=10; for((i=0; i < ${#jobs[@]}; i+=g)); do part=( "${jobs[@]:i:g}" ); job_ids=(); for job_fname in "${part[@]}"; do echo "[dexseq.junctions] submitted $job_fname"; job_id=$(sbatch --mem=' + config.dexseq_memory + ' --parsable ${job_fname}); job_ids+=($job_id); done; for job_id in "${job_ids[@]}"; do scontrol show job $job_id | grep -q "JobState=COMPLETED" || scontrol wait job $job_id; done; echo "[dexseq.junctions] processing next 10"; done; echo "[dexseq.junctions] processing complete"')
        splicekit.core.report.dexseq_feature('junctions')
        splicekit.core.patterns.process() # adds patterns

    if run=="exons" or run==None:
        os.system(f"rm -f results/dexseq/exons/*.tab.gz > /dev/null 2>&1")
        if splicekit.config.platform=="cluster":
            os.system('export BSUB_QUIET=Y; jobs=( $(ls jobs/dexseq/exons/*.job) ); g=10; for((i=0; i < ${#jobs[@]}; i+=g)); do part=( "${jobs[@]:i:g}" ); for job_fname in ${part[*]}; do echo "[dexseq.exons] submitted $job_fname"; bsub -M ' + config.dexseq_memory + ' -K < ${job_fname} & done; wait; echo "[dexseq.exons] processing next 10"; done; echo "[dexseq.exons] processing complete"')
        if splicekit.config.platform=="desktop":
            os.system(f". jobs/dexseq/exons/process.sh")
        if splicekit.config.platform=="SLURM":
            os.system('jobs=( $(ls jobs/dexseq/exons/*.job) ); g=10; for((i=0; i < ${#jobs[@]}; i+=g)); do part=( "${jobs[@]:i:g}" ); job_ids=(); for job_fname in "${part[@]}"; do echo "[dexseq.exons] submitted $job_fname"; job_id=$(sbatch --mem=' + config.dexseq_memory + ' --parsable ${job_fname}); job_ids+=($job_id); done; for job_id in "${job_ids[@]}"; do scontrol show job $job_id | grep -q "JobState=COMPLETED" || scontrol wait job $job_id; done; echo "[dexseq.exons] processing next 10"; done; echo "[dexseq.exons] processing complete"')
        splicekit.core.report.dexseq_feature('exons')

    if run=="genes" or run==None:
        os.system(f"rm -f results/dexseq/genes/*.tab.gz > /dev/null 2>&1")
        if splicekit.config.platform=="cluster":
            os.system('export BSUB_QUIET=Y; jobs=( $(ls jobs/dexseq/genes/*.job) ); g=10; for((i=0; i < ${#jobs[@]}; i+=g)); do part=( "${jobs[@]:i:g}" ); for job_fname in ${part[*]}; do echo "[dexseq.genes] submitted $job_fname"; bsub -M ' + config.dexseq_memory + ' -K < ${job_fname} & done; wait; echo "[dexseq.genes] processing next 10"; done; echo "[dexseq.exons] processing complete"')
        if splicekit.config.platform=="desktop":
            os.system(f". jobs/dexseq/genes/process.sh")
        if splicekit.config.platform=="SLURM":
            os.system('jobs=( $(ls jobs/dexseq/genes/*.job) ); g=10; for((i=0; i < ${#jobs[@]}; i+=g)); do part=( "${jobs[@]:i:g}" ); job_ids=(); for job_fname in "${part[@]}"; do echo "[dexseq.genes] submitted $job_fname"; job_id=$(sbatch --mem=${config.dexseq_memory} --parsable ${job_fname}); job_ids+=($job_id); done; for job_id in "${job_ids[@]}"; do scontrol show job $job_id | grep -q "JobState=COMPLETED" || scontrol wait job $job_id; done; echo "[dexseq.genes] processing next 10"; done; echo "[dexseq.genes] processing complete"')
        splicekit.core.report.dexseq_feature('genes')

    if run=="anchors" or run==None:
        os.system(f"rm -f results/dexseq/donor_anchors/*.tab.gz > /dev/null 2>&1")
        if splicekit.config.platform=="cluster":
            os.system('export BSUB_QUIET=Y; jobs=( $(ls jobs/dexseq/donor_anchors/*.job) ); g=10; for((i=0; i < ${#jobs[@]}; i+=g)); do part=( "${jobs[@]:i:g}" ); for job_fname in ${part[*]}; do echo "[dexseq.donor_anchors] submitted $job_fname"; bsub -M ' + config.dexseq_memory + ' -K < ${job_fname} & done; wait; echo "[dexseq] processing next 10"; done; echo "[dexseq.donor_anchors] processing complete"')
        if splicekit.config.platform=="desktop":
            os.system(f". jobs/dexseq/donor_anchors/process.sh")
        if splicekit.config.platform=="SLURM":
            os.system('jobs=( $(ls jobs/dexseq/donor_anchors/*.job) ); g=10; for((i=0; i < ${#jobs[@]}; i+=g)); do part=( "${jobs[@]:i:g}" ); job_ids=(); for job_fname in "${part[@]}"; do echo "[dexseq.donor_anchors] submitted $job_fname"; job_id=$(sbatch --mem=${config.dexseq_memory} --parsable ${job_fname}); job_ids+=($job_id); done; for job_id in "${job_ids[@]}"; do scontrol show job $job_id | grep -q "JobState=COMPLETED" || scontrol wait job $job_id; done; echo "[dexseq] processing next 10"; done; echo "[dexseq.donor_anchors] processing complete"')
        splicekit.core.report.dexseq_feature('donor_anchors')
        os.system(f"rm -f results/dexseq/acceptor_anchors/*.tab.gz > /dev/null 2>&1")
        if splicekit.config.platform=="cluster":
            os.system('export BSUB_QUIET=Y; jobs=( $(ls jobs/dexseq/acceptor_anchors/*.job) ); g=10; for((i=0; i < ${#jobs[@]}; i+=g)); do part=( "${jobs[@]:i:g}" ); for job_fname in ${part[*]}; do echo "[dexseq.acceptor_anchors] submitted $job_fname"; bsub -M ' + config.dexseq_memory + ' -K < ${job_fname} & done; wait; echo "[dexseq] processing next 10"; done; echo "[dexseq.acceptor_anchors] processing complete"')
        if splicekit.config.platform=="desktop":
            os.system(f". jobs/dexseq/acceptor_anchors/process.sh")
        if splicekit.config.platform=="SLURM":
            os.system('jobs=( $(ls jobs/dexseq/acceptor_anchors/*.job) ); g=10; for((i=0; i < ${#jobs[@]}; i+=g)); do part=( "${jobs[@]:i:g}" ); job_ids=(); for job_fname in "${part[@]}"; do echo "[dexseq.acceptor_anchors] submitted $job_fname"; job_id=$(sbatch --mem=${config.dexseq_memory} --parsable ${job_fname}); job_ids+=($job_id); done; for job_id in "${job_ids[@]}"; do scontrol show job $job_id | grep -q "JobState=COMPLETED" || scontrol wait job $job_id; done; echo "[dexseq] processing next 10"; done; echo "[dexseq.acceptor_anchors] processing complete"')
        splicekit.core.report.dexseq_feature('acceptor_anchors')

def juan():
    splicekit.core.juan.append_results() # reads in results_edgeR_junctions.tab and appends anchor info from results/results_edgeR_anchors/{comparison}_altsplice.tab

def promisc():
    splicekit.core.promisc.process()

def dreme():
    import splicekit.core.motifs
    splicekit.core.motifs.dreme()

def motifs():
    import splicekit.core.motifs
    os.system("rm -f results/motifs/* > /dev/null 2>&1")
    os.system("rm -rf results/motifs/dreme/* > /dev/null 2>&1")
    os.system("rm -rf results/motifs/fasta/* > /dev/null 2>&1")
    os.system("rm -f results/motifs/scanRBP/* > /dev/null 2>&1")
    os.system("rm -f results/motifs/scanRBP/data/* > /dev/null 2>&1")
    os.system("rm -f results/motifs/scanRBP/fasta/* > /dev/null 2>&1")
    splicekit.core.motifs.process()

def judge_process():
    splicekit.judge.process()

def clusterlogfc_process():
    splicekit.clusterlogfc.process()

def jbrowse2_process(force_samples=False, force_annotation=False):
    splicekit.core.jbrowse2.process(force_samples, force_annotation)

def rmats():
    if splicekit.config.platform=="cluster":
        os.system('export BSUB_QUIET=Y; jobs=( $(ls jobs/rmats/*.job) ); g=10; for((i=0; i < ${#jobs[@]}; i+=g)); do part=( "${jobs[@]:i:g}" ); for job_fname in ${part[*]}; do echo "[splicekit.rmats] submitted $job_fname"; bsub -M 8GB -K < ${job_fname} & done; wait; echo "[splicekit.rmats] processing next 10"; done; echo "[splicekit.rmats] processing complete"')
    if splicekit.config.platform=="desktop":
        splicekit.core.mprocess("jobs/rmats/process.sh")
    if splicekit.config.platform=="SLURM":
        os.system('jobs=( $(ls jobs/rmats/*.job) ); g=10; for((i=0; i < ${#jobs[@]}; i+=g)); do part=( "${jobs[@]:i:g}" ); job_ids=(); for job_fname in "${part[@]}"; do echo "[splicekit.rmats] submitted $job_fname"; job_id=$(sbatch --mem=8G --parsable ${job_fname}); job_ids+=($job_id); done; for job_id in "${job_ids[@]}"; do scontrol show job $job_id | grep -q "JobState=COMPLETED" || scontrol wait job $job_id; done; echo "[splicekit.rmats] processing next 10"; done; echo "[splicekit.rmats] processing complete"')

def june_process():
    splicekit.june.process()

def process(force=False):
    setup()
    annotation()
    features()
    edgeR()
    juan()
    judge_process()
    motifs()
    promisc()
    clusterlogfc_process()
    june_process()
    splicekit.report.process()
    jbrowse2_process(force_samples=force, force_annotation=force)

def main():
    help_0 = """
    Usage: splicekit <command> [options]

    Commands:
    -- Processing
        process           run all available analyses

    -- Initialization
        setup             initialize folder structure
        annotation        download sample annotation from Moose
        features          create junction/exon/anchor count tables from .bam files
    
    -- Splicing analyses
        edgeR             run edgeR analyses on junctions, exons and anchors
        judge             run juDGE plot analysis
        promisc           promiscuity analysis: genes and their junctions changes across conditions
        clusterlogfc      log fold-change clusters of significant changes (FDR<0.05) on junction, exon, gene level
        rmats             process data with rMATS
        june              junction-event (junE) analysis

    -- scanRBP, motif and sequence analyses
        motifs            run motif and scanRBP analysis

    -- JBrowse2 genomic browser
        jbrowse2 [process | start] processes (creates) all required JBrowse2 files and/or starts JBrowse2 web server

    Options

    -- hostname
        If you are submitting 'splicekit process' to the cluster, you can provide a custom hostname. This will be used to create JBrowse2 URLs.

    """

    help_edgeR = """
    Usage: splicekit edgeR [sub-command]

    Sub-commands:
        junctions      run edgeR on junctions table data
        exons          run edgeR on exons table data
        anchors        run edgeR on anchors table data
        genes          run edgeR on genes table data

    If no sub-command is given, all analyses will be run (junctions, exons, anchors, genes).
    """

    help_dexseq = """
    Usage: splicekit dexseq [sub-command]

    Sub-commands:
        junctions      run DEXSeq on junctions table data
        exons          run DEXSeq on exons table data
        anchors        run DEXSeq on anchors table data
        genes          run DEXSeq on genes table data

    If no sub-command is given, all analyses will be run (junctions, exons, anchors, genes).
    """

    help_motifs = """
    Usage: splicekit motifs [sub-command]

    Sub-commands:
        dreme          run DREME on produced FASTA files

    If no sub-command is given, all analyses will be run (DREME, scanRBP).
    """

    # lookup version without loading splicekit
    import importlib.util
    splicekit_init = importlib.util.find_spec("splicekit")
    splicekit_path = splicekit_init.origin
    splicekit_folder = os.path.dirname(splicekit_path)
    version = open(os.path.join(splicekit_folder, "version"), "rt").readlines()[0].replace("\n", "").replace("\r", "")
    verbose = args.verbose

    print(f"splicekit | v{version} | https://github.com/bedapub/splicekit")

    if args.version:
        sys.exit(0)

    if len(args.command)==0 and args.help:
        print(help_0)
        sys.exit(0)

    if len(args.command)==0 and not args.help:
        print()
        print("splicekit | please provide an action to perform or specify --help to get help on available actions; running 'splicekit process' will perform all available analyses")
        sys.exit(0)

    # loading splicekit modules takes time, only load if the user is not asking to display help
    if not args.help:
        import splicekit
        if args.hostname!=None:
            splicekit.config.jbrowse2_config(args.hostname)
        else:
            splicekit.config.jbrowse2_config(None)

    if args.command[0]=="basic":
        splicekit.setup()
        splicekit.annotation()
        splicekit.features()
        splicekit.edgeR()
    elif args.command[0]=="setup":
        splicekit.setup()
    elif args.command[0]=="annotation":
        splicekit.annotation()
    elif args.command[0]=="junctions":
        splicekit.junctions()
    elif args.command[0]=="junctions_master":
        splicekit.junctions_master()
    elif args.command[0]=="exons":
        splicekit.exons()
    elif args.command[0]=="genes":
        splicekit.genes()
    elif args.command[0]=="anchors":
        splicekit.anchors()
    elif args.command[0]=="features":
        splicekit.features()
    elif args.command[0]=="edgeR":
        if args.help:
            print(help_edgeR)
            sys.exit(0)
        if len(args.command)>1 and args.command[1] in ["junctions", "exons", "anchors", "genes"]:
            splicekit.core.features.load_genes()
            splicekit.edgeR(run=args.command[1])
        elif len(args.command)==1:
            splicekit.edgeR()
        else:
            print(help_edgeR)
            sys.exit(0)
    elif args.command[0]=="dexseq":
        if args.help:
            print(help_dexseq)
            sys.exit(0)
        if len(args.command)>1 and args.command[1] in ["junctions", "exons", "anchors", "genes"]:
            splicekit.core.features.load_genes()
            splicekit.dexseq(run=args.command[1])
        elif len(args.command)==1:
            splicekit.dexseq()
        else:
            print(help_dexseq)
            sys.exit(0)
    elif args.command[0]=="motifs":
        if args.help:
            print(help_motifs)
            sys.exit(0)
        if len(args.command)>1 and args.command[1]=="dreme":
            splicekit.dreme()
        else:
            splicekit.motifs()
    elif args.command[0]=="promisc":
        splicekit.promisc()
    elif args.command[0]=="cassettes":
        splicekit.cassettes()
    elif args.command[0]=="patterns":
        splicekit.patterns()
    elif args.command[0]=="clusterlogfc":
        splicekit.clusterlogfc_process()
    elif args.command[0]=="process":
        splicekit.process(force=args.force)
    elif args.command[0]=="judge":
        splicekit.judge_process()
    elif args.command[0]=="juan":
        splicekit.juan()
    elif args.command[0]=="june":
        splicekit.june_process()
    elif args.command[0]=="rmats":
        splicekit.rmats()
    elif args.command[0]=="web":
        splicekit.core.jbrowse2.start()
    elif args.command[0] in ["jbrowse2", "jbrowse"]:
        if len(args.command)>1:
            if args.command[1]=="process":
                splicekit.jbrowse2_process(force_samples=args.force, force_annotation=args.force)
            if args.command[1]=="start":
                splicekit.core.jbrowse2.start()
        else:
            splicekit.jbrowse2_process(force_samples=args.force, force_annotation=args.force)
            splicekit.core.jbrowse2.start()
    elif args.command[0]=="report":
        splicekit.report.process()
    else:
        print(help_0)

def splicecompare():
    import seaborn as sns
    import matplotlib.pyplot as plt
    import gzip

    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('splicekit_run1', help="Main folder with splicekit run 1")
    parser.add_argument('comparison1', help="Name of comparison to use in splicekit run 1")
    parser.add_argument('comparison1_fname', help="")
    parser.add_argument('splicekit_run2', help="Main folder with splicekit run 2")
    parser.add_argument('comparison2', help="Name of comparison to use in splicekit run 2")
    parser.add_argument('comparison2_fname', help="")
    parser.add_argument('fname_out', help="")
    parser.add_argument("-fdr", "--fdr", help="FDR threshold, default no filter", default=None, type=float)
    parser.add_argument("-logFC", "--logFC", help="logFC threshold, default=None (do not use as filter)", default=None, type=float)
    parser.add_argument("-dpfi", "--dpfi", help="delta(percentage feature inclusion) threshold, default=None (do not use as filter)", default=None, type=float)
    parser.add_argument("-intersect", "--intersect", help="Report only features that are hit in both comparisons", action="store_true")
    parser.add_argument("-features", "--features", help="Which features to process, comma separated (default: 'exons,junctions,donor_anchors,acceptor_anchors')", default="exons,junctions,donor_anchors,acceptor_anchors")
    args = parser.parse_args()

    def print_version():
        print("splicecompare v0.3")
        print("FDR threshold = ", args.fdr)
        print("logFC threshold = ", args.logFC)
        print("dpfi threshold = ", args.dpfi)
        print("----")

    print_version()

    fdr_key = None
    data_all = {}
    for feature_type in args.features.split(","):
        fout = open(f"{args.fname_out}_{feature_type}.tab", "wt")
        results = {}
        for splicekit_run, comparison, comparison_fname in [(args.splicekit_run1, args.comparison1, args.comparison1_fname), (args.splicekit_run2, args.comparison2, args.comparison2_fname)]:
            f = gzip.open(os.path.join(splicekit_run, "results", "edgeR", feature_type, f"{comparison_fname}.gz"), "rt")
            header = f.readline().replace("\r", "").replace("\n", "").split("\t")
            r = f.readline()
            while r:
                r = r.replace("\r", "").replace("\n", "").split("\t")
                data = dict(zip(header, r))
                # store all the data
                temp = data_all.get(comparison, {})
                temp[data["feature_id"]] = data
                data_all[comparison] = temp
                if data.get("fdr", None) != None:
                    fdr_key = "fdr"
                if data.get("FDR", None) != None:
                    fdr_key = "FDR"
                if args.fdr!=None:
                    if fdr_key!=None:
                        if float(data[fdr_key])>args.fdr:
                            r = f.readline()
                            continue
                if args.logFC!=None:
                    if abs(float(data["logFC"]))<args.logFC:
                        r = f.readline()
                        continue
                if args.dpfi!=None:
                    if abs(float(data["delta_pfi"]))<args.dpfi:
                        r = f.readline()
                        continue
                results_comparison = results.get(comparison, {})
                results_comparison[data["feature_id"]] = data
                results[comparison] = results_comparison
                r = f.readline()
            f.close()
        
        list1 = list(results.get(args.comparison1, {}).keys())
        list2 = list(results.get(args.comparison2, {}).keys())

        fout.write("# " + str(len(list1)) + f" {feature_type} detected for {args.comparison1}" + "\n")
        fout.write("# " + str(len(list2)) + f" {feature_type} detected for {args.comparison2}" + "\n")
        fout.write("# intersect = " + str(len(list(set(list2).intersection(set(list1))))) + "\n")

        feature_keys = list(set(list(results.get(args.comparison1, {}).keys()) + list(results.get(args.comparison2, {}).keys())))
        header = ["feature_id", "chr", "strand", "gene_name", f"fdr_{args.comparison1}", f"logFC_{args.comparison1}", f"delta_pfi_{args.comparison1}", f"fdr_{args.comparison2}", f"logFC_{args.comparison2}", f"delta_pfi_{args.comparison2}"]
        fout.write("\t".join(header) + "\n")

        final_results = []
        feature_intersect = list(set(list2).intersection(set(list1))) # list of features for which both comparison1 and comparison2 satisfy condition (FDR<thr, logFC<thr etc)
        for feature_id in feature_keys:
            row = [feature_id]
            data1 = data_all.get(args.comparison1, {}).get(feature_id, {})
            data2 = data_all.get(args.comparison2, {}).get(feature_id, {})
            if args.intersect and feature_id not in feature_intersect:
                continue
            if data1=={} or data2=={}:
                continue
            data_common = data1 if data1!={} else data2
            row.append(data_common["chr"])
            row.append(data_common["strand"])
            row.append(data_common["gene_name"])
            
            for data in [data1, data2]:
                row.append(data.get(fdr_key, ""))
                row.append(data.get("logFC", ""))
                row.append(data.get("delta_pfi", ""))

            combined_logFC = abs(float(data1.get("logFC", 0))) + abs(float(data2.get("logFC", 0)))
            final_results.append((combined_logFC, row))
        final_results.sort(reverse=True)
        for combined_logFC, row in final_results:
            fout.write("\t".join(row) + "\n")
        fout.close()    