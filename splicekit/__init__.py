import os
import sys
import splicekit.core as core

# load version
splicekit_path = os.path.abspath(__file__)
splicekit_folder = os.path.dirname(splicekit_path)
version = open(os.path.join(splicekit_folder, "version"), "rt").readlines()[0].replace("\n", "").replace("\r", "")

if not os.path.exists("splicekit.config"):
    print("splicekit | please run splicekit in a folder with splicekit.config present")
    sys.exit(0)

print("splicekit | loading splicekit.config")
import splicekit.config as config
print("splicekit | loading splicekit.core")
import splicekit.core as core
print("splicekit | loading splicekit.core.annotation")
import splicekit.core.annotation
print("splicekit | loading splicekit.core.features")
import splicekit.core.features
print("splicekit | loading splicekit.core.exons")
import splicekit.core.exons
print("splicekit | loading splicekit.core.genes")
import splicekit.core.genes
print("splicekit | loading splicekit.core.report")
import splicekit.core.report
print("splicekit | loading splicekit.core.patterns")
import splicekit.core.patterns
print("splicekit | loading splicekit.core.promisc")
import splicekit.core.promisc
print("splicekit | loading splicekit.core.anchors ")
import splicekit.core.anchors 
print("splicekit | loading splicekit.core.junctions")
import splicekit.core.junctions
print("splicekit | loading splicekit.core.jbrowse2")
import splicekit.core.jbrowse2
print("splicekit | loading splicekit.core.juan")
import splicekit.core.juan
print("splicekit | loading splicekit.judge")
import splicekit.judge
print("splicekit | loading splicekit.core.delta_dar")
import splicekit.core.delta_dar
print("splicekit | loading splicekit.clusterlogfc")
import splicekit.clusterlogfc
print("splicekit | loading splicekit.report")
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
        os.system("jobs=( $(ls jobs/count_junctions/*.job) ); g=10; " "for((i=0; i < ${#jobs[@]}; i+=g)); do " "part=( \"${jobs[@]:i:g}\" ); " "for job_fname in ${part[*]}; do " "echo \"splicekit | features | junctions | submitted $job_fname\"; " "sbatch --mem=8G --parsable ${job_fname} & " "done; wait; " "echo \"splicekit | features | junctions | processing next 10\"; " "done; " "echo \"splicekit | features | junctions | processing complete\"")
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
    jbrowse2_process(force_samples=force, force_annotation=force)
