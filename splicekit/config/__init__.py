# config module of splicekit

import os
import sys
import splicekit
import splicekit.core as core
import socket

module_desc = "splicekit | config |"
splicekit.verbose and print(f"{module_desc} loading")

if os.path.exists("splicekit.config"):
    config_lines = open("splicekit.config").readlines()
    for cline in config_lines:
        exec(cline.replace("\r", "").replace("\n", ""))

jbrowse2_port = 8007
jbrowse2_url = None

try:
    hostname
    ip_addr=socket.gethostbyname(hostname)
except:
    hostname="localhost"
    ip_addr=socket.gethostbyname(hostname)

def jbrowse2_config():
    global jbrowse2_url
    port = jbrowse2_port
    jbrowse2_url = f"http://{ip_addr}:{port}/jbrowse2/?config=splicekit_data/config.json"
    splicekit.verbose and print(f"{module_desc} JBrowse2 URL = {jbrowse2_url}")

# read in location of gtf and fasta files
if genome_version!=None:
    temp = os.popen(f"pybio path {species} -genome_version {genome_version}").read()
else:
    temp = os.popen(f"pybio path {species}").read()
temp = temp.split("\n")
for line in temp:
    if line.endswith(".fasta"):
        fasta_path = line
    if line.endswith(".gtf.gz"):
        gtf_path = line
    if line.endswith(".gff3.gz"):
        gff3_path = line

# memory parameters
        
try:
    edgeR_memory
except:
    edgeR_memory = "8GB"

try:
    dexseq_memory
except:
    dexseq_memory = "8GB"

try:
    dexseq_scripts
except:
    dexseq_scripts = ""

try:
    dexseq_FDR_thr
except:
    dexseq_FDR_thr = 0.05

try:
    cluster_queue
except:
    cluster_queue = "short"

try:
    clip
except:
    clip = None

try:
    bam_path
    bam_path_defined = True
except NameError:
    bam_path_defined = False

try:
    bam_column
    bam_column_defined = True
except NameError:
    bam_column_defined = False

if not bam_path_defined and not bam_column_defined:
    print(f"{module_desc} ERROR: neither bam_path nor bam_column is set in splicekit.config")
    print(f"{module_desc} Set at least one of the following in splicekit.config:")
    print(f"{module_desc}   bam_path = \"/path/to/bam/files\"")
    print(f"{module_desc}     BAMs are then expected at {{bam_path}}/{{sample_id}}.bam")
    print(f"{module_desc}   bam_column = \"bam_file\"")
    print(f"{module_desc}     samples.tab must then have a column with this name containing the full BAM path for each sample")
    sys.exit(1)

try:
    bam_column
except:
    bam_column = "bam_file"