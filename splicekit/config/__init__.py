import os
import sys
import splicekit.core as core
import splicekit.core.jbrowse2 as jbrowse2
import socket

if os.path.exists("splicekit.config"):
    config_lines = open("splicekit.config").readlines()
    for cline in config_lines:
        exec(cline.replace("\r", "").replace("\n", ""))

jbrowse2_url = ""

def jbrowse2_config(hostname):
    global jbrowse2_url
    # JBrowse2 URL
    port = jbrowse2.port
    if hostname==None: # get hostname?
        hostname=socket.gethostname()
    ip_addr=socket.gethostbyname(hostname)
    jbrowse2_url = f"http://{ip_addr}:{port}/jbrowse2/?config=splicekit_data/config.json"

jbrowse2_config(None)

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

