import os
import sys
import gzip
import splicekit.config as config
import pybio
import copy

module_name = "splicekit | patterns |"

pattern_area = (-2, 6)

def process():
    def process_file(fname):
        print(f"{module_name} +patterns at file {fname}")
        fin = gzip.open(f"results/edgeR/{fname}.tab.gz", "rt")
        header = fin.readline().replace("\r", "").replace("\n", "").split("\t")
        header_new = header.copy()
        for el in ["donor_pattern", "acceptor_pattern"]:
            if el not in header:
                header_new.append(el)
        fout = gzip.open(f"results/edgeR/{fname}.tab.gz.temp", "wt")
        fout.write("\t".join(header_new)+"\n")
        r = fin.readline()
        while r:
            r = r.replace("\r", "").replace("\n", "").split("\t")
            data_new = dict(zip(header_new, r))
            coords = data_new["feature_id"].split('_')
            start = int(coords[-2])
            stop = int(coords[-1])
            strand = coords[-3][-1]
            chr = '_'.join(coords[:-2])[:-1]
            if strand=="+":
                donor_site, acceptor_site = start, stop
            else:
                donor_site, acceptor_site = stop, start
            donor_seq = pybio.core.genomes.seq(config.species, chr, strand, donor_site, pattern_area[0], pattern_area[1], genome_version=config.genome_version)
            acceptor_seq = pybio.core.genomes.seq(config.species, chr, strand, acceptor_site, pattern_area[0], pattern_area[1], genome_version=config.genome_version)
            data_new["donor_pattern"] = donor_seq
            data_new["acceptor_pattern"] = acceptor_seq
            fout.write("\t".join(str(data_new[h]) for h in header_new) + "\n")
            r = fin.readline()
        fout.close()
        fin.close()
        os.system(f"mv results/edgeR/{fname}.tab.gz.temp results/edgeR/{fname}.tab.gz")
    process_file(f"junctions_results_fdr005")
    process_file(f"junctions_results_complete")
