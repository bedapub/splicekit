import os
import sys
import gzip
import splicekit.config as config
import pybio
import copy

pattern_area = (-2, 6)

def process():
    def process_file(fname):
        fin = gzip.open(f"results/edgeR/{fname}.tab.gz", "rt")
        header = fin.readline().replace("\r", "").replace("\n", "").split("\t")
        header_new = header.copy()
        for el in ["donor_pattern", "acceptor_pattern"]:
            if el not in header:
                header_new.append(el)
        fout = gzip.open(f"results/{fname}.tab.gz.temp", "wt")
        fout.write("\t".join(header_new)+"\n")
        r = fin.readline()
        while r:
            r = r.replace("\r", "").replace("\n", "").split("\t")
            data_out = dict(zip(header_new, r))
            coords = data_out["feature_id"].split('_')
            start = int(coords[-2])
            stop = int(coords[-1])
            strand = coords[-3][-1]
            chr = '_'.join(coords[:-2])[:-1]
            if strand=="+":
                donor_site, acceptor_site = start, stop
            else:
                donor_site, acceptor_site = stop, start
            donor_seq = pybio.core.genomes.seq(config.species, chr, strand, donor_site, pattern_area[0], pattern_area[1])
            acceptor_seq = pybio.core.genomes.seq(config.species, chr, strand, acceptor_site, pattern_area[0], pattern_area[1])
            data_out["donor_pattern"] = donor_seq
            data_out["acceptor_pattern"] = acceptor_seq
            fout.write("\t".join(str(data_out[h]) for h in header_new) + "\n")
            r = fin.readline()
        fout.close()
        os.system(f"mv results/edgeR/{fname}tab.gz.temp results/edgeR/{fname}.tab.gz")
        fin.close()
    process_file(f"junctions_results")
    process_file(f"junctions_results_complete")
