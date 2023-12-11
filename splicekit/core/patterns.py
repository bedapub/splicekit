import os
import sys
import gzip
import splicekit.config as config
import pybio
import copy

pattern_area = (-2, 6)

def process(version=""):
    def process_file(fname):
        fin = gzip.open(f"results/{fname}.tab.gz", "rt")
        header = fin.readline().replace("\r", "").replace("\n", "").split("\t")
        header_out = header.copy()
        if "donor_pattern" not in header_out:
            header_out.append("donor_pattern")
        if "acceptor_pattern" not in header_out:
            header_out.append("acceptor_pattern")
        fout = gzip.open(f"results/{fname}_new.tab.gz", "wt")
        fout.write("\t".join(header_out)+"\n")
        r = fin.readline()
        while r:
            r = r.replace("\r", "").replace("\n", "").split("\t")
            data_out = dict(zip(header_out, r))
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
            fout.write("\t".join(str(data_out[h]) for h in header_out) + "\n")
            r = fin.readline()
        fout.close()
        os.system(f"mv results/{fname}_new.tab.gz results/{fname}.tab.gz")
        fin.close()
    process_file(f"results_edgeR{version}_junctions")
    process_file(f"results_edgeR{version}_junctions_all")
