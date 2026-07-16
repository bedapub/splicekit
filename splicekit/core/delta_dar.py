"""
# Description
Reads in results_edgeR_jucntions.tab and appends three columns: 
    donor_acceptor_ratio (dar_i) = 100* donor_i/(donor_i + acceptor_i) for i = test
    donor_acceptor_ratio (dar_i) = 100* donor_i/(donor_i + acceptor_i) for i = control
    delta_dar = dar_test - dar_control
"""

import os
import sys
import pandas as pd

# define paths, read in comparison and results file
sample_annotation_dir =  'annotation' 
results_dir = 'results'
donor_data_dir = "data/comparison_donor_anchors_data/"
acceptor_data_dir = "data/comparison_acceptor_anchors_data/"
comparison_fname = '{comparison}.tab'

def compute():

    comparisons_df = pd.read_csv(f"{sample_annotation_dir}/comparisons.tab", sep = '\t')
    results_df = pd.read_csv(f"{results_dir}/results_edgeR_junctions.tab", sep = '\t')

    results_df['DAR_test'] = None
    results_df['DAR_control'] = None
    results_df['delta_DAR'] = None

    results_df.index = results_df.result_id
    comparisons_df.index = comparisons_df.comparison

    # do calculations for the entire results table --------------------------------------------------------
    comparisons = list(set(results_df.comparison))

    for comp in comparisons:

        # focus on comparison only
        results_df_subset = results_df[ (results_df.comparison == comp)].copy()

        # read in feature data
        donor_df = pd.read_csv(donor_data_dir+ comparison_fname.format(comparison=comp), sep = '\t')
        acceptor_df = pd.read_csv(acceptor_data_dir+ comparison_fname.format(comparison=comp), sep = '\t')

        # iterate over detected junctions in comparison and compute delta_pct_anchors
        for id in results_df_subset.index:
            
            # get detected junction
            junction = results_df_subset.loc[id,]

            # get number of samples per condition 
            number_test_samples = len(comparisons_df.loc[comp, 'compound_samples'].split(','))
            number_control_samples = len(comparisons_df.loc[comp, 'DMSO_samples'].split(','))
            
            # for a given junction we have to identify the donor and acceptor exon
            strand = junction['strand']
            chr = junction['chr']
            gene_id = junction['gene_id']
            gene_name = junction['gene_name']

            feature_start = junction['feature_start']
            feature_stop = junction['feature_stop']

            if strand == '+':
                anchor1_df = donor_df[(donor_df.chr==chr)&(donor_df.gene_id==gene_id)&(donor_df.feature_stop==feature_start)]
                anchor2_df = acceptor_df[(acceptor_df.chr==chr)&(acceptor_df.gene_id==gene_id)&(acceptor_df.feature_start==feature_stop)]
            else:
                anchor1_df = donor_df[(donor_df.chr==chr)&(donor_df.gene_id==gene_id)&(donor_df.feature_start==feature_stop)]
                anchor2_df = acceptor_df[(acceptor_df.chr==chr)&(acceptor_df.gene_id==gene_id)&(acceptor_df.feature_stop==feature_start)]


            # actual computation
            anchors_cov_test = (anchor1_df.sum_feature_test.to_list()[0]/number_test_samples,anchor2_df.sum_feature_test.to_list()[0]/number_test_samples)
            anchors_cov_control = (anchor1_df.sum_feature_control.to_list()[0]/number_control_samples, anchor2_df.sum_feature_control.to_list()[0]/number_control_samples)


            if sum(anchors_cov_test) != 0:
                dar_test = 100 * anchors_cov_test[0]/sum(anchors_cov_test)
            else:
                dar_test = 0
            if sum(anchors_cov_control) != 0:
                dar_control = 100 * anchors_cov_control[0]/sum(anchors_cov_control)
            else:
                dar_control = 0

            delta_dar = dar_test - dar_control
            results_df.loc[id,['delta_DAR']] = delta_dar
            results_df.loc[id,['DAR_control']] = dar_control
            results_df.loc[id,['DAR_test']] = dar_test 

    results_df.to_csv(f"{results_dir}/results_edgeR_junctions.tab", sep='\t', index=False)

