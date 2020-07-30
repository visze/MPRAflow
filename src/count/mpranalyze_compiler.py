#contributed by Tal Ashuach

import os
import re
import sys
import gzip
import pandas as pd
import numpy as np
import click
from collections import defaultdict

# options
@click.command()
@click.option('--input',
              'input_file',
              required=True,
              type=click.Path(exists=True, readable=True),
              help='Merged DNA and RNA counts of all replicates.')
@click.option('--rna-counts-output',
              'rna_counts_output_file',
              required=True,
              type=click.Path(writable=True),
              help='RNA counts output for MPRAnalyze')
@click.option('--dna-counts-output',
              'dna_counts_output_file',
              required=True,
              type=click.Path(writable=True),
              help='DNA counts output for MPRAnalyze.')
@click.option('--rna-annotation-output',
              'rna_annotation_output_file',
              required=True,
              type=click.Path(writable=True),
              help='RNA annotation output for MPRAnalyze.')
@click.option('--dna-annotation-output',
              'dna_annotation_output_file',
              required=True,
              type=click.Path(writable=True),
              help='DNA annotation output for MPRAnalyze.')



def cli(input_file, rna_counts_output_file, dna_counts_output_file, rna_annotation_output_file, dna_annotation_output_file):

    annot_pattern = re.compile("^([DR]NA).*\(condition (.*), replicate (.*)\)$")
    def get_annot(head):
        m = annot_pattern.match(head)
        if m is not None:
            return m.group(1,2,3)

    # read input
    df = pd.read_csv(input_file,sep="\t", header='infer')

    # remove unknown BC
    df = df[df['label']!='no_BC'].set_index('label').fillna(0)

    # get replicates and DNA/RNAs
    header_annot = [get_annot(h) for h in list(df.columns)[2:]]
    dna_annot = [h for h in header_annot if h[0] == 'DNA']
    rna_annot = [h for h in header_annot if h[0] == 'RNA']
    n_dna_obs = len(dna_annot)
    n_rna_obs = len(rna_annot)

    # counts for observation
    # rna_dict = defaultdict(list)
    # dna_dict = defaultdict(list)

    dna_df = df.iloc[:,2:(2+n_dna_obs)].applymap(np.int64)
    rna_df = df.iloc[:,(2+n_dna_obs):].applymap(np.int64)


    n_bc = df.groupby('label').Barcode.agg(len).max()

    ## write output DNA annotations
    dna_colnames = []
    with gzip.open(dna_annotation_output_file, 'wt') as ofile:
        ofile.write('\t'.join(["sample", "type", "condition", "replicate", "barcode"]) + '\n')
        for i in range(1, n_bc + 1):
            for x in dna_annot:
                sample_name = '_'.join(list(x) + [str(i)])
                dna_colnames.append(sample_name)
                ofile.write('\t'.join([sample_name] + list(x) + [str(i)]) + '\n')

    ## write output DNA counts
    with gzip.open(dna_counts_output_file, 'wt') as ofile:
        ofile.write('\t'.join(['seq_id'] + dna_colnames) + '\n')
        for seq_id in dna_df.index.unique():
            ofile.write('\t'.join(
                [seq_id] + # sequence ID
                list(map(str,dna_df.loc[seq_id].values.flatten())) + # flattened list of observations
                (['0'] * n_dna_obs * (n_bc - dna_df.loc[seq_id].index.size)) # zero padding
                ) + '\n')

    ## write output RNA annotations
    rna_colnames = []
    with gzip.open( rna_annotation_output_file, 'wt') as ofile:
        ofile.write('\t'.join(["sample", "type", "condition", "replicate", "barcode"]) + '\n')
        for i in range(1, n_bc + 1):
            for x in rna_annot:
                sample_name = '_'.join(list(x) + [str(i)])
                rna_colnames.append(sample_name)
                ofile.write('\t'.join([sample_name] + list(x) + [str(i)]) + '\n')

    ## write output RNA counts
    with gzip.open(rna_counts_output_file, 'wt') as ofile:
        ofile.write('\t'.join(['seq_id'] + rna_colnames) + '\n')
        for seq_id in rna_df.index.unique():
            ofile.write('\t'.join(
                [seq_id] + # sequence ID
                list(map(str,rna_df.loc[seq_id].values.flatten())) + # flattened list of observations
                (['0'] * n_dna_obs * (n_bc - rna_df.loc[seq_id].index.size)) # zero padding
                ) + '\n')

if __name__ == '__main__':
    cli()
