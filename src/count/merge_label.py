#Author: Gracie Gordon 2019
#merge tables

import sys
import pandas as pd
import numpy as np

import math
import pickle

from Bio import SeqIO

import click


# options
@click.command()
@click.option('--counts',
              'counts_file',
              required=True,
              type=click.Path(exists=True, readable=True),
              help='Merged DNA and RNA counts.')
@click.option('--assignment',
              'assignment_file',
              required=True,
              type=click.Path(exists=True, readable=True),
              help='Assignment pickle file.')
@click.option('--design',
              'design_file',
              required=True,
              type=click.Path(exists=True, readable=True),
              help='Design file in fasta format.')
@click.option('--minRNACounts',
              'minRNACounts',
              default=1,
              show_default=True,
              type=int,
              help='Minimum number of RNA counts required per barcode. If 0 pesudocounts are used.')
@click.option('--minDNACounts',
              'minDNACounts',
              default=0,
              show_default=True,
              type=int,
              help='Minimum number of DNA counts required per barcode. If 0 pesudocounts are used.')
@click.option('--output',
              'output_file',
              required=True,
              type=click.Path(writable=True),
              help='Output file.')
def cli(counts_file, assignment_file, design_file, minRNACounts, minDNACounts, output_file):
    pseudocount = 1 if minRNACounts == 0 or minDNACounts == 0 else 0

    # process fastq
    design=open(design_file)
    fasta_dict = {rec.id : rec.seq for rec in SeqIO.parse(design, "fasta")}

    assoc=pickle.load(open(assignment_file,'rb'))

    BC_key = {}
    for k,v in assoc.items():
        for x in v:
             BC_key.setdefault(x,k)


    #get count df
    counts=pd.read_csv(counts_file, sep='\t', header=None,names=['Barcode','dna_count','rna_count'])

    #fill in labels from dictionary
    label=[]
    for i in counts.Barcode:
       try:
               label.append(BC_key[i])
       except:
               label.append('no_BC')

    # INFO: sequence is not needed becaus enot printed out
    # seqs=[]
    # for l in label:
    #     try:
    #         seqs.append(str(fasta_dict[l]).upper())
    #     except:
    #         seqs.append('NA')
    # counts.insert(0,'Sequence',seqs)

    # Insert the label of oligo
    counts.insert(0, 'name', label)

    # remove Barcorde. Not needed anymore
    counts.drop(['Barcode'], axis=1, inplace=True)

    # number of DNA and RNA counts
    total_dna_counts=sum(counts['dna_count'])
    total_rna_counts=sum(counts['rna_count'])

    # filter
    counts = counts[counts['dna_count']>=minDNACounts]
    counts = counts[counts['rna_count']>=minRNACounts]



    grouped_label = counts.groupby(counts['name']).agg({'dna_count' : ['sum', 'count'], 'rna_count' : ['sum', 'count']})

    grouped_label.reset_index(inplace=True)

    output = pd.DataFrame()

    click.echo(grouped_label.head())

    output['name'] = grouped_label['name']
    output['dna_counts'] = grouped_label.dna_count['sum']
    output['rna_counts'] = grouped_label.rna_count['sum']



    scaling = 10**min([len(str(total_dna_counts))-1,len(str(total_rna_counts))-1])
    scaling = 10**6

    output['dna_normalized']=(grouped_label.dna_count['sum']+pseudocount)/((grouped_label.dna_count['count']+pseudocount))/total_dna_counts*scaling



    output['rna_normalized']=(grouped_label.rna_count['sum']+pseudocount)/((grouped_label.rna_count['count']+pseudocount))/total_rna_counts*scaling

    output['ratio'] = output.rna_normalized/output.dna_normalized
    output['log2'] = np.log2(output.ratio)

    output['n_obs_bc'] = grouped_label.dna_count['count']


    click.echo(output_file)

    output.to_csv(output_file, index=False,sep='\t', compression='gzip')

if __name__ == '__main__':
    cli()
