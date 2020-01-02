[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A519.0.7-brightgreen.svg)](https://www.nextflow.io/)

# MPRAflow

This pipeline processes sequencing data from Massively Parallel Reporter Assays (MPRA) to create count tables for candidate sequences tested in the experiment.

This package contains two utilities:

## ASSOCIATION:
This utility takes in library association sequencing data (FASTQ) and a design file (FASTA) to assign barcodes to the corresponding elements tested. Functionality includes filtering for quality and coverage of barcodes. This utility must be run before the COUNT utility.

## COUNT:
This utility processes sequence data (FASTQ) of barcodes from the DNA and RNA fractions of the MPRA experiment and outputs count tables labeled with the element tested and a label provided in the design file. This utility can process multiple replicates and conditions in a parallelized manner. Based on a user specified flag, the pipeline will either output normalized activity for each tested sequence, or will combine the results into a single count matrix compatible with MPRAnalyze.


## Installation

### Required packages

- conda

Download here: `https://docs.conda.io/en/latest/miniconda.html`



### Clone repository

```bash
git clone https://github.com/shendurelab/MPRAflow.git
```

### Set up conda environment:
This pipeline uses python2.7 and python3.6 and is set up to run on a Linux system. Two .yml files are provided to create the appropriate environments. The general environment with nextflow located in the home directory called `environment.yml` and a specific python 2.7 environment in the `conf` folder: `mpraflow_py27.yml`.

The different environments are handled internally by nextflow. Therefore your compute node, where you start MPRAflow, have to have access to the internet.

Install the the conda environment. The general conda environment is called `MPRAflow`.
```bash
cd MPRAflow
conda env create -n MPRAflow -f environment.yml
```

If you do not have access to the internet, you have to run the previous command on a node with internet. Afterwards you need to start nextflow too (see `Steps to run the pipeline`). After creation of the second conda environment by nextflow you can cancel it and start it on your internal node. Be aware that folders must have access on all nodes.

Nextflow has problems using conda 4.7 and highet, because the `source activate` command is replaced by `conda activate`. If you get error messages after running you can make a symbolik link of the `activate` command from you `bin` folder of the `conda` or `miniconda` folder to your `MPRAflow` environment `bin` folder. E.g. like:

```bash
ln -s ~/miniconda3/bin/activate ~/miniconda3/envs/MPRAflow/bin/activate
```



## Running the pipeline

#### Steps to run the pipeline

This pipeline comes with a `nextflow.config` file set up to run on HPC clusters, allowing each process to be run as a separate 'qsub' command. The config contains example code for SGE, LSF, and SLURM architectures. The default is SGE.
Please remove the `\\` for the architecture you would like to use and place `\\` in front of any architectures not currently in use. A '\\' in front of all of them runs the pipeline on your local machine. If you run MPRAflow on a cluster system make sure be that you export all environment variables. E.g. this can be done with the `-V` option by SGE.

**NOTE:**: Please consult your cluster's wiki page for cluster specific commands and change `clusterOptions = ` to reflect these specifications. Additionally, for large libraries, more memory can be specified in this location.

Please use a submit script for steps 2 and 3. For full details of mandatory and optional arguments run:

 ```bash
 conda activate MPRAflow
 nextflow run count.nf --help
 nextflow run association.nf --help
 ```

This pipeline expects the FASTQ files to be demultiplexed and trimmed to only include sequence from the insert, barcodes, and/or UMIs.

## Quick Start

1. Create an 'experiment' csv in the format below, including the header. `DNA_R1` or `RNA_R1` is name of the gziped fastq of the forward read of the DNA or RNA from the defined condition and replicate. `DNA_R2` or `RNA_R2` is the corresponding index read with UMIs (excluding sample barcodes) and `DNA_R3` or `RNA_R3` of the reverse read. If you do not have UMIs remove the columns `DNA_R2` and `RNA_R2` or leave them empty.

   ```
   Condition,Replicate,DNA_R1,DNA_R2,DNA_R3,RNA_R1,RNA_R2,RNA_R3
   condition1,1,cond1_rep1_DNA_FWD_reads.fastq.gz,cond1_rep1_DNA_IDX_reads.fastq.gz,cond1_rep1_DNA_REV_reads.fastq.gz,cond1_rep1_RNA_FWD_reads.fastq.gz,cond1_rep1_RNA_IDX_reads.fastq.gz,cond1_rep1_RNA_REV_reads.fastq.gz
   condition1,2,cond1_rep2_DNA_FWD_reads.fastq.gz,cond1_rep2_DNA_IDX_reads.fastq.gz,cond1_rep2_DNA_REV_reads.fastq.gz,cond1_rep2_RNA_FWD_reads.fastq.gz,cond1_rep2_RNA_IDX_reads.fastq.gz,cond1_rep2_RNA_REV_reads.fastq.gz
   condition2,1,cond2_rep1_DNA_FWD_reads.fastq.gz,cond2_rep1_DNA_IDX_reads.fastq.gz,cond2_rep1_DNA_REV_reads.fastq.gz,cond2_rep1_RNA_FWD_reads.fastq.gz,cond2_rep1_RNA_IDX_reads.fastq.gz,cond2_rep1_RNA_REV_reads.fastq.gz
   condition2,2,cond2_rep2_DNA_FWD_reads.fastq.gz,cond2_rep2_DNA_IDX_reads.fastq.gz,cond2_rep2_DNA_REV_reads.fastq.gz,cond2_rep2_RNA_FWD_reads.fastq.gz,cond2_rep2_RNA_IDX_reads.fastq.gz,cond2_rep2_RNA_REV_reads.fastq.gz
   ```

2. If you would like each insert to be colored based on different user-specified categories, such as "positive control", "negative control", "shuffled control", and "putative enhancer", to assess the overall quality the user can create a 'label' tsv in the format below that maps the name to category:

   ```
   insert1_name insert1_label
   insert2_name insert2_label
   ```
   The insert names must exactly match the names in the design FASTA file.

3. Run Association if using a design with randomly paired candidate sequences and barcodes

   ```bash
   conda activate MPRAflow
   nextflow run association.nf --fastq-insert "${fastq_prefix}_R1_001.fastq.gz" --design "ordered_candidate_sequences.fa" --fastq-bc "${fastq_prefix}_R2_001.fastq.gz"
   ```
    **NOTE:** This will run in local mode, please submit this command to your cluster's queue if you would like to run a parallelized version.

4. Run Count

   ```bash
   conda activate MPRAflow
   nextflow run count.nf --dir "bulk_FASTQ_directory" --e "experiment.csv" --design "ordered_candidate_sequences.fa" --association "dictionary_of_candidate_sequences_to_barcodes.p"
   ```
   Be sure that the `experiment.csv` is correct. All fastq files must be in the same folder given by the `--dir` option. If you do not have UMIs please use the option `--no-umi`.



5. Run Saturation mutagenesis

    ```bash
    conda activate MPRAflow
    nextflow run saturationMutagenesis.nf --dir "directory_of_DNA/RNA_counts" --e "satMutexperiment.csv" --assignment "yourSpecificAssignmentFile.variants.txt.gz"
    ```

    **Note** The experiment file is different from the count workflow. It should contain the condition, replicate and filename of the counts, like:

    ```
    Condition,Replicate,COUNTS
    condition1,1,cond1_1_counts.tsv.gz
    condition1,2,cond1_2_counts.tsv.gz
    condition1,3,cond1_3_counts.tsv.gz
    condition2,1,cond2_1_counts.tsv.gz
    condition2,2,cond2_2_counts.tsv.gz
    condition2,3,cond2_3_counts.tsv.gz
    ```

    The count files can be generated by the count workflow, are named: `<condition>_<replicate>_counts.tsv.gz` and can be found in the `outs/<condition>/<replicate>` folder. THey have to be copied or linked into the `--dir` folder.


## Example files can be found in the example folder in this repository

![MPRA_nextflow](https://github.com/shendurelab/MPRAflow/blob/master/MPRA_nextflow.png)
