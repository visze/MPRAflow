# MPRAflow Changelog

## development

### count.nf

- moving count scripts frpm `src` to `src/count`
- removing `merge_counts.py` script. merging is done by the process `dna_rna_merge_counts` and it is a simple but very fast bash `join` command.
- complete refactor of `merge_label.py` to speed up (minutes instead of days!). Especially using grouping option together with aggregation instead of an (endless) loop gives performance boost. Removing unnecessary filter (BC length == 15). Is already done by the `filter_counts` process. lso adding possibility to filter for min. number of counts in DNA or/and RNA (basically to avoid not informative zero counts in RNA).
- using command line options for ALL count scripts to make them better readable and usabale in other context.
- refactor of `plot_perInsertCounts_correlation.R`. Now a min threshold can be used and files are saved with and without min threshold. Downsample instances to (default) `10000` (if larger) for boxplot picture `all_barcodesPerInsert_box.png`. Otherwise the script takes days on larger data. Bringing correlation plots together reduce number of files.
- adding new count statistic (like number of barcodes, barcodes shared between RNA/DNA, etc) in `statistic_raw_count.tsv` or `statistic_filtered_count.tsv`
- compressing output files to save disk space

## v2.2

No workflow changes. Only a few fixes and some restructuring of configs. Using nextflow version 20.01 now!

### global changes

* nextflow version 20.01 is needed because of multiMap() function
* introducing new config file `conf/global.config` with global variables like the min. required nextflow version and the actual MPRAflow version.
* moving cluster config to a sepArate file: `conf/cluster.config`. Try to adapt the times to the sort and longtime labels. Modify SLURM queue to SLUM not SGE options.
* improved documentation

### saturationMutagenesis.nf

* Bugfix of default out dir. It was not set to `params.outdir = "outs"` so it tries to create a folder `null`. Now in `params.outdir`
* removing default `params.version` and `params.nf_required_version`. Now in `conf/global.config`
* Catching cases when barcode/p-value filtering produces 0 variants
* Change update depricated fork method. Now works with nextflow 20.01

### count.nf

* removing default `params.version`, `params.nf_required_version` and `params.outdir`. Now in `conf/global.config`.

### association.nf

* removing default `params.version`, `params.nf_required_version` and `params.outdir`. Now in `conf/global.config`.


## v2.1

Initial MPRAflow version for publication.
