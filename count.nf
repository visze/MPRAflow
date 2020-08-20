#!/usr/bin/env nextflow
/*
========================================================================================
                         MPRAflow
========================================================================================
MPRA Analysis Pipeline.
Started 2019-07-29.
Count Utility
#### Homepage / Documentation
https://github.com/shendurelab/MPRAflow
#### Authors
Vikram Agarwal <thefinitemachine@gmail.com>
Gracie Gordon <gracie.gordon@ucsf.edu>
Martin Kircher <martin.kircher@bihealth.de>
Max Schubach <max.schubach@bihealth.de>

----------------------------------------------------------------------------------------
*/

def helpMessage() {
    log.info"""
    =========================================
    shendurelab/MPRAflow v${params.version}
    =========================================
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run main.nf
    Mandatory arguments:
      --dir                         Fasta directory (must be surrounded with quotes)
      --e, --experiment-file        Experiment csv file

    Recommended (Needed for the full count workflow. Can be neglected when only counts for saturation mutagenesis are needed):
      --design                      Fasta of ordered insert sequences.
      --association                 Pickle dictionary from library association process.

    Options:
      --labels                      tsv with the oligo pool fasta and a group label (ex: positive_control), a single label will be applied if a file is not specified
      --outdir                      The output directory where the results will be saved (default outs)
      --bc-length                   Barcode length (default 15)
      --umi-length                  UMI length (default 10)
      --no-umi                      Use this flag if no UMI is present in the experiment (default with UMI)
      --merge-intersect             Only retain barcodes in RNA and DNA fraction, and not 0 counts; inner join instead of full join (TRUE/FALSE, default: FALSE)
      --mpranalyze                  Only generate MPRAnalyze outputs
      --thresh                      minimum number of observed barcodes to retain insert (default 10)

    Extras:
      --h, --help                   Print this help message
      --name                        Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.
    """.stripIndent()
}

/*
* SET UP CONFIGURATION VARIABLES
*/

// Show help message
if (params.containsKey('h') || params.containsKey('help')){
    helpMessage()
    exit 0
}


// Configurable variables
params.name = false
params.email = false
params.plaintext_email = false
output_docs = file("$baseDir/docs/output.md")

//defaults
results_path = params.outdir

params.mpranalyze=false
params.thresh=10

// Validate Inputs

// experiment file saved in params.experiment_file
if (params.containsKey('e')){
    params.experiment_file=file(params.e)
} else if (params.containsKey("experiment-file")) {
    params.experiment_file=file(params["experiment-file"])
} else {
    exit 1, "Experiment file not specified with --e or --experiment-file"
}
if( !params.experiment_file.exists()) exit 1, "Experiment file ${params.experiment_file} does not exist"

// design file saved in params.design_file
if ( params.containsKey("design")){
    params.design_file=file(params.design)
    if( !params.design_file.exists() ) exit 1, "Design file ${params.design} does not exist"
    if( !params.containsKey("association") ) exit 1, "Association file has to be specified with --association when using option --design"
} else if (params.mpranalyze) {
    exit 1, "Design file has to be specified with --association when using flag --mpranalyze"
} else {
    log.info("Running MPRAflow count without design file.")
}

// Association file in params.association_file
if (params.containsKey("association")){
    params.association_file=file(params.association)
    if( !params.association_file.exists() ) exit 1, "Association pickle ${params.association_file} does not exist"
} else if (params.mpranalyze) {
    exit 1, "Association file has to be specified with --association when using flag --mpranalyze"
} else {
    log.info("Running MPRAflow count without association file.")
}

// label file saved in label_file
if (params.containsKey("labels")){
    label_file=file(params.labels)
    if (!label_file.exists()) {
      println "Label file ${label_file} does not exist. Use NA as label!"
      label_file=file("NA")
    }
} else {
    label_file=file("NA")
}

// merge intersect (inner join)
if (params.containsKey("merge-intersect")){
    params.merge_intersect = true
} else {
    params.merge_intersect = false
}

// check if UMIs or not are present
if (params.containsKey("no-umi")){
    params.no_umi = true
} else {
    params.no_umi = false
}

// BC length
if (params.containsKey("bc-length")){
    params.bc_length = params["bc-length"]
} else {
    params.bc_length = 15
}
// UMI length
if (params.containsKey("umi-length")){
    params.umi_length = params["umi-length"]
} else {
    params.umi_length = 10
}

// Create FASTQ channels
if (params.no_umi) {
  reads_noUMI = Channel.fromPath(params.experiment_file).splitCsv(header: true).flatMap{
    row -> [
      tuple(row.Condition,row.Replicate,"DNA",
        [row.Condition,row.Replicate,"DNA"].join("_"),
        file([params.dir,"/",row.dna,row.DNA_BC_F].join()),
        file([params.dir,"/",row.dna,row.DNA_BC_R].join())
      ),
      tuple(row.Condition,row.Replicate,"RNA",
        [row.Condition,row.Replicate,"RNA"].join("_"),
        file([params.dir,"/",row.rna,row.RNA_BC_F].join()),
        file([params.dir,"/",row.rna,row.RNA_BC_R].join())
      )
    ]
  }
} else {
  reads = Channel.fromPath(params.experiment_file).splitCsv(header: true).flatMap{
    row -> [
      tuple(row.Condition, row.Replicate, "DNA",
        [row.Condition,row.Replicate,"DNA"].join("_"),
        file([params.dir,"/",row.DNA_BC_F].join()),
        file([params.dir,"/",row.DNA_UMI].join()),
        file([params.dir,"/",row.DNA_BC_R].join()),
      ),
      tuple(row.Condition, row.Replicate, "RNA",
        [row.Condition,row.Replicate,"RNA"].join("_"),
        file([params.dir,"/",row.RNA_BC_F].join()),
        file([params.dir,"/",row.RNA_UMI].join()),
        file([params.dir,"/",row.RNA_BC_R].join()),
      ),
    ]
  }
}


// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
custom_runName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
  custom_runName = workflow.runName
}



// Header log info
log.info """=======================================================
                                          ,--./,-.
          ___     __   __   __   ___     /,-._.--~\'
    |\\ | |__  __ /  ` /  \\ |__) |__         }  {
    | \\| |       \\__, \\__/ |  \\ |___     \\`-._,-`-,
                                          `._,._,\'
MPRAflow v${params.version}"
======================================================="""
def summary = [:]
summary['Pipeline Name']    = 'shendurelab/MPRAflow'
summary['Pipeline Version'] = params.version
summary['Run Name']         = custom_runName ?: workflow.runName

//summary['Thread fqdump']    = params.threadfqdump ? 'YES' : 'NO'
summary['Output dir']       = params.outdir
summary['Working dir']      = workflow.workDir
//summary['Container Engine'] = workflow.containerEngine
//if(workflow.containerEngine) summary['Container'] = workflow.container
summary['Current home']     = "$HOME"
summary['Current user']     = "$USER"
summary['Current path']     = "$PWD"
summary['Working dir']      = workflow.workDir
summary['Output dir']       = params.outdir
summary['Script dir']       = workflow.projectDir
summary['Config profile']   = workflow.profile
summary['Experiment file']  = params.experiment_file
summary['Reads']            = (params.no_umi ? reads_noUMI : reads.print())
summary['Design file']      = if (params.design) params.design_file
if (params.association) summary['Association file'] = params.association_file
summary['UMIs']             = (params.no_umi ? "Reads without UMI" : "Reads with UMI")
if (!params.no_umi) summary['UMI length']       = params.umi_length 
summary['BC length']        = params.bc_length
summary['BC threshold']     = params.thresh
summary['Non zero counts (merge-intersect)']     = params.merge_intersect
summary['mprAnalyze']       = params.mpranalyze

if(params.email) summary['E-mail Address'] = params.email
log.info summary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "========================================="

// Check that Nextflow version is up to date enough
// try / throw / catch works for NF versions < 0.25 when this was implemented
try {
    if( ! nextflow.version.matches(">= $params.nf_required_version") ){
        throw GroovyException('Nextflow version too old')
    }
} catch (all) {
    log.error "====================================================\n" +
              "  Nextflow version $params.nf_required_version required! You are running v$workflow.nextflow.version.\n" +
              "  Pipeline execution will continue, but things may break.\n" +
              "  Please run `nextflow self-update` to update Nextflow.\n" +
              "============================================================"
}

println 'start analysis'

/*
* STEP 1: Create BAM files
* contributions: Martin Kircher, Max Schubach, & Gracie Gordon
*/
//if UMI
if (!params.no_umi) {
    process 'create_BAM' {
        tag "make idx"
        label 'longtime'

        conda 'conf/mpraflow_py27.yml'

        input:
            tuple val(cond), val(rep), val(type), val(datasetID), file(fw_fastq), file(umi_fastq), file(rev_fastq) from reads
            val(bc_length) from params.bc_length
        output:
            tuple val(cond), val(rep), val(type),val(datasetID),file("${datasetID}.bam") into clean_bam
        shell:
            """
            #!/bin/bash
            echo $datasetID

            echo $fw_fastq
            echo $umi_fastq
            echo $rev_fastq

            umi_length=`zcat $umi_fastq | head -2 | tail -1 | wc -c`
            umi_length=\$(expr \$((\$umi_length-1)))

            fwd_length=`zcat $fw_fastq | head -2 | tail -1 | wc -c`
            fwd_length=\$(expr \$((\$fwd_length-1)))

            rev_start=\$(expr \$((\$fwd_length+1)))

            rev_length=`zcat $rev_fastq | head -2 | tail -1 | wc -c`
            rev_length=\$(expr \$((\$rev_length-1)))

            minoverlap=`echo \${fwd_length} \${fwd_length} $bc_length | awk '{print (\$1+\$2-\$3-1 < 11) ? \$1+\$2-\$3-1 : 11}'`

            echo \$rev_start
            echo \$umi_length
            echo \$minoverlap

            paste <( zcat $fw_fastq ) <( zcat $rev_fastq  ) <( zcat $umi_fastq ) | awk '{if (NR % 4 == 2 || NR % 4 == 0) {print \$1\$2\$3} else {print \$1}}' | python ${"$baseDir"}/src/count/FastQ2doubleIndexBAM.py -p -s \$rev_start -l 0 -m \$umi_length --RG ${datasetID} | python ${"$baseDir"}/src/MergeTrimReadsBAM.py --FirstReadChimeraFilter '' --adapterFirstRead '' --adapterSecondRead '' -p --mergeoverlap --minoverlap \$minoverlap > ${datasetID}.bam
            """
    }
}

//if no UMI
/*
* contributions: Martin Kircher, Max Schubach, & Gracie Gordon
*/
if (params.no_umi) {
    process 'create_BAM_noUMI' {
        tag "make idx"
        label 'longtime'

        conda 'conf/mpraflow_py27.yml'

        input:
            tuple val(cond), val(rep),val(type),val(datasetID),file(fw_fastq), file(rev_fastq) from reads_noUMI
            val(bc_length) from params.bc_length
        output:
            tuple val(cond), val(rep),val(type),val(datasetID),file("${datasetID}.bam") into clean_bam
        when:
            params.no_umi
        shell:
            """
            #!/bin/bash
            echo $datasetID

            echo $fw_fastq
            echo $rev_fastq

            fwd_length=`zcat $fw_fastq | head -2 | tail -1 | wc -c`
            fwd_length=\$(expr \$((\$fwd_length-1)))

            rev_start=\$(expr \$((\$fwd_length+1)))

            rev_length=`zcat $rev_fastq | head -2 | tail -1 | wc -c`
            rev_length=\$(expr \$((\$rev_length-1)))

            minoverlap=`echo \${fwd_length} \${fwd_length} $bc_length | awk '{print (\$1+\$2-\$3-1 < 11) ? \$1+\$2-\$3-1 : 11}'`

            echo \$rev_start
            echo \$minoverlap

            paste <( zcat $fw_fastq ) <(zcat $rev_fastq  ) | \
            awk '{
                if (NR % 4 == 2 || NR % 4 == 0) {
                  print \$1\$2
                } else {
                  print \$1
                }}' | python ${"$baseDir"}/src/count/FastQ2doubleIndexBAM.py -p -s \$rev_start -l 0 -m 0 --RG ${datasetID} | python ${"$baseDir"}/src/MergeTrimReadsBAM.py --FirstReadChimeraFilter '' --adapterFirstRead '' --adapterSecondRead '' -p --mergeoverlap  --minoverlap \$minoverlap> ${datasetID}.bam
              """
    }
}


/*
* STEP 2: create raw counts
* contributions: Martin Kircher, Max Schubach, & Gracie Gordon
*/

process 'raw_counts'{
    label 'shorttime'

    conda 'conf/mpraflow_py36.yml'

    publishDir "$params.outdir/$cond/$rep"

    input:
        tuple val(cond), val(rep), val(type), val(datasetID), file(bam) from clean_bam
        val(umi_length) from params.umi_length
    output:
        tuple val(cond), val(rep), val(type), val(datasetID), val("raw"), file("${datasetID}_raw_counts.tsv.gz") into raw_ct, raw_ct_stats
    script:
        if(params.no_umi)
            """
            #!/bin/bash

            samtools view -F 1 -r $datasetID $bam | \
            awk '{print \$10}' | \
            sort | \
            gzip -c > ${datasetID}_raw_counts.tsv.gz
            """

        else
            """
            #!/bin/bash

            samtools view -F 1 -r $datasetID $bam | \
            awk -v 'OFS=\\t' '{ for (i=12; i<=NF; i++) {
              if (\$i ~ /^XJ:Z:/) print \$10,substr(\$i,6,$umi_length)
            }}' | \
            sort | uniq -c | \
            awk -v 'OFS=\\t' '{ print \$2,\$3,\$1 }' | \
            gzip -c > ${datasetID}_raw_counts.tsv.gz
            """

}

/*
* STEP 3: Filter counts for correct barcode length
* contributions: Martin Kircher, Max Schubach, & Gracie Gordon
*/

process 'filter_counts'{
    label 'shorttime'
    publishDir "$params.outdir/$cond/$rep"

    conda 'conf/mpraflow_py27.yml'

    input:
        tuple val(cond), val(rep), val(type), val(datasetID), val(rawType), file(rc) from raw_ct
        val(bcLength) from params.bc_length
    output:
        tuple val(cond), val(rep), val(type), val(datasetID), val("filtered"), file("${datasetID}_filtered_counts.tsv.gz") into filter_ct, filter_ct_stats
    shell:
        """
        bc=$bcLength
        echo \$bc
        zcat $rc | grep -v "N" | \
        awk -v var="\$bc" -v 'OFS=\t' '{ if (length(\$1) == var) { print } }' | \
        sort | \
        gzip -c > ${datasetID}_filtered_counts.tsv.gz
        """
}

/*
* Statistic
*/
raw_ct_stats.concat(filter_ct_stats).into{ ct_stats_1; ct_stats_2 }

process 'statistic_counts' {
      label 'shorttime'

      input:
          tuple val(cond), val(rep), val(type), val(datasetID), val(countType),file(rc) from ct_stats_1
      output:
          tuple val(cond), val(rep), val(type), val(datasetID), val(countType), file("${datasetID}_${countType}_count.stats") into count_stats
      shell:
          """
          paste <( echo "$cond") <( echo "$rep") <( echo "$type") \
          <( zcat $rc | \
            awk -v OFS='\t' 'BEGIN{pbar="NA"}{ count += \$NF; lines+=1; if (pbar != \$1) { barcodes+=1 }; pbar=\$1 }END{ print count,lines,barcodes }' ) \
          <( zcat $rc | cut -f 2 | sort -u | wc -l ) \
          > ${datasetID}_${countType}_count.stats
          """
}

process 'count_stats_merge'{
    label 'shorttime'

    result = count_stats.groupTuple(by: 4).multiMap{i ->
                              countType: i[4]
                              files: i[5]
                            }

    input:
        file(countStatFiles) from result.files
        val(countType) from result.countType
    output:
        tuple val(countType),file("${countType}_count_stats.tsv") into count_stats_merge
    script:
        countStat = countStatFiles.collect{"$it"}.join(' ')
    shell:
        """
        cat $countStat | sort -k1,1 -k3,3 -k2,2 > ${countType}_count_stats.tsv
        """
}

process 'statistic_BC_in_RNA_DNA' {
      label 'shorttime'

      input:
          tuple val(cond),val(rep),val(typeA),val(typeB),val(datasetIDA),val(datasetIDB),val(typeCount),file(countA),file(countB) from ct_stats_2.groupTuple(by: [0,1,4]).map{i -> i.flatten()}
      output:
          tuple val(cond), val(rep), val(typeCount), file("${cond}_${rep}_${typeCount}_BC_in_RNA_DNA.stats") into BC_in_RNA_DNA_stats
      script:
          def dna = typeA == 'DNA' ? countA : countB
          def rna = typeA == 'DNA' ? countB : countA

          """
          paste <( echo "$cond") <( echo "$rep") \
            <( join <( zcat $dna | cut -f 1 | sort | uniq ) \
            <( zcat $rna | cut -f 1 | sort | uniq ) | wc -l ) \
            > ${cond}_${rep}_${typeCount}_BC_in_RNA_DNA.stats
          """
}

process 'stats_BC_in_RNA_DNA_merge'{
    label 'shorttime'

    result = BC_in_RNA_DNA_stats.groupTuple(by: 2).multiMap{i ->
                              countType: i[2]
                              files: i[3]
                            }

    input:
        file(countStatFiles) from result.files
        val(countType) from result.countType
    output:
        tuple val(countType),file("${countType}_BC_in_RNA_DNA.tsv") into stats_BC_in_RNA_DNA_merge
    script:
        countStat = countStatFiles.collect{"$it"}.join(' ')
    shell:
        """
        cat $countStat | sort -k1,1 -k2,2 > ${countType}_BC_in_RNA_DNA.tsv
        """
}


process 'stats_final'{
    label 'shorttime'
    publishDir "$params.outdir", mode:'copy'

    conda 'conf/mpraflow_r.yml'

    result = count_stats_merge.concat(stats_BC_in_RNA_DNA_merge).groupTuple(by: 0).multiMap{i ->
                              countType: i[0]
                              files: i[1]
                            }

    input:
        file(statFiles) from result.files
        val(countType) from result.countType
    output:
        tuple val(countType),file("statistic_${countType}_count.tsv") into count_stats_final
    script:
        def stat = statFiles.toList()
        def count = stat[0]
        def shared = stat[1]
        """
        Rscript ${"$baseDir"}/src/count/combine_count_stats.R --count $count --shared $shared --output statistic_${countType}_count.tsv
        """
}

/*
* STEP 4: Record overrepresended UMIs and final count table
* contributions: Martin Kircher, Max Schubach, & Gracie Gordon
*/
if (params.no_umi) {
  process 'final_counts_no_UMI'{
      label 'shorttime'
      publishDir "$params.outdir/$cond/$rep", mode:'copy'

      input:
          tuple val(cond), val(rep),val(type),val(datasetID),val(filteredType),file(fc) from filter_ct
      output:
          tuple val(cond), val(rep),val(type),val(datasetID),file("${datasetID}_counts.tsv.gz") into final_count
      script:
          """
          #!/bin/bash

          zcat $fc | awk '{print \$1}' | \
          uniq -c | \
          gzip -c > ${datasetID}_counts.tsv.gz
          """
  }
} else {
  process 'final_counts'{
      label 'shorttime'
      publishDir "$params.outdir/$cond/$rep", mode:'copy'

      input:
          tuple val(cond), val(rep),val(type),val(datasetID),val(filteredType),file(fc) from filter_ct
      output:
          file("${datasetID}_freqUMIs.txt") into frequent_UMI
          tuple val(cond), val(rep),val(type),val(datasetID),file("${datasetID}_counts.tsv.gz") into final_count
      script:
          """
          #!/bin/bash

          for i in $fc; do
            echo \$(basename \$i);
            zcat \$i | cut -f 2 | sort | uniq -c | sort -nr | head;
            echo;
          done > ${datasetID}_freqUMIs.txt

          zcat $fc | awk '{print \$1}' | \
          uniq -c | \
          gzip -c > ${datasetID}_counts.tsv.gz
          """
  }
}

/*
* STEP 5: Merge each DNA and RNA file
* contributions: Gracie Gordon, Max Schubach
*/
process 'dna_rna_merge_counts'{
    publishDir "$params.outdir/$cond/$rep", mode:'copy'
    label 'shorttime'

    conda 'conf/mpraflow_py36.yml'

    input:
        tuple val(cond),val(rep),val(typeA),val(typeB),val(datasetIDA),val(datasetIDB),file(countA),file(countB) from final_count.groupTuple(by: [0,1]).map{i -> i.flatten()}
    output:
        tuple val(cond), val(rep), file("${cond}_${rep}_counts.tsv.gz") into merged_dna_rna, merged_dna_rna2
    script:
        def dna = typeA == 'DNA' ? countA : countB
        def rna = typeA == 'DNA' ? countB : countA
        if (params.merge_intersect) {
            """
            join -1 1 -2 1 -t"\$(echo -e '\\t')" \
            <( zcat  $dna | awk 'BEGIN{ OFS="\\t" }{ print \$2,\$1 }' | sort ) \
            <( zcat $rna | awk 'BEGIN{ OFS="\\t" }{ print \$2, \$1 }' | sort) | \
            gzip -c > ${cond}_${rep}_counts.tsv.gz
            """
        } else {
            """
            join -e 0 -a1 -a2 -t"\$(echo -e '\\t')" -o 0 1.2 2.2 \
            <( zcat  $dna | awk 'BEGIN{ OFS="\\t" }{ print \$2,\$1 }' | sort ) \
            <( zcat $rna | awk 'BEGIN{ OFS="\\t" }{ print \$2, \$1 }' | sort) | \
            gzip -c > ${cond}_${rep}_counts.tsv.gz
            """
        }
            
}


process 'correlate_BC_counts'{
    publishDir "$params.outdir/$cond", mode:'copy'
    label 'longtime'

    conda 'conf/mpraflow_r.yml'

    result = merged_dna_rna2.groupTuple(by: 0, sort: true).multiMap{i ->
                              cond: i[0]
                              replicate: i[1]
                              files: i[2]
                            }
    input:
        file(pairlistFiles) from result.files
        val(replicate) from result.replicate
        val(cond) from result.cond
    output:
        file "${cond}_barcode_DNA_pairwise.png"
        file "${cond}_barcode_RNA_pairwise.png"
        file "${cond}_barcode_Ratio_pairwise.png"
        file "${cond}_barcode_correlation.tsv" into bc_correlation
        file "${cond}_DNA_perBarcode.png"
        file "${cond}_RNA_perBarcode.png"
    script:
        collected=pairlistFiles.collect{"$it"}
        pairlist = []
        for (i = 0; i < collected.size(); i++) {
          pairlist << "${cond}_${replicate[i]}_counts.tsv.gz"
        }
        pairlist = pairlist.join(',')
        replicates= replicate.join(',')
        """
        Rscript ${"$baseDir"}/src/count/plot_perBCCounts_correlation.R \
        --condition $cond \
        --files $pairlist --replicates $replicates
        """
}

process 'combine_bc_correlation' {
    label 'shorttime'
    publishDir "$params.outdir", mode:'copy'

    conda 'conf/mpraflow_py36.yml'

    input:
        file(bc_correlation_collected) from bc_correlation.collect()
    output:
        file("statistic_bc_correlation.tsv")
    script:
        collected = bc_correlation_collected.collect{"$it"}.join(' ')
    shell:
        """
        (
          cat $collected | head -n 1;
          for i in $collected; do
            cat \$i | tail -n +2;
          done;
        ) > statistic_bc_correlation.tsv
        """
}

//MPRAnalyze option
//contributions: Gracie Gordon
if(params.mpranalyze){
    /*
    * STEP 6: Merge all DNA/RNA counts into one big file
    * contributions: Gracie Gordon, Max Schubach
    */

    process 'final_merge'{
        label 'longtime'
        publishDir "$params.outdir/$cond", mode:'copy'

        conda 'conf/mpraflow_py36.yml'

        result = merged_dna_rna.groupTuple(by: 0).multiMap{i ->
                                  cond: i[0]
                                  replicates: i[1]
                                  files: i[2]
                                }


        input:
            file(pairlistFiles) from result.files
            val(replicate) from result.replicates
            val(cond) from result.cond
        output:
            tuple val(cond),file("${cond}_count.csv.gz") into merged_out
        script:
            collected=pairlistFiles.collect{"$it"}
            res=[]
            for (i = 0; i < collected.size(); i++) {
               r=replicate[i]
               f="${cond}_${r}_counts.tsv.gz"
               res.add("--counts $r $f")
            }
            counts=res.join(' ')
        shell:
            """
            python ${"$baseDir"}/src/count/merge_all.py --condition $cond --output ${cond}_count.csv.gz $counts
            """
    }


    /*
    * STEP 7: Add label to outfile
    * contributions: Gracie Gordon, Max Schubach
    */

    process 'final_label'{
        label 'shorttime'
        publishDir "$params.outdir/$cond", mode:'copy'

        conda 'conf/mpraflow_py36.yml'

        input:
            tuple val(cond), file(table) from merged_out
            file(des) from params.design_file
            file(association) from params.association_file
        output:
            tuple val(cond), file("${cond}_final_labeled_counts.tsv.gz") into labeled_out
        shell:
            """
            python ${"$baseDir"}/src/count/label_final_count_mat.py --counts $table --assignment $association --output ${cond}_final_labeled_counts.tsv.gz --design $des
            """
    }

    /*
    * STEP 8: Generate inputs
    * contributions: Tal Ashuach, Max Schubach
    */

    process 'generate_mpranalyze_inputs'{
        label 'shorttime'
        publishDir "$params.outdir/$cond", mode:'copy'

        conda 'conf/mpraflow_py36.yml'

        input:
            tuple val(cond),file(labeled_file) from labeled_out
        output:
            file("rna_counts.tsv.gz") into mpranalyze_rna_counts
            file("dna_counts.tsv.gz") into mpranalyze_dna_counts
            file("rna_annot.tsv.gz") into mpranalyze_rna_annotation
            file("dna_annot.tsv.gz") into mpranalyze_fna_annotation
        shell:
            """
            python ${"$baseDir"}/src/count/mpranalyze_compiler.py --input $labeled_file --rna-counts-output rna_counts.tsv.gz --dna-counts-output dna_counts.tsv.gz --rna-annotation-output rna_annot.tsv.gz --dna-annotation-output dna_annot.tsv.gz
            """
    }

}


/*
* STEP 5: Merge each DNA and RNA file label with sequence and insert and normalize
* contributions: Gracie Gordon Max Schubach
*/
//merge and normalize
if(!params.mpranalyze && params.containsKey("association")){

    process 'dna_rna_merge'{
        label 'longtime'
        publishDir "$params.outdir/$cond/$rep", mode:'copy', pattern: '*_*_assigned_counts.tsv.gz'

        conda 'conf/mpraflow_py36.yml'

        input:
            tuple val(cond), val(rep), file(counts) from merged_dna_rna
            file(des) from params.design_file
            file(association) from params.association_file
        output:
             tuple val(cond), val(rep), file("${cond}_${rep}_assigned_counts.tsv.gz") into merged_ch, merged_ch2
             tuple val(cond), val(rep), file("${cond}_${rep}_assigned_counts.statistic.tsv.gz") into assigned_stats
        script:
            if (params.merge_intersect) {
                filter = "--minRNACounts 1 --minDNACounts 1"
            } else {
                filter = "--minRNACounts 0 --minDNACounts 0"
            }
        shell:
            """
            python ${"$baseDir"}/src/count/merge_label.py --counts ${counts} \
            ${filter} \
            --assignment $association --design $des \
            --output ${cond}_${rep}_assigned_counts.tsv.gz \
            --statistic ${cond}_${rep}_assigned_counts.statistic.tsv.gz
            """

    }

    // collect statistics per condition
    process 'combine_stats_dna_rna_merge' {
        label 'shorttime'
        publishDir "$params.outdir/$cond", mode:'copy'

        conda 'conf/mpraflow_py36.yml'

        result = assigned_stats.groupTuple(by: 0, sort: true).multiMap{i ->
                                  cond: i[0]
                                  replicate: i[1]
                                  files: i[2]
                                }

          input:
              file(pairlistFiles) from result.files
              val(replicate) from result.replicate
              val(cond) from result.cond
        output:
             file("${cond}_assigned_counts.statistic.tsv.gz") into statistic_assignment_per_cond
       script:
           collected = pairlistFiles.collect{"$it"}
           res=[]
           for (i = 0; i < collected.size(); i++) {
              r=replicate[i]
              f="${cond}_${r}_assigned_counts.statistic.tsv.gz"
              res.add("--statistic $r $f")
           }
           statistic=res.join(' ')
        shell:
            """
            python ${"$baseDir"}/src/count/merge_statistic_tables.py \
            --condition $cond \
            $statistic \
            --output ${cond}_assigned_counts.statistic.tsv.gz
            """
    }

    // combine statistics of all conditions
    process 'combine_stats_dna_rna_merge_all' {
        label 'shorttime'
        publishDir "$params.outdir", mode:'copy'

        conda 'conf/mpraflow_py36.yml'

        input:
            file(statistic_per_cond) from statistic_assignment_per_cond.collect()
        output:
            file("statistic_assigned_counts.tsv")
        script:
            collected = statistic_per_cond.collect{"$it"}.join(' ')
        shell:
            """
            (zcat $collected | head -n 1;
              for i in $collected; do
                zcat \$i | tail -n +2
              done) > statistic_assigned_counts.tsv
            """
    }

    /*
    * STEP 6: Calculate correlations between Replicates
    * contributions: Vikram Agarwal, Gracie Gordon, Max Schubach
    */
    process 'calc_correlations'{
        label 'shorttime'
        publishDir "$params.outdir/$cond", mode:'copy'

        conda 'conf/mpraflow_r.yml'

        result = merged_ch.groupTuple(by: 0, sort: true).multiMap{i ->
                                  cond: i[0]
                                  replicate: i[1]
                                  files: i[2]
                                }

        input:
            file(pairlistFiles) from result.files
            val(replicate) from result.replicate
            val(cond) from result.cond
            file(lab) from label_file
        output:
            file "${cond}_all_barcodesPerInsert_box.png"
            file "${cond}_all_barcodesPerInsert_box_minThreshold.png"
            file "${cond}_DNA_pairwise.png"
            file "${cond}_DNA_pairwise_minThreshold.png"
            file "${cond}_group_barcodesPerInsert_box.png"
            file "${cond}_group_barcodesPerInsert_box_minThreshold.png"
            file "${cond}_Ratio_pairwise.png"
            file "${cond}_Ratio_pairwise_minThreshold.png"
            file "${cond}_RNA_pairwise.png"
            file "${cond}_RNA_pairwise_minThreshold.png"
            file "${cond}_correlation.tsv" into oligo_correlation
            file "${cond}_correlation_minThreshold.tsv" into oligo_correlation_thresh
        script:
            collected=pairlistFiles.collect{"$it"}
            pairlist = []
            for (i = 0; i < collected.size(); i++) {
              pairlist << "${cond}_${replicate[i]}_assigned_counts.tsv.gz"
            }
            pairlist = pairlist.join(',')
            replicates= replicate.join(',')
            def label = lab.exists() ? "--label $lab" : ""
            """
            Rscript ${"$baseDir"}/src/count/plot_perInsertCounts_correlation.R \
            --condition $cond $label \
            --files $pairlist \
            --replicates $replicates \
            --threshold $params.thresh
            """
    }
    /*
    * combine oligo correlations
    */
    process 'combine_oligo_correlation' {
        label 'shorttime'
        publishDir "$params.outdir", mode:'copy'

        conda 'conf/mpraflow_py36.yml'

        input:
            file(oligo_correlation) from oligo_correlation.collect()
            file(oligo_correlation_thresh) from oligo_correlation_thresh.collect()
        output:
            file("oligo_correlation.tsv")
        script:
            oligo_correlation_collected = oligo_correlation.collect{"$it"}.join(' ')
            oligo_correlation_thresh_collected = oligo_correlation_thresh.collect{"$it"}.join(' ')
        shell:
            """
            (
              cat $oligo_correlation_collected | head -n 1 | awk -v 'OFS=\\t' '{print \$0,"threshold (min $params.thresh)"}';
              for i in $oligo_correlation_collected; do
                cat \$i | tail -n +2 | awk -v 'OFS=\\t' '{print \$0,"False"}'
              done;
              for i in $oligo_correlation_thresh_collected; do
                cat \$i | tail -n +2 | awk -v 'OFS=\\t' '{print \$0,"True"}'
              done;
            ) > oligo_correlation.tsv
            """
    }



    /*
    * contributions: Vikram Agarwal & Gracie Gordon
    */
    process 'make_master_tables' {
        label 'shorttime'
        publishDir "$params.outdir/$cond", mode:'copy'

        conda 'conf/mpraflow_r.yml'

        result = merged_ch2.groupTuple(by: 0, sort: true).multiMap{i ->
                                  cond: i[0]
                                  replicate: i[1]
                                  files: i[2]
                                }

        input:
            file(pairlistFiles) from result.files
            val(replicate) from result.replicate
            val(cond) from result.cond
        output:
            file "average_allreps.tsv.gz"
            file "allreps.tsv.gz"
            file "allreps_minThreshold.tsv.gz"
        script:
            collected=pairlistFiles.collect{"$it"}
            pairlist = []
            for (i = 0; i < collected.size(); i++) {
              pairlist << "${cond}_${replicate[i]}_assigned_counts.tsv.gz"
            }
            pairlist = pairlist.join(',')
            replicates= replicate.join(',')
        shell:
            """
            Rscript ${"$baseDir"}/src/count/make_master_tables.R --condition $cond \
            --threshold $params.thresh \
            --files $pairlist --replicates $replicates \
            --output allreps_minThreshold.tsv.gz \
            --output-all allreps.tsv.gz \
            --statistic average_allreps.tsv.gz
            """
    }

}
