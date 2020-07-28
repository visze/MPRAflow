#Adapted from Vikram Agarwal by Gracie Gordon

cpath <- grep('conda', .libPaths(), value=TRUE, ignore.case=TRUE)
.libPaths(cpath)

library(dplyr)
library(optparse)


option_list <- list(
    make_option(c("-c", "--condition"), type="character",
        help="Condition name"),
    make_option(c("-l", "--label"), type="character",
        help="Label file. (optional)"),
    make_option(c("-f", "--files"), type="character",
        help="Comma separated input files of assigned counts"),
    make_option(c("-r", "--replicates"), type="character",
        help="Comma separated name of the replicates (same order than files)"),
    make_option(c("-o", "--output"), type="character",
        help="Output file of master table"),
    make_option(c("-s", "--statistic"), type="character",
        help="Statistic of master table"),
    make_option(c("-t", "--threshold"), type="integer", default=10,
        help="Number of required barcodes (default 10)")
)

arguments <- parse_args(OptionParser(option_list=option_list), positional_arguments=TRUE)
opt <- arguments$options

if (!"condition" %in% names(opt)) {
  stop("--condition parameter must be provided. See script usage (--help)")
}
if (!"files" %in% names(opt)) {
  stop("--files parameter must be provided. See script usage (--help)")
}
if (!"replicates" %in% names(opt)) {
  stop("--replicates parameter must be provided. See script usage (--help)")
}
if (!"output" %in% names(opt)) {
  stop("--output parameter must be provided. See script usage (--help)")
}
if (!"statistic" %in% names(opt)) {
  stop("--statistic parameter must be provided. See script usage (--help)")
}



#exp=args[1]
cond=opt$cond
thresh=opt$threshold
files <- strsplit(opt$files,",")[[1]]
replicates=strsplit(opt$replicates,",")[[1]]
if (length(files) != length(replicates)) {
    stop("Number of input files must be euqal to number of replicates")
}
outfile=opt$output
avg_outfile=opt$statistic
#out=args[3]

nrm_reps=c()
all_reps=c()
##MAKE MASTER TABLE

for (i in length(files)){
   file=files[i]
   rep=replicates[i]

   tab=as.data.frame(read.table(file,header=TRUE))

   filter_tab=tab[tab$n_obs_bc >= thresh,]

   n_inserts=(dim(filter_tab)[1])

   cond_col=as.data.frame(rep(cond,n_inserts))
   rep_col=as.data.frame(rep(rep,n_inserts))
   colnames(cond_col)='condition'
   colnames(rep_col)='replicate'

   pref=as.data.frame(cbind(cond_col,rep_col))
   lab_tab=as.data.frame(cbind(pref,filter_tab))

   ##NORMALISZE PER REPLCATE TO COMBINE IN NEXT FUNCTION
   nrm_reps1=lab_tab
   nrm_reps1$ratio=(nrm_reps1$ratio)/(median(nrm_reps1$ratio))
   nrm_reps1$log2=round(log2(nrm_reps1$ratio),8)
   nrm_reps=rbind(nrm_reps,nrm_reps1)
   all_reps=rbind(all_reps,lab_tab)
}

write.table(all_reps,file=outfile,quote=FALSE,sep='\t',row.names=FALSE,col.names=TRUE )


##MAKE AVERAGED ACROSS REPLICATES

all_avg <- all_reps %>% group_by(condition) %>% summarize(
                    mean_ratio=mean(ratio),
                    mean_log2=log2(mean(ratio)),
                    mean_n_obs_bc=mean(n_obs_bc)
                  )


write.table(all_avg,file=avg_outfile,quote=FALSE,sep='\t',row.names=FALSE,col.names=TRUE )

print('done')
