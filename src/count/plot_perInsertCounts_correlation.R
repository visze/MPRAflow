
#adapted from Vikram Agarwal by Gracie Gordon

cpath <- grep('conda', .libPaths(), value=TRUE, ignore.case=TRUE)
.libPaths(cpath)


library(ggplot2)
library(optparse)
library(cowplot)
library(dplyr)


option_list <- list(
    make_option(c("-c", "--condition"), type="character",
        help="Condition name"),
    make_option(c("-l", "--label"), type="character",
        help="Label file. (optional)"),
    make_option(c("-f", "--files"), type="character",
        help="Comma separated input files of assigned counts"),
    make_option(c("-r", "--replicates"), type="character",
        help="Comma separated name of the replicates (same order than files)"),
    make_option(c("-t", "--threshold"), type="integer", default=10,
        help="Number of required barcodes (default 10)")
)

parser <- OptionParser(option_list=option_list)
arguments <- parse_args(parser, positional_arguments=TRUE)
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

# condition
cond=opt$condition
# labels
if ("label" %in% names(opt)) {
  label_f=as.data.frame(read.table(opt$label, sep='\t',header=T,stringsAsFactors = F))
  colnames(label_f)=c('name','label')
  useLabels=TRUE
} else {
  useLabels=FALSE
}


# replicates and count files
files <- strsplit(opt$files,",")[[1]]
replicates=strsplit(opt$replicates,",")[[1]]
if (length(files) != length(replicates)) {
    stop("Number of input files must be euqal to number of replicates")
}

num_replicates=length(replicates)

data <- data.frame(File=files,Replicate=replicates)
data$Condition <- cond

print(data)



# pairwise comparison only if more than one replicate
thresh=opt$threshold

plot_correlations <- function(data,condition,r1,r2,name){
  dna_p <- ggplot(data, aes(log2(dna_normalized.x), log2(dna_normalized.y))) +
              geom_point(aes(colour = label), show.legend = TRUE) +
              xlim(-5,5) + ylim(-5,5) +
              xlab(sprintf(paste("log2 Normalized DNA count per insert,\n replicate", r1))) +
              ylab(sprintf(paste("log2 Normalized DNA count per insert,\n replicate", r2))) +
              geom_text(x=0, y=4.5,label=sprintf("   r = %.2f", cor(log2(data$dna_normalized.x),log2(data$dna_normalized.y),method="pearson")),size=10) +
              geom_text(x=0, y=4, label=sprintf("rho = %.2f", cor(data$dna_normalized.x,data$dna_normalized.y,method="spearman")),size=10) +
              geom_abline(intercept = 0, slope = 1) +
              theme_classic(base_size = 30)
  ggsave(sprintf("%s_%s_%s_DNA_%s.png",condition,as.character(r1),as.character(r2),name), plot=dna_p, width=15,height=15)

  rna_p <- ggplot(data, aes(log2(rna_normalized.x), log2(rna_normalized.y))) +
              geom_point(aes(colour = label), show.legend = TRUE) +
              xlim(-5,5) + ylim(-5,5) +
              xlab(sprintf(paste("log2 Normalized RNA count per insert,\n replicate", r1))) +
              ylab(sprintf(paste("log2 Normalized RNA count per insert,\n replicate", r2))) +
              geom_text(x=0, y=4.5,label=sprintf("   r = %.2f", cor(log2(data$rna_normalized.x),log2(data$rna_normalized.y),method="pearson")),size=10) +
              geom_text(x=0, y=4, label=sprintf("rho = %.2f", cor(data$rna_normalized.x,data$rna_normalized.y,method="spearman")),size=10) +
              geom_abline(intercept = 0, slope = 1) +
              theme_classic(base_size = 30)
  ggsave(sprintf("%s_%s_%s_RNA_%s.png",condition,as.character(r1),as.character(r2),name), plot=rna_p,width=15,height=15)

  ratio_p <- ggplot(data, aes(log2(ratio.x), log2(ratio.y))) +
                geom_point(aes(colour = label), show.legend = TRUE) +
                xlim(-5,5) + ylim(-5,5) +
                xlab(sprintf(paste("log2 RNA/DNA per insert,\n replicate", r1))) +
                ylab(sprintf(paste("log2 RNA/DNA per insert,\n replicate", r2))) +
                geom_text(x=0, y=4.5,label=sprintf("   r = %.2f", cor(log2(data$ratio.x),log2(res$ratio.y),method="pearson")),size=10) +
                geom_text(x=0, y=4, label=sprintf("rho = %.2f", cor(data$ratio.x,data$ratio.y,method="spearman")),size=10) +
                geom_abline(intercept = 0, slope = 1) +
                theme_classic(base_size = 30)
  ggsave(sprintf("%s_%s_%s_Ratio_%s.png",condition,as.character(r1),as.character(r2),name), plot=ratio_p,width=15,height=15)
}

write_correlation <- function(data,condition,r1,r2,name){
  outs = as.character(sprintf("%s vs %s:  RNA: %.5f DNA: %.5f ratio: %.5f NormSymmetry: %d",r1,r2,
      cor(data$rna_normalized.x,data$rna_normalized.y,method="spearman"),
      cor(data$dna_normalized.x,data$dna_normalized.y,method="spearman"),
      cor(data$ratio.x,data$ratio.y,method="spearman"),
      abs(length(which((data$ratio.x-data$ratio.y)>0))-length(which((data$ratio.x-data$ratio.y)<0)))
      + abs(length(which((data$ratio.x-data$ratio.y)>0))-length(which((data$ratio.x-data$ratio.y)<0)))
      + abs(length(which((data$ratio.x-data$ratio.y)>0))-length(which((data$ratio.x-data$ratio.y)<0)))
  ))

  write(outs,file=sprintf("%s_%s_%s_%s.txt",condition,as.character(r1),as.character(r2),name),append=TRUE)
}

if(data %>% nrow >1){

  # make pairwise combinations
  selected <- combn(data$Replicate,2)
  print('sel')
  print(selected)


  print('reps')
  for(i in seq(1,dim(selected)[2])){
    print(selected[,i])
    r1=selected[1,i]
    r2=selected[2,i]
    data1 <- read.table(as.character((data %>% filter(Replicate == r1))$File),as.is=T,sep="\t",header=T,stringsAsFactors = F)
    data2 <- read.table(as.character((data %>% filter(Replicate == r2))$File),as.is=T,sep="\t",header=T,stringsAsFactors = F)


    # Remove unassigned barcodes
    data1<-data1 %>% filter(name != 'no_BC')
    data2<-data2 %>% filter(name != 'no_BC')

    res <- data1 %>% inner_join(data2,by=c('name'))
    if (useLabels){
      res <- res %>% inner_join(label_f, by=c('name'))
    } else {
      res$label = 'NA'
    }

    plot_correlations(res,cond,r1,r2,"pairwise_minThreshold")
    write_correlation(res,cond,r1,r2,"correlation_minThreshold")
    res <- res %>% filter(n_obs_bc.x >= thresh, n_obs_bc.y >= thresh)
    plot_correlations(res,cond,r1,r2,"pairwise_minThreshold")
    write_correlation(res,cond,r1,r2,"correlation_minThreshold")
  }
}


print('hist')

all=data.frame()
plots = list()

for(n in 1:(data%>%nrow)){
    assigned_counts <- read.table(as.character(data[n,]$File),as.is=T,sep="\t",header=T,stringsAsFactors = F) %>% filter(name != 'no_BC')
    intercept <- median(assigned_counts$n_obs_bc)
    plots[[n]] <- ggplot(assigned_counts, aes(x=n_obs_bc)) +
      geom_histogram(bins=300) +
      geom_vline(xintercept=intercept, colour="red") +
      xlim(0,300) + ggtitle(paste("replicate", data[n,]$Replicate, sep=' '))
    all <- bind_rows(assigned_counts)

}

hist_plot <- do.call("plot_grid", c(plots))

ggsave(sprintf("%s_barcodesPerInsert.png",cond),hist_plot)

#
print('boxplot')

 if (useLabels){
   all <- all %>% inner_join(label_f, by=c('name'))
 } else {
   all$label = 'NA'
 }

print('merged')
print(head(all))
all$name <- factor(all$name)
all$log2 <- as.numeric(as.character(all$log2))
all$label <- as.factor(all$label)

all=all[-1,]
all=all[order(all$log2),]

bymedian=with(all,reorder(name,-log2,median,order=TRUE))
all$name=factor(all$name,levels=levels(bymedian))

all_subsample <- all %>% group_by(label) %>% sample_n(min(1000,n)) %>% ungroup()

bp <- ggplot(all_subsample, aes(x=name, y=log2, color=label)) +
      geom_boxplot() +
      xlab('insert') +
      ylab('log2 fold change') +
      theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank(), axis.line=element_line(colour="black"), axis.title.x=element_text(size=15), axis.title.y=element_text(size=15), axis.text.y=element_text(size=15), legend.text=element_text(size=15))

ggsave(sprintf("%s_all_barcodesPerInsert_box.png",cond),bp)

bp <- ggplot(all, aes(x=label, y=log2, fill=label)) +
      geom_violin() +
      geom_boxplot(width=0.1,fill='white') +
      xlab('insert') +
      ylab('log2 fold change') +
      theme(axis.text.x=element_text(angle=90,hjust=1,size=15), panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank(), axis.line=element_line(colour="black"), axis.title.x=element_text(size=15), axis.title.y=element_text(size=15), axis.text.y=element_text(size=15), legend.text=element_text(size=15))

ggsave(sprintf("%s_group_barcodesPerInsert_box.png",cond),bp)
