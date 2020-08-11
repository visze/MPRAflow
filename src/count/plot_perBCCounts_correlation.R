
cpath <- grep('conda', .libPaths(), value=TRUE, ignore.case=TRUE)
.libPaths(cpath)

library(optparse)
library(cowplot)
library(tidyverse)
#library(viridis)


option_list <- list(
    make_option(c("-c", "--condition"), type="character",
        help="Condition name"),
    make_option(c("-f", "--files"), type="character",
        help="Comma separated input files of assigned counts"),
    make_option(c("-r", "--replicates"), type="character",
        help="Comma separated name of the replicates (same order than files)")
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

plot_correlations_dna <- function(data,condition,r1,r2,name){
  max <- max(log2(data$`DNA.y`))
  dna_p <- ggplot(data, aes(log2(DNA.x), log2(DNA.y))) +
              stat_density_2d(aes(fill = stat(level)), geom = 'polygon') +
              scale_fill_viridis_c(name = "density") +
              geom_point() +
              xlab(sprintf(paste("log2 Normalized DNA count per barcode,\n replicate", r1))) +
              ylab(sprintf(paste("log2 Normalized DNA count per barcode,\n replicate", r2))) +
              geom_text(x=0, y=max-0.5,label=sprintf("   r = %.2f", cor(log2(data$DNA.x),log2(data$DNA.y),method="pearson")),size=10) +
              geom_text(x=0, y=max-1.0, label=sprintf("rho = %.2f", cor(data$DNA.x,data$DNA.y,method="spearman")),size=10) +
              geom_abline(intercept = 0, slope = 1) +
              theme_classic(base_size = 30)
  return(dna_p)
}
plot_correlations_rna <- function(data,condition,r1,r2,name){
  max <- max(log2(data$`RNA.y`))
  rna_p <- ggplot(data, aes(log2(RNA.x), log2(RNA.y))) +
              stat_density_2d(aes(fill = stat(level)), geom = 'polygon') +
              scale_fill_viridis_c(name = "density") +
              geom_point() +
              xlab(sprintf(paste("log2 Normalized RNA count per barcode,\n replicate", r1))) +
              ylab(sprintf(paste("log2 Normalized RNA count per barcode,\n replicate", r2))) +
              geom_text(x=0, y=max-0.5,label=sprintf("   r = %.2f", cor(log2(data$RNA.x),log2(data$RNA.y),method="pearson")),size=10) +
              geom_text(x=0, y=max-1.0, label=sprintf("rho = %.2f", cor(data$RNA.x,data$RNA.y,method="spearman")),size=10) +
              geom_abline(intercept = 0, slope = 1) +
              theme_classic(base_size = 30)
  return(rna_p)
}
plot_correlations_ratio <- function(data,condition,r1,r2,name){
  max <- max(log2(data$`Ratio.y`))
  ratio_p <- ggplot(data, aes(log2(Ratio.x), log2(Ratio.y))) +
                stat_density_2d(aes(fill = stat(level)), geom = 'polygon') +
                scale_fill_viridis_c(name = "density") +
                geom_point() +
                xlab(sprintf(paste("log2 RNA/DNA per barcode,\n replicate", r1))) +
                ylab(sprintf(paste("log2 RNA/DNA per barcode,\n replicate", r2))) +
                geom_text(x=0, y=max-0.5,label=sprintf("   r = %.2f", cor(log2(data$Ratio.x),log2(res$Ratio.y),method="pearson")),size=10) +
                geom_text(x=0, y=max-1.0, label=sprintf("rho = %.2f", cor(data$Ratio.x,data$Ratio.y,method="spearman")),size=10) +
                geom_abline(intercept = 0, slope = 1) +
                theme_classic(base_size = 30)
  return(ratio_p)
}

getCorrelationStats <- function(data,condition,r1,r2,name){

  norm <- abs(length(which((data$Ratio.x-data$Ratio.y)>0)) - length(which((data$Ratio.x-data$Ratio.y)<0)))
  + abs(length(which((data$Ratio.x-data$Ratio.y)>0))-length(which((data$Ratio.x-data$Ratio.y)<0)))
  + abs(length(which((data$Ratio.x-data$Ratio.y)>0))-length(which((data$Ratio.x-data$Ratio.y)<0)))
  outs <- data.frame(
                Comparison = sprintf("%s vs %s",r1,r2),
                DNA = sprintf("%.5f", cor(data$DNA.x,data$DNA.y,method="spearman")),
                RNA = sprintf("%.5f", cor(data$RNA.x,data$RNA.y,method="spearman")),
                Ratio = sprintf("%.5f", cor(data$Ratio.x,data$Ratio.y,method="spearman")),
                NormSymmetry=norm, stringsAsFactors=FALSE)
  return(outs)

}

writeCorrelationPlots <- function(plots, name){
  correlation_plots <- cowplot::plot_grid(plotlist = plots, ncol = 1)
  # correlation_plots <- do.call("grid.arrange", c(plots))

  ggsave(name,correlation_plots, width=15,height=10*length(plots))
}

writeCorrelation <- function(correlations, name){
    write.table(correlations,file=name,quote=FALSE,sep='\t',row.names=FALSE,col.names=TRUE )
}

readData <- function(file) {
  data <- read.table(file,as.is=T,sep="\t",header=F,stringsAsFactors = F)
  colnames(data) <- c("Barcode", "DNA", "RNA")
  data <- data %>% filter(DNA>0,RNA>0) %>%
          mutate(DNA=DNA/sum(DNA)*10**6,RNA=RNA/sum(RNA)*10**6, Ratio=log2(DNA/RNA))
  return(data)
}

if(data %>% nrow >1){

  # make pairwise combinations
  selected <- combn(data$Replicate,2)
  print('sel')
  print(selected)

  plots_correlations_rna = list()
  plots_correlations_dna = list()
  plots_correlations_ratio = list()
  stats_correlations = data.frame()
  print('reps')
  for(i in seq(1,dim(selected)[2])){
    print(selected[,i])
    r1=selected[1,i]
    r2=selected[2,i]
    data1 <- readData(as.character((data %>% filter(Replicate == r1))$File))
    data2 <- readData(as.character((data %>% filter(Replicate == r2))$File))


    res <- data1 %>% inner_join(data2,by=c('Barcode'))

    plots_correlations_dna[[i]] <- plot_correlations_dna(res,cond,r1,r2,"pairwise")
    plots_correlations_rna[[i]] <- plot_correlations_rna(res,cond,r1,r2,"pairwise")
    plots_correlations_ratio[[i]] <- plot_correlations_ratio(res,cond,r1,r2,"pairwise")

    stats_correlations <- stats_correlations %>% bind_rows(getCorrelationStats(res,cond,r1,r2,"correlation"))

  }

  writeCorrelationPlots(plots_correlations_dna, sprintf("%s_barcode_DNA_pairwise.png",cond))
  writeCorrelationPlots(plots_correlations_rna, sprintf("%s_barcode_RNA_pairwise.png",cond))
  writeCorrelationPlots(plots_correlations_ratio, sprintf("%s_barcode_Ratio_pairwise.png",cond))

  writeCorrelation(stats_correlations, sprintf("%s_barcode_correlation.tsv",cond))
}
