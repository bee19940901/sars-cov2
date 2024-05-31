library(argparse)
library(tidyverse)

get.args <- function(){
  parser <- ArgumentParser(description="Coverage Distribution Volin.")
  parser$add_argument("-if", "--in_file", type="character", required=TRUE)
  parser$add_argument("-of", "--out_pdf", type="character", required=TRUE)
  parser$add_argument("-og", "--out_png", type="character", required=TRUE)
  args <- parser$parse_args()
  if(!is.null(args$help) || is.null(args$in_file) || is.null(args$out_pdf) || is.null(args$out_png)){
    parser$parse_args()
    quit(status = 1)
  } else {
    return(args)
  }
}

main <- function(){
  
  args <- get.args()
  df <- read_tsv(args$in_file)
  
  mgi <- df %>%
    filter(platform=="MGI") %>%
    select(qc.overallScore)
  
  np <- df %>%
    filter(platform=="NanoPore") %>%
    select(qc.overallScore)
  
  p.value <- 
    wilcox.test(
      mgi$qc.overallScore,
      np$qc.overallScore,
      alternative='two.sided'
    )$p.value
  
  if(p.value > 0.05){
    sub.title <- sprintf("P-value=ns")
  } else if (p.value <= 0.05 ) {
    sub.title <- sprintf("P-value = %s", p.value)
  }
  
  p <-
    df %>%
    ggplot(mapping = aes(x = platform, y=log2(qc.overallScore + 1), fill=platform)) +
    geom_violin() +
    geom_boxplot(width=0.1, position=position_dodge(0.9), outlier.shape = NA) +
    geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.9), shape = 16, size = 1.5, alpha=.5) +
    ggtitle("Consensus Genomes QC-value Distribution", subtitle = sub.title) +
    xlab("Platform") +
    ylab("Log2(QC-value+1)") +
    theme_classic() +
    theme(legend.position = "none")
  ggsave(args$out_png)
  ggsave(args$out_pdf)
}


main()
quit(status = 0)


