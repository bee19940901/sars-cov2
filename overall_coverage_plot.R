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
  df <- read_table(args$in_file)
  mgi.coverages <- df %>%
    filter(platform=="MGI") %>%
    select(coverage)
  np.coverages <- df %>%
    filter(platform=="MGI") %>%
    select(coverage)
  if(shapiro.test(mgi.coverages$coverage)$p.value >0.05 && shapiro.test(np.coverages$coverage)$p.value > 0.05){
    p.value <- t.test(mgi.coverages$coverage, np.coverages$coverage, alternative = 'two.sided')$p.value
  } else {
    p.value <- wilcox.test(mgi.coverages$coverage, np.coverages$coverage, alternative = 'two.sided')$p.value
  }
  if(p.value >= 0.05){
    main.title <- sprintf("Coverage Distribution, P-value > 0.05")
  } else {
    main.title <- sprintf("Coverage Distribution, P-value = %s", p.value)
  }
  p <- df %>%
    ggplot(mapping = aes(y=coverage, x=platform, fill=platform)) +
    geom_violin() +
    geom_boxplot(width=0.2,position=position_dodge(0.9)) +
    ggtitle( main.title) +
    xlab("Platform") +
    ylab("Coverage(%)") +
    theme_classic() +
    theme(legend.position = "none")
  print(args$out_pdf)
  ggsave(args$out_png)
  ggsave(args$out_pdf)
}

main()
quit(status = 0)
  










