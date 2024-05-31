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
  mgi.CV <- df %>%
    filter(platform=="MGI") %>%
    select(CV)
  np.CV <- df %>%
    filter(platform=="MGI") %>%
    select(CV)
  if(shapiro.test(mgi.CV$CV)$p.value >0.05 && shapiro.test(np.CV$CV)$p.value > 0.05){
    p.value <- t.test(mgi.CV$CV, np.CV$CV, alternative = 'two.sided')$p.value
  } else {
    p.value <- wilcox.test(mgi.CV$CV, np.CV$CV, alternative = 'two.sided')$p.value
  }
  if(p.value > 0.05){
    main.title <- sprintf("CV Distribution, P-value > 0.05")
  } else {
    main.title <- sprintf("CV Distribution, P-value = %s", p.value)
  }
  p <- df %>%
    ggplot(mapping = aes(x = platform, y=CV, fill=platform)) +
    geom_violin() +
    geom_boxplot(width=0.1, position=position_dodge(0.9), outlier.shape = NA) +
    geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.9), shape = 16, size = 1.5, alpha=.5) +
    ggtitle(main.title) +
    xlab("Platform") +
    ylab("CV") +
    theme_classic() +
    theme(legend.position = "none")
  ggsave(args$out_png)
  ggsave(args$out_pdf)
}

main()
quit(status = 0)


