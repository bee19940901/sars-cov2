library(tidyverse)
library(argparse)

get.args <- function(){
  pa <- ArgumentParser()
  pa$add_argument("-if", "--in_file", type="character", required=TRUE, help="测序深度文件")
  pa$add_argument("-of", "--out_pdf", type="character", required=TRUE, help="pdf结果文件")
  pa$add_argument("-og", "--out_png", type="character", required=TRUE, help="结果文件")
  args <- pa$parse_args()
  if(!is.null(args$help) || is.null(args$in_file) || is.null(args$out_pdf) || is.null(args$out_png)){
    pa$print_usage()
    quit(status = 1)
  } else{
    return(args)
  }
}

main <- function(){
  args <- get.args()
  p <- 
    read_tsv(args$in_file, col_names = FALSE) %>%
    ggplot(mapping = aes(x=X2, y=log2(X3+1))) +
    geom_area(fill="pink", linewidth=0.1) +
    ggtitle("Depth Distribution") +
    xlab("Position") +
    ylab("Log2(Depth)") +
    theme_classic()
  ggsave(args$out_pdf)
  ggsave(args$out_png)
}

main()
quit(status = 0)



