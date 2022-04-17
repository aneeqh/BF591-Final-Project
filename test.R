#import libraries
library(tidyverse)
library(shiny)
library(ggplot2)

#set filepaths
meta_path = 'data/sample_metadata.csv'
norm_counts_path = 'data/GSE64810_mlhd_DESeq2_norm_counts_adjust.txt'
deseq_res_path = 'data/GSE64810_mlhd_DESeq2_diffexp_DESeq2_outlier_trimmed_adjust.txt'

#import data
meta = read_csv(meta_path)

#function to process metadata table into summary table-
summ_table <- function(meta_tib){
  summ_tib <- tibble("Columns" = colnames(meta_tib), "Type" = sapply(meta_tib, class),
                    "Mean" = sapply(meta_tib, mean, na.rm = TRUE), "SD" = sapply(meta_tib, sd, na.rm = TRUE))
  return(summ_tib)
}
