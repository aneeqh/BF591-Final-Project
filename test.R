#import libraries
library(tidyverse)
library(shiny)
library(ggplot2)

#set filepaths
meta_path = 'data/sample_metadata.csv'
norm_counts_path = 'data/norm_counts.csv'
deseq_res_path = 'data/GSE64810_mlhd_DESeq2_diffexp_DESeq2_outlier_trimmed_adjust.txt'

#import data
meta = read_csv(meta_path)

#function to process metadata table into summary table-
summ_table <- function(meta_tib){
  summ_tib <- tibble("Columns" = colnames(meta_tib), "Type" = sapply(meta_tib, class),
                    "Mean" = sapply(meta_tib, mean, na.rm = TRUE), "SD" = sapply(meta_tib, sd, na.rm = TRUE))
  return(summ_tib)
}

plot_aod <- function(meta_tib){
  histo <- ggplot(meta_tib, aes(Age_of_death))+
    geom_histogram(bins = 10)+
    labs(title = 'Histogram of Age of Death')+
    theme_bw()
  return(histo)
}

#normalize the counts data
counts <- read_csv(norm_counts_path)
