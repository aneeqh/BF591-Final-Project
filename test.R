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
#plotting histograms function
plot_aod <- function(meta_tib){
  histo <- ggplot(meta_tib, aes(Age_of_death))+
    geom_histogram(bins = 10)+
    labs(title = 'Histogram of Age of Death')+
    theme_bw()
  return(histo)
}

#normalize the counts data
counts <- read_csv(norm_counts_path)

counts_summ <- function(counts_tib, perc_var, nz_genes){
  tot_samp <- ncol(counts_tib)-1  #store original number of samples and genes
  tot_genes <- nrow(counts_tib)
  counts_tib <- counts %>% mutate(variance = apply(counts[-1], MARGIN = 1, FUN = var))  #calculate variance and then percentile
  perc_val <- quantile(counts_tib$variance, probs = perc_var/100)   #calculate percentile
  counts_tib <- filter(counts_tib, variance >= perc_val)  #filter by percentile
  counts_tib <- na_if(counts_tib, 0)    #make zeroes NA's
  counts_tib$non_zero <- tot_samp-rowSums(is.na(counts_tib))  #calculate no. of genes with enough non zero samples
  counts_tib <- filter(counts_tib, non_zero >= nz_genes)  #filter by non-zero samples
  filt_genes <- nrow(counts_tib)    #calculate the number and % of genes passing the filters
  perc_pass_genes <- filt_genes/tot_genes*100
  fail_genes <- tot_genes-filt_genes
  perc_fail_genes <- fail_genes/tot_genes*100
  #produce the summary tibble
  summ_tib <- tibble('Measure' = c('Number of Samples', 'Number of Rows', 'No. of genes passing filters', "Percentage (%) passing filter", 'No. of genes not passing filters', 'Percentafe (%) not passing filter'),
                     'Value' = c(tot_samp, tot_genes, filt_genes, perc_pass_genes, fail_genes, perc_fail_genes))
  return(summ_tib)
}
counts_summ(counts, 10, 65)


