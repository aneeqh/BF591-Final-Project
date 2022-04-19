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

counts_summ <- function(counts_tib, perc_var, ngenes){
  #original number of samples and genes-
  nsamp <- ncol(counts_tib)-1
  ngenes <- nrows(counts_tib)
  #calculate variance and then percentile-
  counts_tib <- counts_tib %>% rowwise() %>%
    mutate(variance = var(c_across(2:69)))
  #calculate percentile -
  perc_val <- quantile(counts_tib$variance, probs = perc_var/100)
  #filter by percentile-
  counts_tib <- filter(counts_tib, variance >= perc_val)
  #calculate no. of genes with enough non zero samples
  
  #then produce the summary tibble
  summ_tib <- tibble('Measure' = c('Number of Samples', 'Number of Rows') , 'Value' = c(nsamp, ngenes))
    
  return(summ_tib)
}
counts_summ(counts)