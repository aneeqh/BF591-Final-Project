#import libraries
library(tidyverse)
library(shiny)
library(ggplot2)

ah#set filepaths
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

med_vs_var <- function(counts_tib, perc_var){
  #make a plot tibble
  plot_tib <- counts_tib%>%
    mutate(Median = apply(counts_tib[-1], MARGIN = 1, FUN = median), 
           Variance = apply(counts_tib[-1], MARGIN = 1, FUN = var))
  perc_val <- quantile(plot_tib$Variance, probs = perc_var/100)   #calculate percentile
  plot_tib <- plot_tib %>% mutate(thresh = case_when(Variance >= perc_val ~ "TRUE", TRUE ~ "FALSE"))
  #plot scatter plot
  cols <- c("FALSE" = "red", "TRUE" = "black")
  scatter <- ggplot(plot_tib, aes(Median, Variance))+
    geom_point(aes(color=thresh))+
    scale_color_manual(values = cols)+
    labs(title = 'Plot of Median vs Variance.', subtitle = "Genes filtered out are in red.")+
    scale_y_log10()+
    scale_x_log10()+
    theme_bw()
  return(scatter)
}
med_vs_var(counts, 10)

df <- counts %>%   
  mutate(Median = apply(counts[-1], MARGIN = 1, FUN = median)) %>% na_if(0)  #calc median, convert 0 to NA
df$no_zeros <- rowSums(is.na(df))  #make new col, with counts.
df <- df %>% mutate(thresh = case_when(no_zeros <= 6  ~ "TRUE", TRUE ~ "FALSE"))


#plotting heatmap
plot_heatmap <- function(counts_tib, perc_var){
  #produce plot_tib
  plot_tib <- counts_tib %>% 
    mutate(variance = apply(counts_tib[-1], MARGIN = 1, FUN = var), .after = gene)
  perc_val <- quantile(plot_tib$variance, probs = perc_var/100, na.rm = TRUE)   #calculate percentile
  plot_tib <- filter(plot_tib, variance >= perc_val) #filter the tibble
  hmap <- pheatmap::pheatmap(as.matrix(plot_tib[-c(1,2)]), scale = "row")
  return(hmap)
}

#pca plot
plot_pca <- function(counts_tib, perc_var, comp1, comp2){
  #make plot tib-
  filt_tib <- counts_tib %>% 
    mutate(variance = apply(counts_tib[-1], MARGIN = 1, FUN = var), .after = gene)
  perc_val <- quantile(filt_tib$variance, probs = perc_var/100, na.rm = TRUE)   #calculate percentile
  filt_tib <- filter(filt_tib, variance >= perc_val) #filter the tibble
  pca_res <- prcomp(t(filt_tib[-c(1,2)]), scale = FALSE) #transpose the data and perform PCA
  #extract variance
  variance <- summary(pca_res)$importance[2,]
  x <- round(variance[comp1]*100, 2)
  y <- round(variance[comp2]*100, 2)
  #produce PCA plot
  plot_tib <- tibble(PC1 = tr_pca$x[,1], PC2=tr_pca$x[,2])
  pca <- ggplot(plot_tib, aes(PC1, PC2))+
    geom_point()+
    labs(title="Princple Component Analysis Plot")+
    xlab(str_c(comp1, x, "% variance", sep=" "))+
    ylab(str_c(comp2, y, "% variance", sep=" "))+
    theme_bw()
  return(pca)
}