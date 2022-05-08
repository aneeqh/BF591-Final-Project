#import libraries
library(tidyverse)
library(shiny)
library(ggplot2)
library(DT)

ui <- fluidPage(
  titlePanel("Counts Matrix Exploration"),
  p("On this page you can upload the normalized counts matrix and explore the data."),
  sidebarLayout(
    sidebarPanel(
      fileInput(inputId = 'countsFP', label = 'Load normalized counts matrix CSV'),
      sliderInput(inputId = 'percvar', label = 'Select the minimum percentile variance of genes', min = 0, max = 100, value = 80),
      sliderInput(inputId = 'nonzero', label = 'Select the minimum number of non-zero samples', min = 0, max = 69, value = 5),
      submitButton(text='Submit !')
    ),
    mainPanel(
      tabsetPanel(
        tabPanel('Summary', tableOutput(outputId = 'normsumm')),
        tabPanel('Diagnostic Plots', p('Please wait 10-15 seconds after submitting for the plots to load'), plotOutput(outputId = 'medvar'), plotOutput(outputId = 'medzero')),
        tabPanel('Heatmap', p('Please wait 20 seconds after submitting for the heatmap to load'), plotOutput(outputId = "hmap")),
        tabPanel('PCA', selectInput(inputId = "comp1", label="Select X-axis", choices = c("PC1", "PC2", "PC3", "PC4", "PC5", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")),
          selectInput(inputId = "comp2", label="Select Y-axis", choices = c("PC1", "PC2", "PC3", "PC4", "PC5", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10"), selected = "PC2"),
          plotOutput(outputId = "PCAplot"))
      ))
  ))
server <- function(input, output, session){
  options(shiny.maxRequestSize=30*1024^2)
  #function to take norm counts input file
  load_counts <- reactive({
    if (!is.null(input$countsFP)){
      meta <- read_csv(input$countsFP$datapath)
      return(meta)}
    else{
      return(NULL)
    }
  })
  #function to produce summary table
  counts_summ <- function(counts_tib, perc_var, nz_genes){
    if (!is.null(input$countsFP)){
      tot_samp <- ncol(counts_tib)-1  #store original number of samples and genes
      tot_genes <- nrow(counts_tib)
      counts_tib <- counts_tib %>% mutate(variance = apply(counts_tib[-1], MARGIN = 1, FUN = var))  #calculate variance and then percentile
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
      summ_tib <- tibble('Measure' = c('Number of Samples', 'Number of Genes', 'No. of genes passing filters', "Percentage (%) passing filter", 'No. of genes not passing filters', 'Percentafe (%) not passing filter'),
                         'Value' = c(tot_samp, tot_genes, filt_genes, perc_pass_genes, fail_genes, perc_fail_genes))
      return(summ_tib)}
    else{return(NULL)}
  }
  #function to produce median vs variance plot
  med_vs_var <- function(counts_tib, perc_var){
    if (!is.null(input$countsFP)){
      #make a plot tibble
      plot_tib <- counts_tib%>%
        mutate(Median = apply(counts_tib[-1], MARGIN = 1, FUN = median), 
               Variance = apply(counts_tib[-1], MARGIN = 1, FUN = var))
      perc_val <- quantile(plot_tib$Variance, probs = perc_var/100)   #calculate percentile
      plot_tib <- plot_tib %>% mutate(thresh = case_when(Variance >= perc_val ~ "TRUE", TRUE ~ "FALSE"))
      #plot scatter plot
      cols <- c("FALSE" = "red", "TRUE" = "black")
      scatter <- ggplot(plot_tib, aes(Median, Variance))+
        geom_point(aes(color=thresh), alpha=0.75)+
        scale_color_manual(values = cols)+
        labs(title = 'Plot of Median vs Variance.', subtitle = "Genes filtered out are in red.")+
        scale_y_log10()+
        scale_x_log10()+
        theme_bw()+
        theme(legend.position = 'bottom')
      return(scatter)}
      else{return(NULL)}
  }
  #function to produce median vs non-zero samples plot
  med_vs_nz <- function(counts_tib, nz_genes){
    if (!is.null(input$countsFP)){
      tot_samp <- ncol(counts_tib)-1  #store original number of samples and genes
      #make a plot tibble
      plot_tib <- counts_tib %>%   
        mutate(Median = apply(counts_tib[-1], MARGIN = 1, FUN = median)) %>% na_if(0)  #calc median, convert 0 to NA
      plot_tib$no_zeros <- rowSums(is.na(plot_tib))  #make new col, with counts.
      plot_tib <- plot_tib %>% mutate(thresh = case_when(no_zeros <= nz_genes ~ "TRUE", TRUE ~ "FALSE"))
      #plot scatter plot
      cols <- c("FALSE" = "red", "TRUE" = "black")
      scatter <- ggplot(plot_tib, aes(Median, no_zeros))+
        geom_point(aes(color=thresh), alpha=0.75)+
        scale_color_manual(values = cols)+
        scale_x_log10()+
        labs(title = 'Plot of Median vs Number of Non-Zero genes', subtitle = "Genes filtered out are in red.")+
        theme_bw()+
        ylab('Number of samples with zero count')+
        theme(legend.position = 'bottom')
      return(scatter)}
    else{return(NULL)}
  }
  #function to produce heatmap
  plot_heatmap <- function(counts_tib, perc_var){
    if (!is.null(input$countsFP)){
      counts_tib <- log10(counts_tib[-1])
      #produce plot_tib
      plot_tib <- counts_tib %>% 
        mutate(variance = apply(counts_tib, MARGIN = 1, FUN = var))
      perc_val <- quantile(plot_tib$variance, probs = perc_var/100, na.rm = TRUE)   #calculate percentile
      plot_tib <- filter(plot_tib, variance >= perc_val) #filter the tibble
      hmap <- heatmap(as.matrix(plot_tib[-ncol(plot_tib)]), scale = "row")
      return(hmap)}
    else{return(NULL)}
  }
  #function to produce PCA plot
  plot_pca <- function(counts_tib, perc_var, comp1, comp2){
    if (!is.null(input$countsFP)){
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
      plot_tib <- tibble(PC1 = pca_res$x[,comp1], PC2=pca_res$x[,comp2])
      pca <- ggplot(plot_tib, aes(PC1, PC2))+
        geom_point()+
        labs(title="Princple Component Analysis Plot")+
        xlab(str_c(comp1, x, "% variance", sep=" "))+
        ylab(str_c(comp2, y, "% variance", sep=" "))+
        theme_bw()
      return(pca)}
    else{return(NULL)}
  }
  #methods to render tables and plots
  output$normsumm <- renderTable({
    counts_summ(load_counts(), input$percvar, input$nonzero)
  })
  output$medvar <- renderPlot({
    med_vs_var(load_counts(), input$percvar)
  })
  output$medzero <- renderPlot({
    med_vs_nz(load_counts(), input$nonzero)
  })
  output$hmap <- renderPlot({
    plot_heatmap(load_counts(), input$percvar)
  })
  output$PCAplot <- renderPlot({
    plot_pca(load_counts(), input$percvar, input$comp1, input$comp2)
  })
}

shinyApp(ui=ui, server = server)