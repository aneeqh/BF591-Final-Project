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
      sliderInput(inputId = 'percvar', label = 'Select the minimum percentile variance of genes', min = 0, max = 100, value = 10),
      sliderInput(inputId = 'nonzero', label = 'Select the number of non-zero samples', min = 0, max = 69, value = 5),
      submitButton(text='Submit !')
    ),
    mainPanel(
      tabsetPanel(
        tabPanel('Summary', tableOutput(outputId = 'normsumm')),
        tabPanel('Diagnostic Plots'),
        tabPanel('Heatmap'),
        tabPanel('PCA')
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
      print('Here!')
      tot_samp <- ncol(counts_tib)-1  #store original number of samples and genes
      tot_genes <- nrow(counts_tib)
      counts_tib <- counts %>% mutate(variance = apply(counts[-1], MARGIN = 1, FUN = var))  #calculate variance and then percentile
      perc_val <- quantile(counts_tib$variance, probs = perc_var/100)   #calculate percentile
      counts_tib <- filter(counts_tib, variance >= perc_val)  #filter by percentile
      print('here2')
      counts_tib <- na_if(counts_tib, 0)    #make zeroes NA's
      counts_tib$non_zero <- tot_samp-rowSums(is.na(counts_tib))  #calculate no. of genes with enough non zero samples
      counts_tib <- filter(counts_tib, non_zero >= nz_genes)  #filter by non-zero samples
      print('here3')
      filt_genes <- nrow(counts_tib)    #calculate the number and % of genes passing the filters
      perc_pass_genes <- filt_genes/tot_genes*100
      fail_genes <- tot_genes-filt_genes
      perc_fail_genes <- fail_genes/tot_genes*100
      #produce the summary tibble
      print('here4')
      summ_tib <- tibble('Measure' = c('Number of Samples', 'Number of Rows', 'No. of genes passing filters', "Percentage (%) passing filter", 'No. of genes not passing filters', 'Percentafe (%) not passing filter'),
                         'Value' = c(tot_samp, tot_genes, filt_genes, perc_pass_genes, fail_genes, perc_fail_genes))
      return(summ_tib)}
    else{return(NULL)}
  }
  
  #methods to render tables and plots
  output$normsumm <- renderTable({
    counts_summ(load_counts(), input$percvar, input$nonzero)
  })
}

shinyApp(ui=ui, server = server)