#import libraries
library(tidyverse)
library(shiny)
library(ggplot2)
library(DT)

ui <- fluidPage(
  titlePanel("Sample Information Exploration"),
  p("On this page you can upload the metadata and explore the sample information matrix."),
  sidebarLayout(
    sidebarPanel(
      fileInput(inputId = 'metaFP', label = 'Load sample information matrix CSV')
    ),
    mainPanel(
      tabsetPanel(
       tabPanel('Summary', tableOutput(outputId = 'summ')),
       tabPanel('Metadata', dataTableOutput(outputId = "metadata")),
       tabPanel('Plots', plotOutput(outputId = 'aod'), plotOutput(outputId = 'rin'), plotOutput(outputId = 'pmi'))
      ))
  )
)

server <- function(input, output, session){
  #function to take metadata input file
  load_meta <- reactive({
    if (!is.null(input$metaFP)){
      meta <- read_csv(input$metaFP$datapath)
      return(meta)}
    else{
      return(NULL)
    }
  })
  
  #function to produce table
  summ_table <- function(meta_tib){
    if (!is.null(input$metaFP)){
    summ_tib <- tibble("Columns" = colnames(meta_tib), "Type" = sapply(meta_tib, class),
                       "Mean" = sapply(meta_tib, mean, na.rm = TRUE), "SD(+/-)" = sapply(meta_tib, sd, na.rm = TRUE))
    return(summ_tib)}
    else{return(NULL)}
  }
  #function to produce AOD histogram
  plot_aod <- function(meta_tib){
    if (!is.null(input$metaFP)){
      histo <- ggplot(meta_tib, aes(Age_of_death))+
      geom_histogram(bins = 10, color = "black", fill = "coral1")+
      labs(title = 'Histogram of Age of Death')+
      xlab('Age of Death')+
      ylab('Count')+
      theme_bw()
      return(histo)}
    else{return(NULL)}
  }
  #function to produce RIN histogram
  plot_rin <- function(meta_tib){
    if (!is.null(input$metaFP)){
      histo <- ggplot(meta_tib, aes(RIN))+
        geom_histogram(bins = 10, color = "black", fill = "cadetblue2")+
        labs(title = 'Histogram of RIN')+
        xlab('RIN')+
        ylab('Count')+
        theme_bw()
      return(histo)}
    else{return(NULL)}
  }
  #function to produce PMI histogram
  plot_pmi <- function(meta_tib){
    if (!is.null(input$metaFP)){
      histo <- ggplot(meta_tib, aes(PMI))+
        geom_histogram(bins = 10, color = "black", fill = "orchid3")+
        labs(title = 'Histogram of PMI')+
        xlab('PMI')+
        ylab('Count')+
        theme_bw()
      return(histo)}
    else{return(NULL)}
  }
  #methods to render tables and graphs
  output$summ <- renderTable({
    summ_table(load_meta())
    })
  output$metadata <- DT::renderDataTable({
    load_meta()
    })
  output$aod <- renderPlot({
    plot_aod(load_meta())
    })
  output$rin <- renderPlot({
    plot_rin(load_meta())
  })
  output$pmi <- renderPlot({
    plot_pmi(load_meta())
  })
}

shinyApp(ui=ui, server = server)