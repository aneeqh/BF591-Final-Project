#import libraries
library(tidyverse)
library(shiny)
library(ggplot2)

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
       tabPanel('Metadata'),
       tabPanel('Plots')
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
  #methods to render tables and graphs
  output$summ <- renderTable({
    summ_table(load_meta())
    })
}

shinyApp(ui=ui, server = server)
