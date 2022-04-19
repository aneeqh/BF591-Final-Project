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
  #function to take norm counts input file
  load_counts <- reactive({
    if (!is.null(input$countsFP)){
      meta <- read_csv(input$countsFP$datapath)
      return(meta)}
    else{
      return(NULL)
    }
  })
  
}

shinyApp(ui=ui, server = server)