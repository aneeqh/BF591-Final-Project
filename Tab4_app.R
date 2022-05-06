#import libraries
library(tidyverse)
library(shiny)
library(ggplot2)
library(DT)

#code for front end
ui <- fluidPage(
  titlePanel("Individual Gene Expresion Visaulization"),
  p("On this page you can visualize data about each gene individually."),
  sidebarLayout(
    sidebarPanel(
      fileInput(inputId = 'counts4FP', label = 'Load normalized counts matrix CSV'),
      fileInput(inputId = 'meta4FP', label = 'Load sample information matrix CSV'),
      selectInput("metachoice", choices = c("PMI","Age of death"="Age_of_death", "RIN", "Sequence Reads"="Seq_reads")),
      #insert gene search box here
      selectInput("plotType", label = "Choose what type of plot to make", choices = c("bar", "box", "violin", "beeswarm")),
      submitButton(text='Plot !')
    ),
    mainPanel(
      tabsetPanel(
        
      )
    )
  )
)

#back end code-
server <- function(input, output, session){
  
}

shinyApp(ui=ui, server = server)