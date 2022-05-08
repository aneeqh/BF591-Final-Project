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
      selectInput("metachoice", choices = c("Diagnosis"="Diagnosis", "PMI"="PMI","Age of death"="Age_of_death", "RIN"="RIN", "Sequence Reads"="Seq_reads"),
                  label = "Select metadata category", selected = "Age_of_death"),
      #insert gene search box here
      textInput("gene", label = "Enter gene to search for", placeholder = "ENSG00000000003.10"),
      #selectInput("plotType", label = "Choose what type of plot to make", choices = c("Bar", "Box", "Violin", "Beeswarm")),
      submitButton(text='Plot !')
    ),
    mainPanel(
      plotOutput("distroplot")
    )
  )
)

#back end code-
server <- function(input, output, session){
  #function to take norm counts input file
  load_counts <- reactive({
    if (!is.null(input$counts4FP)){
      counts <- read_csv(input$counts4FP$datapath)
      return(counts)}
    else{
      return(NULL)
    }
  })
  #function to take metadata input file
  load_meta <- reactive({
    if (!is.null(input$meta4FP)){
      meta <- read_csv(input$meta4FP$datapath)
      return(meta)}
    else{
      return(NULL)
    }
  })
  #function to make distribution plots-
  plot_distro <- function(counts_tib, meta_tib, meta_cat, selectgene){
    if (!is.null(input$meta4FP) & !is.null(input$counts4FP)){
      counts_tib <- column_to_rownames(counts_tib, var = "gene")
      gene_counts <- as.numeric(as.vector(counts_tib[selectgene,]))
      plot_tib <- tibble(Gene_Counts = gene_counts, meta_value = pull(meta_tib, meta_cat))
      if (meta_cat == "Diagnosis"){
        plot <- ggplot(plot_tib, aes(meta_value))+
          geom_bar()+
          theme_bw()+
          labs(title = "Plot of gene counts vs Diagnosis")
        return(plot)
      }
      else {
        plot <- ggplot(plot_tib, aes(meta_value,Gene_Counts))+
          geom_point()+
          theme_bw()+
          labs(title = str_c("Plot of gene counts vs ", meta_cat))
        return(plot)
      }}
    else{return(NULL)}
  }
  #methods to return outputs
  output$distroplot <- renderPlot({
    plot_distro(load_counts(), load_meta(), input$metachoice, input$gene)
  })
}

shinyApp(ui=ui, server = server)