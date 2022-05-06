#import libraries
library(tidyverse)
library(shiny)
library(ggplot2)
library(DT)

ui <- fluidPage(
  titlePanel("Differential Expression results"),
  p("On this page you can upload the results of differential expression analysis and explore them."),
  p("To use this application, download the CSV 'deseq_res.csv' from the data directory of this app's repository."),
  sidebarLayout(
    sidebarPanel(
      fileInput(inputId = 'DEFP', label = 'Load differential expression results CSV'),
      br(),
      p('A volcano plot can be generated with "log2 fold-change" on the x-axis and "p-adjusted" on the y-axis'),
      br(),
      radioButtons(inputId = 'xaxis', label = 'Choose the column for the x-axis', choices = c('baseMean','log2FoldChange', 'lfcSE','stat','pvalue','padj'), 
                   selected = 'log2FoldChange'),
      radioButtons(inputId = 'yaxis', label = 'Choose the column for the x-axis', choices = c('baseMean','log2FoldChange', 'lfcSE','stat','pvalue','padj'), 
                   selected = 'padj'),
      colourpicker::colourInput(inputId = 'basecol', label='Base point color', value='#22577A'),
      colourpicker::colourInput(inputId = 'highcol', label='Highlight point color', value='#FFCF56'),
      sliderInput(inputId = 'padjmag', label = 'Select the magnitude of the p adjusted coloring:', min = -50, max = 0, value = -10,),
      submitButton(text='Plot')
    ),
    mainPanel(
      tabsetPanel(
        tabPanel('Table', dataTableOutput("difftable")),
        tabPanel('Volcano Plot', plotOutput("volcano")),
        tabPanel('Plot table', dataTableOutput("plotTable"))
    )
  )
))
  
server <- function(input, output, session){
  options(shiny.maxRequestSize=30*1024^2)
  #function to load input
  load_de <- reactive({
    if (!is.null(input$DEFP)){
      defp <- read_csv(input$DEFP$datapath)
      return(defp)}
    else{return(NULL)}
    })
  #function to produce volcano plot
  volcano_plot <-
    function(dataf, x_name, y_name, slider, color1, color2) {
      if (!is.null(input$DEFP)){
        #make plot df-
        cols <- c("FALSE" = color1, "TRUE" = color2)
        threshold <- 10^slider
        plot_data <- dataf %>%
          dplyr::mutate(thresh = case_when(padj <= threshold ~ "TRUE", TRUE ~ "FALSE"))
        #produce volcano plot
        volc_plot <- ggplot( plot_data, aes(get(x_name), -log10(get(y_name))))+
          geom_point(aes(color=thresh))+
          scale_color_manual(values = cols)+
          theme_bw()+
          theme(legend.position = "bottom", aspect.ratio=1.1)
        return(volc_plot)}
      else{return(NULL)}
    }
  #function to produce plot table
  draw_table <- function(dataf, slider) {
    if (!is.null(input$DEFP)){
      df <- dataf %>% 
        dplyr::filter(padj <= 10^slider)
      df[] <- lapply(df, formatC, format="f", digits = 4)
      return(df)}
    else{return(NULL)}
  }
  #section to return outputs to App
  output$difftable <- DT::renderDataTable({
    load_de()
  })
  output$volcano <- renderPlot({
    volcano_plot(load_de(), input$xaxis, input$yaxis, input$padjmag, input$basecol, input$highcol)
  })
  output$plotTable <- renderDataTable({
    draw_table(load_de(), input$padjmag)
  })
}

shinyApp(ui=ui, server = server)