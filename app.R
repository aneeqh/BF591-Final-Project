#import libraries
library(tidyverse)
library(shiny)
library(ggplot2)
library(DT)
library(bslib)

#front-end
ui <- fluidPage(
  theme = bs_theme(version = 4, bootswatch = "minty"),
  titlePanel("BF 591-R Final project"),
  p("On this website you can run analyses on the Huntington's disease dataset. Use the data in the data folder to run various analyses."),
  p("Developed by Aneeq Husain"),
  tabsetPanel(
    tabPanel("Samples", 
             h3("Sample Information Exploration"),
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
             )),
    tabPanel("Counts",
             h3("Counts Matrix Exploration"),
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
             )),
    tabPanel("Diff. Expr.",
             h3("Differential Expression Results"),
             p("On this page you can upload the results of differential expression analysis and explore them."),
             p("To use this application, upload a CSv file contaainig the differential expression data."),
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
             )),
    tabPanel("Indiv. Gene Expr.",
             h3("Individual Gene Expresion Visaulization"),
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
             ))
  )
)

server <- function(input, output, session){
  options(shiny.maxRequestSize=30*1024^2)  #increase file upload limit
  #function to take metadata input file for first tab
  load_meta <- reactive({
    if (!is.null(input$metaFP)){
      meta <- read_csv(input$metaFP$datapath)
      return(meta)}
    else{return(NULL)}
  })
  #function to take norm counts input file for second tab
  load_counts <- reactive({
    if (!is.null(input$countsFP)){
      counts <- read_csv(input$countsFP$datapath)
      return(counts)}
    else{return(NULL)}
  })
  #function to load Diff. Expr. results file
  load_de <- reactive({
    if (!is.null(input$DEFP)){
      defp <- read_csv(input$DEFP$datapath)
      return(defp)}
    else{return(NULL)}
  })
  #function to take inputs for fourth tab-
  load_4meta <- reactive({
    if (!is.null(input$meta4FP)){
      meta <- read_csv(input$meta4FP$datapath)
      return(meta)}
    else{return(NULL)}
  })
  load_4counts <- reactive({
    if (!is.null(input$counts4FP)){
      counts <- read_csv(input$counts4FP$datapath)
      return(counts)}
    else{return(NULL)}
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
  #function to produce volcano plot
  volcano_plot <- function(dataf, x_name, y_name, slider, color1, color2) {
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
  #methods to display output to front-end
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
  output$difftable <- DT::renderDataTable({
    load_de()
  })
  output$volcano <- renderPlot({
    volcano_plot(load_de(), input$xaxis, input$yaxis, input$padjmag, input$basecol, input$highcol)
  })
  output$plotTable <- renderDataTable({
    draw_table(load_de(), input$padjmag)
  })
  output$distroplot <- renderPlot({
    plot_distro(load_4counts(), load_4meta(), input$metachoice, input$gene)
  })
}

shinyApp(ui=ui, server = server)