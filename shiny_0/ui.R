library("shiny")

datasetList <- gsub("_final_domains.txt", "", list.files("~/Documents/CSO/shinyTrial/data/domains/"))

# Define UI for application that draws a histogram
shinyUI(fluidPage(
  
  # Application title
  titlePanel("Compare TAD callers"),
  
  # Sidebar with a slider input for the number of bins
  sidebarLayout(
    sidebarPanel(
      ### SELECT INPUT DATASETS
      checkboxGroupInput("av_dataset", 
                         label = h3("Dataset to select"), 
                         # choices = as.list(setNames(c(1:length(datasetList)), datasetList)), # the values can be arbitrary (stored in vector input$av_dataset)
                         choices = as.list(setNames(datasetList, paste(datasetList, "domains"))), # the values can be arbitrary (stored in vector input$av_dataset)
                         selected = NULL), # should correspond to one value of the above vector
      
      checkboxGroupInput("av_plot", 
                         label = h3("Plots to display"), 
                         # choices = as.list(setNames(c(1:length(datasetList)), datasetList)), # the values can be arbitrary (stored in vector input$av_dataset)
                         choices = list("Chromosome coverage" = "cover", "Number of TADs" = "nbr"),
                         selected = NULL), # should correspond to one value of the above vector
      
      
      HTML("<br><h3> MoC for pair of datasets </h3>"),
      
      selectizeInput("selectForMoC1", 
                         label = NULL, 
                         # choices = as.list(setNames(c(1:length(datasetList)), datasetList)), # the values can be arbitrary (stored in vector input$av_dataset)
                         choices = as.list(setNames(datasetList, paste(datasetList, "domains"))),
                         options = list(
                                    placeholder = 'Please select a caller from the list below',
                                    onInitialize = I('function() { this.setValue(""); }')
                                  ),
                         selected = NULL), # should correspond to one value of the above vector
      selectizeInput("selectForMoC2", 
                  label = NULL, 
                  # choices = as.list(setNames(c(1:length(datasetList)), datasetList)), # the values can be arbitrary (stored in vector input$av_dataset)
                  choices = as.list(setNames(datasetList, paste(datasetList, "domains"))),
                  options = list(
                    placeholder = 'Please select a caller from the list below',
                    onInitialize = I('function() { this.setValue(""); }')
                  ),                  
                  selected = NULL), # should correspond to one value of the above vector
      actionButton("computeMoC", 
                   label = "Get MoC"
                   ),
      HTML("<br>"),
      textOutput("MoCpending"),
      HTML("<br>"),
      uiOutput("MoCvalue")
      
      
      
      
    ), #-close sidebarPanel
    
    # Show a plot of the generated distribution
    mainPanel(
      #*************************************** CHROMOSOME COVERAGE
      ### HEADER COVERAGE
      uiOutput("headerCoverage"),
      ### PLOT COVERAGE
      plotOutput("coveragePlot"),
      
      #*************************************** NUMBER
      ### HEADER NUMBER
      uiOutput("headerNbr"),
      ### PLOT NUMBER
      plotOutput("nbrPlot")
      
      
    ) #-close mainPanel
  )
))