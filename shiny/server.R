#######################################################################################
############################# CARPnTrees web application ############################## 
############# building consensus tree from ancient and recent phylogenies #############  
# Marie Zufferey - UNIL - December 2015 ###############################################
# License: Open Source AGPL v3 ########################################################
#######################################################################################

###########################################
##### Server implementation - web app #####
###########################################

# server implementation: we call shinyServer and pass it
# a function that accepts 2 parameters: input and output

library(shiny)
library(shinyBS)
library(ape)
library(shinyjs)
library(graphics)
source("app_consensus.R")
source("app_recent.R")
source("app_ancient.R")
shinyServer(function(input, output){
  
  ###############################-----------------------------------------------------------------------------------
  ### Tab with help html file ###-----------------------------------------------------------------------------------
  ###############################-----------------------------------------------------------------------------------
  
  # function to get the file with help information
  getHelp <- function(){
    return(includeHTML("help.html"))
  }
  # display the help in the tab
  output$help <- renderUI({
    getHelp()
  })
  
  ##################################-----------------------------------------------------------------------------------
  ### Tab with contact html file ###-----------------------------------------------------------------------------------
  ##################################-----------------------------------------------------------------------------------
  
  # function to get the file with help information
  getContact <- function(){
    return(includeHTML("contact.html"))
  }
  # display the help in the tab
  output$contact <- renderUI({
    getContact()
  })
  
  ##############################################-----------------------------------------------------------------------------------
  ### Tab with data information              ###-----------------------------------------------------------------------------------
  ##############################################-----------------------------------------------------------------------------------
  
  output$showDataTableP <- renderTable({
    if(input$demo){
      inFile <- "demo_seq.fasta"
      analysis <- "seq"  
    }else{
      inFile <- input$polyFile$datapath
      analysis <-input$analysis
    }
    if(is.null(inFile)){
      return(NULL)
    }else{
      if(analysis=="seq"){
        head(read.table(inFile, header=F,
                        col.names="FASTA_file"))          
      }else if(analysis=="str"){
        head(read.table(inFile, header=input$headP,
                        sep=input$sepP))
      }
    }
  })
  
  output$showDataTableH <- renderTable({
    if(input$demo){
      inFile <- "demo_HG.csv"
      analysis <- "seq"  
      sepH <- ","
      headH <- TRUE
    }else{
      inFile <- input$haploFile$datapath
      analysis <-input$analysis
      headH <- input$headH
      sepH <- input$sepH
    }
    if(is.null(inFile)){
      return(NULL)
    }else{
      head(read.table(inFile, header=headH,
                      sep=sepH))	    
    }
  })
  
  ######################################-----------------------------------------------------------------------------------
  ### Tab with the output trees      ###-----------------------------------------------------------------------------------
  ######################################-----------------------------------------------------------------------------------
  
  ##### Info about the analyses
  #*************************************************
  output$your_analysis <- renderUI({
    
    if(input$MorP=="mtDNA" | input$demo){
      t1 <- ("for mtDNA.")
    } else if (input$MorP=="NRY"){
      t1 <- ("for NRY.")
    }
    T1 <- "Data: "
    
    if(input$analysis=="seq" | input$demo){
      t2 <- ("aligned sequences (FASTA format).")
    }else if(input$analysis=="str"){
      t2 <- ("microsatellites (data frame).")
    }
    T2 <- "Recent polymorphism data: "
    
    if(input$haploData | input$demo){
      t3 <- ("in separate file.")
    }else{
      t3 <- ("in same file as recent polymorphism data.")
    }
    T3 <- "Ancient polymorphism data: "
    
    if(input$allHaploData){
      t4 <- ("own haplogroup tree.")
    }else if(!input$allHaploData | input$demo){
      t4 <- ("provided haplogroup tree.")
    }
    T4 <- "Haplogroup tree: "
    
    HTML(paste("<ul><li><u>",T1,"</u>&nbsp&nbsp",t1,"</li>",        
               "<li><u>",T2,"</u>&nbsp&nbsp",t2,"</li>",        
               "<li><u>",T3,"</u>&nbsp&nbsp",t3,"</li>",
               "<li><u>",T4,"</u>&nbsp&nbsp",t4,"</li></ul>"))
  })
  ###
  ##### Functions needed to set parameters before compute trees #####
  ###
  getHGtree <- function(){
    if(input$allHaploData){
      allHaploTreeName <- input$allHaploFile$datapath
    }else if(!input$allHaploData){
      allHaploTreeName <- ifelse(input$MorP=="mtDNA", "allHaploTreeM.newick", NA)
    }
    return(allHaploTreeName)
  }
  
  getHeadP <- function(){
    headP <- ifelse(input$analysis=="seq"|input$demo, NA, input$headP)    
    return(headP)
  }
  
  getSepP <- function(){
    sepP <- ifelse(input$analysis=="seq"|input$demo, NA, input$sepP)
    return(sepP)
  }
  
  getHaplo <- function(){
    haploFile <- ifelse(is.null(input$haploFile), NA, input$haploFile$datapath)
    if(input$demo){
      haploFile <- "demo_HG.csv"
    }
    return(haploFile)
  }
  
  getHeadH <- function(){
    headH <- ifelse(!input$haploData, getHeadP(), input$headH)
    if(input$demo){
      headH <- TRUE
    }
    return(headH)
  }
  
  getSepH <- function(){
    sepH <- ifelse(!input$haploData, getSepP(), input$sepH)
    if(input$demo){
      sepH <- ","
    }
    return(sepH)
  }
  
  getPoly <- function(){
    polyFile <- ifelse(is.null(input$polyFile), NA, input$polyFile$datapath)
    if(input$demo){
      polyFile <- "demo_seq.fasta"
    }
    return(polyFile)
  }
  
  ##### The consensus tree
  #*************************************************
  consensusTree <- reactive({
    headP <- getHeadP()
    sepP <- getSepP()
    headH <- getHeadH()
    sepH <- getSepH()
    haploFile <- getHaplo()
    polyFile <- getPoly()
    allHaploTree <- getHGtree()
    if(is.na(polyFile)){
      tree <- "Waiting for your data !"
    }else if(input$analysis=="seq" & is.na(haploFile) ){
      tree <- "Waiting for your haplogroup data !"
    } else if(input$allHaploData & is.null(input$allHaploFile)){
      tree <- "Waiting for your data !"
    }else if(!input$allHaploData & is.na(allHaploTree) & !input$demo){
      tree <- "Provided haplogroup tree for Y-chromosome not yet available..."
    } else if(input$haploData & is.na(haploFile)){
      tree <- "Waiting for your data !"
    }else{
      if(input$demo){
        tree <- try(CARP_main_consensus("demo_seq.fasta", "demo_HG.csv", "seq",  "allHaploTreeM.newick", NA, NA,TRUE,","))
      }else{
        tree <- try(CARP_main_consensus(polyFile, haploFile, input$analysis, allHaploTree, headP, sepP,headH, sepH))
      }
    }
    if(class(tree)=="try-error"){
      tree <- "Error during tree building ! Did you select the right parameters ?"
    }
    return(tree)
  })
  
  ### Ouput text
  output$newick_consensus <- renderUI({
    HTML(paste0("<style> textarea#output{width:500px;height:65px;}</style>
                <textarea id=\"output\">", consensusTree(),"</textarea>"))
  })
  ### Plot in popup window  
  output$consensusTree <- renderPlot({
    t <- consensusTree()
    tree <- read.tree(text=t)
    detach(package:shinyjs)
    a <- try(plot.phylo(tree) )
    if(class(a)=="try-error"){
      print("not available")
    }
    require(shinyjs)    
  })
  
  
  ##### Recent polymorphism tree (optional)
  #*************************************************
  recentTree <- reactive({
    headP <- getHeadP()
    sepP <- getSepP()
    haploFile <- getHaplo()
    polyFile <- getPoly()
    if(is.na(polyFile)){
      tree <- "Waiting for your data !"
    }
    else{
      if(input$demo){
        tree <- try(CARP_main_recent("demo_seq.fasta", "demo_HG.csv", "seq", NA, NA))
      }else{
        tree <- try(CARP_main_recent(polyFile, haploFile, input$analysis, headP, sepP))
      }
    }
    if(class(tree)=="try-error"){
      tree <- "Error during tree building ! Did you select the right parameters ?"
    } 
    return(tree)
  })
  ### Ouput text
  output$newick_recent <- renderUI({     
    HTML(paste0("<style> textarea#output{width:500px;height:65px;}</style>
                <textarea id=\"output\">",recentTree(),"</textarea>"))
  })  
  
  ### Plot in popup window
  output$recentTree <- renderPlot({
    t <- recentTree()
    tree <- read.tree(text=t)
    detach(package:shinyjs)
    plot.phylo(tree)  
    require(shinyjs)
  })
  
  ##### Ancient polymorphism tree (optional)
  #*************************************************
  ancientTree <- reactive({
    headP <- getHeadP()
    sepP <- getSepP()
    headH <- getHeadH()
    sepH <- getSepH()
    haploFile <- getHaplo()
    polyFile <- getPoly()
    allHaploTree <- getHGtree()
    if(is.na(polyFile) & is.na(haploFile)){
      tree <- "Waiting for your data !"
    }else if(input$analysis=="seq" & is.na(haploFile)){
      tree <- "Waiting for your haplogroup data !"
    }else if(input$allHaploData & is.null(input$allHaploFile)){
      tree <- "Waiting for your data !"
    }else if(!input$allHaploData & is.na(allHaploTree) & !input$demo){
      tree <- "Provided haplogroup tree for Y-chromosome not yet available..."
    } else{
      if(input$demo){
        tree <- try(CARP_main_ancient("demo_seq.fasta", "demo_HG.csv", "seq",  "allHaploTreeM.newick", NA, NA,TRUE,","))
      }else{
        tree <- try(CARP_main_ancient(polyFile, haploFile, input$analysis, allHaploTree, headP, sepP, headH, sepH))
      }
    }
    if(class(tree)=="try-error"){
      tree <- "Error during tree building ! Did you select the right parameters ?"
    }
    return(tree)
  })
  
  ancientTreePlot <- reactive({   # use this function if want to plot with diff. length and bold HG name
    headP <- getHeadP()
    sepP <- getSepP()
    headH <- getHeadH()
    sepH <- getSepH()
    haploFile <- getHaplo()
    polyFile <- getPoly()
    allHaploTree <- getHGtree()
    CARP_plot_ancient(polyFile, haploFile, input$analysis, allHaploTree, headP, sepP, headH, sepH)
  })
  
  ### Ouput text 
  output$newick_ancient <- renderUI({  
    HTML(paste0("<style> textarea#output{width:500px;height:65px;}</style>
                <textarea id=\"output\">", ancientTree(),"</textarea>"))
  })  
  
  ### Plot in popup window
  output$ancientTree <- renderPlot({
    #other possibility:
    #ancientTreePlot()  # call the function directly plot ancient polymorphism phylogeny with bold HG and diff. length
    t <- ancientTree()
    tree <- read.tree(text=t)
    detach(package:shinyjs)
    plot.phylo(tree)  
    require(shinyjs)
  })
  
  ###
  ##### Functions to download the output NEWICK #####
  ###
  ##### For consensus tree
  output$outputDwld <- downloadHandler(
    filename="consensus_tree.newick",
    content = function(file){
      writeLines(consensusTree(), file)
    }
  )
  
  ##### For recent polymorphsim tree
  output$outputDwld_recent <- downloadHandler(
    filename="recent_tree.newick",
    content = function(file){
      writeLines(recentTree(), file)
    }
  )
  
  ##### For ancient polymorphsim tree
  output$outputDwld_ancient <- downloadHandler(
    filename="ancient_tree.newick",
    content = function(file){
      writeLines(ancientTree(), file)
    }
  )
  
  ###
  ##### Functions to download the output PNG #####
  ###
  ##### For consensus tree
  
  output$outputPlot <- downloadHandler(
    filename="consensus_tree.png",
    content = function(file){
      png(file)
      t <- consensusTree()
      tree <- read.tree(text=t)
      detach(package:shinyjs)
      plot.phylo(tree)  
      dev.off()
      require(shinyjs)
    }
  )
  
  ##### For recent polymorphsim tree
  output$outputPlot_recent <- downloadHandler(
    filename="recent_tree.png",
    content = function(file){
      png(file)
      t <- recentTree()
      tree <- read.tree(text=t)
      plot.phylo(tree)  
      dev.off()
    }
  )
  
  ##### For ancient polymorphsim tree
  output$outputPlot_ancient <- downloadHandler(
    filename="ancient_tree.png",
    content = function(file){
      png(file)
      t <- ancientTree()
      tree <- read.tree(text=t)
      plot.phylo(tree)  
      dev.off()    }
  )
  
  ##### Reset button
  observeEvent(input$reset_button, {
    reset("analysis")
    reset("MorP")
    #     reset("polyFile")
    #     reset("haploFile")
    #     reset("haploData")
    #     reset("AR")
    #     reset("allHaploData")
    #     reset("newick_consensus")
    #     reset("newick_recent")
    #     reset("newick_ancient")        
  })
})

