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


demo_files <- list.files("data", pattern="_domains.txt", full.names = TRUE)



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
  

  
  
  ##############################################################-----------------------------------------------------------------------------------
  ### Tab with data information (head of the loaded files)   ###-----------------------------------------------------------------------------------
  ##############################################################-----------------------------------------------------------------------------------
  
  output$showDataHead <- renderTable({
    if(input$demo){
      file_heads <- lapply(demo_files, function(x) head(read.table(x, header=F, stringsAsFactors = FALSE, col.names=c("chromo", "start", "end"))))
      names(file_heads) <- gsub("_final_domains.txt", "", basename(demo_files))
      # print(file_heads)
      # for(i in file_heads) {write.table(i, file="");cat("\n")}
      file_heads <- lapply(demo_files, function(x) print(head(read.table(x, header=F, stringsAsFactors = FALSE, col.names=c("chromo", "start", "end")))))
      
    }else{
      inFile <- NULL
      print(inFile)
    }
  })
  
  
  #######################################-----------------------------------------------------------------------------------
  ### Tab with current options        ###-----------------------------------------------------------------------------------
  #######################################-----------------------------------------------------------------------------------
  
  output$showCurrentOptions <- renderUI({
    
    
    HTML(paste("<u>Similarity metric</u>"))
    
    
    # getDomainFiles <- function(){
    #   if(input$demo){
    #     all_domain_files <- demo_files
    #   }else {
    #     all_domain_files <- input$allDomainFiles$datapath
    #   }
    #   return(all_domain_files)
    # }
    # getFileHeader <- function(){
    #   fileHead <- ifelse(input$demo, FALSE, input$fileHeader)
    #   return(fileHead)
    # }
    # getFieldSep <- function(){
    #   sep <- ifelse(input$demo, "\t", input$fieldSep)
    #   return(sep)
    # }
    # 
    # HTML(paste("<u>Input files</u>"))
    # HTML(paste0(getDomainFiles(), collapse=", "))
    # tags$hr()
    # 
    # HTML(paste("<u>Field separator</u>"))
    # HTML(paste0(getFieldSep(), collapse=", "))
    # 
    # tags$hr()
    # 
    # HTML(paste("<u>File with header ?</u>"))
    # HTML(paste0(getFileHeader(), collapse=", "))
    # 
    # tags$hr()
    # 
    # HTML(paste0(input$simMetric, collapse=", "))  
    
    })
  
  
  ######################################-----------------------------------------------------------------------------------
  ### Tab with the desc. stats      ###-----------------------------------------------------------------------------------
  ######################################-----------------------------------------------------------------------------------
  output$showDescStat <- renderUI({
    
    HTML(paste("<u>Plot with desc. stats</u>"))
    
    
  })
  

  ######################################-----------------------------------------------------------------------------------
  ### Tab with the output trees      ###-----------------------------------------------------------------------------------
  ######################################-----------------------------------------------------------------------------------
  
  ##### Info about the analyses
  #*************************************************
  output$your_analysis <- renderUI({
    
    T1 <- "Data: "
    t1 <- ("for NRY.")
    
    HTML(paste("<ul><li><u>",T1,"</u>&nbsp&nbsp",t1,"</li>"))
  })
  ###
  ##### Functions needed to set parameters before compute trees #####
  ###
  getDomainFiles <- function(){
    if(input$demo){
      all_domain_files <- demo_files
    }else {
      all_domain_files <- input$allDomainFiles$datapath
    }
    return(all_domain_files)
  }
  getFileHeader <- function(){
    fileHead <- ifelse(input$demo, FALSE, input$fileHeader)
    return(fileHead)
  }
  getFieldSep <- function(){
    sep <- ifelse(input$demo, "\t", input$fieldSep)
    return(sep)
  }

  ##### The similarity heatmap
  #*************************************************
  heatmapSimilarity <- reactive({
    
    return(list("plot1", "plot2"))
    
    domainFiles <- getDomainFiles()
    
    all_TAD_dt <- lapply(domainFiles, function(x) {
      i_dt <- read.delim(x, header=getFileHeader(), sep=getFieldSep(), stringsAsFactors = FALSE)
      stopifnot(ncol(i_dt) == 3)
      colnames(i_dt) <- c("chromo", "start", "end")
      stopifnot(is.numeric(i_dt$start))
      stopifnot(is.numeric(i_dt$end))
    })

    # if(is.na(domainFiles)){
      if(length(domainFiles) == 0){
      simHeat <- "Waiting for your data !"
    } else {
      simHeat <- try(plot(NULL,xlim=c(0,1),ylim=c(0,1))) # ADD HERE THE PLOT SIMILARITY FUNCTION  
    }

    if(class(simHeat)=="try-error"){
      simHeat <- "Error during similarity calculation and plotting ! Did you select the right parameters ?"
    }
    
    return(simHeat)
  })
  

  
  ### Ouput text
  output$partition_similarity <- renderUI({
    HTML(paste0("<style> textarea#output{width:500px;height:65px;}</style>
                <textarea id=\"output\">", heatmapSimilarity()[[2]],"</textarea>"))
  })
  
  ### Plot in popup window  
  output$heatmapSimilarity <- renderPlot({
    heatPlot <- heatmapSimilarity()[[1]]
    # detach(package:shinyjs)
    a <- try(plot(heatPlot) )
    if(class(a)=="try-error"){
      print("not available")
    }
    # require(shinyjs)    
  })
  
  ###
  ##### Functions to download the similarity values TXT #####
  ###
  ##### 

  output$outputDwld_similarity <- downloadHandler(
    filename="caller_similarity.txt",
    content = function(file){
      write.table(heatmapSimilarity()[[2]], file=file, sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
    }
  )
  
  ###
  ##### Functions to download the similarity heatmap PNG #####
  ###
  ##### 
  output$outputPlot_similarity <- downloadHandler(
    filename="similarity_heatmap",
    content = function(file){
      p <- heatmapSimilarity()[[1]]
      ggsave(p, file = file)
    }
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
    #     reset("partition_similarity")
    #     reset("newick_recent")
    #     reset("newick_ancient")        
  })
})

