library(shiny)
library(scales)
library(ggplot2)

datasetList <- gsub("_final_domains.txt", "", list.files("~/Documents/CSO/shinyTrial/data/domains"))

mainDir <- paste0("~/Documents/CSO/shinyTrial/data")

chrSize <- read.delim(paste0(mainDir, "/input/chr6.size"), header=F)[1,2]

pyexec <- "python"

scriptMoC <- "/media/electron/mnt/ed2/shared/TADcompare/pipeline/src/printMoC.py"

coverColor <- "sienna4"
nbrColor <- "navy"

# Define server logic required to draw a histogram
shinyServer(function(input, output, session) {

  #***************************************************************** CHROMOSOME COVERAGE
  ##### CHROMOSOME COVERAGE TITLE
  output$headerCoverage <- renderUI({
    if(! "cover" %in% input$av_plot)
      return(NULL)
    if(length(input$av_dataset) == 0)
      return(NULL)
    HTML("<h3> Chromosome coverage by TADs</h3>")
  })
  ##### PLOT CHROMOSOME COVERAGE
  output$coveragePlot <- renderPlot({
    if(! "cover" %in% input$av_plot)
      return(NULL)
    if(length(input$av_dataset) == 0)
      return(NULL)
    covDT <- foreach(ds = input$av_dataset, .combine="rbind") %do% {
      # callers <- c("DI", "ICFinder", "TopDom")
      # covDT <- foreach(ds = callers, .combine="rbind") %do% {
      callerDT <- read.delim(paste0(mainDir, "/domains/", ds, "_final_domains.txt"), header=F, stringsAsFactors = F, col.names=c("chromo", "start", "end"))
      data.frame(caller = ds, coverage = sum(callerDT$end-callerDT$start+1)/chrSize)
    }
    barplot(covDT$coverage, axes=F, col =coverColor, names.arg = covDT$caller, axisnames=T, main = "Ratio of chromosome covered by TADs", ylab="chromosome ratio")
    axis(2, at=seq(0,1,0.2))
    box(bty="l")
  })
  
  
  #***************************************************************** NUMBER
  ##### CHROMOSOME COVERAGE TITLE
  output$headerNbr <- renderUI({
    if(! "nbr" %in% input$av_plot)
      return(NULL)
    if(length(input$av_dataset) == 0)
      return(NULL)
    HTML("<h3> Number of TADs on the chromosome</h3>")
  })
  ##### PLOT CHROMOSOME COVERAGE
  output$nbrPlot <- renderPlot({
    if(! "nbr" %in% input$av_plot)
      return(NULL)
    if(length(input$av_dataset) == 0)
      return(NULL)
    covDT <- foreach(ds = input$av_dataset, .combine="rbind") %do% {
      # callers <- c("DI", "ICFinder", "TopDom")
      # covDT <- foreach(ds = callers, .combine="rbind") %do% {
      callerDT <- read.delim(paste0(mainDir, "/domains/", ds, "_final_domains.txt"), header=F, stringsAsFactors = F, col.names=c("chromo", "start", "end"))
      data.frame(caller = ds, number = nrow(callerDT))
    }
    barplot(covDT$number, axes=F, col =nbrColor, names.arg = covDT$caller, axisnames=T, main = "Number of TADs on chromosome",ylab="# of TADs")
    axis(2)
    box(bty="l")
  })
  
  
  #***************************************************************** MOC CALCULATION
  
  output$MoCvalue <- renderPrint({
    if(input$selectForMoC1 == "" | input$selectForMoC2 == "")
      return(cat(""))
    if(input$computeMoC == 0)
      return(cat(""))
    if(input$computeMoC == 1) {
      file1 <- paste0(mainDir, "/domains/", input$selectForMoC1, "_final_domains.txt")
      file2 <- paste0(mainDir, "/domains/", input$selectForMoC2, "_final_domains.txt")
      cmd <- paste(pyexec, scriptMoC, "-f1", file1, "-f2", file2, "-c", chrSize)
      MoC  <- round(as.numeric(as.character(system(cmd, intern=T))), 2)
      return(cat(paste0(input$selectForMoC1, " - ", input$selectForMoC2, " MoC = ", MoC, "\n")))
    }
  })
  
  # increase the counter
  observe({
    output$MoCpending <- renderText({
      if(input$computeMoC == 1)
        return("*******************")
    })
  })
    
  # increase the counter
  observe({
    if(input$computeMoC == 0){
      # updateActionButton(session, "computeMoC")
    }else if(input$computeMoC == 1) {
      updateActionButton(session, "computeMoC", label = "... done")
    }
  })

  
}) #end-Shiny server


