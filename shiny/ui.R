#######################################################################################
############################# CARPnTrees web application ############################## 
############# building consensus tree from ancient and recent phylogenies #############  
# Marie Zufferey - UNIL - December 2015 ###############################################
# License: Open Source AGPL v3 ########################################################
#######################################################################################

#####################################
##### User interface - web app  #####
#####################################
#minimal code 
#headerPanel(), sidebarPanel(), mainPanel()

library(shiny)
library(shinythemes)
library(shinyBS)
library(shinyjs)
library(graphics)

shinyUI(fluidPage(    
  theme = shinytheme("cosmo"),
  HTML('<title>Compare TAD lists</title>
       <h1 style="color:#2780E3;font-size:400%;" align="center"><b>Compare TAD partitions </b></h1>
       <h4 align="center"><i>comparison of TAD caller outputs </i></h4>
       <hr>'),
  sidebarLayout(  
    sidebarPanel(
      
      checkboxInput('demo', 'Example with demo data', FALSE),
      
      tags$hr(),
      
      
      selectInput("simMetric", 
                  label = "Choose a comparison metric",
                  choices = list("Measure of Concordance (MoC)" = "moc", 
                                 "Jaccard Index (bins)" = "ji_bin",
                                 "Jaccard Index (boundaries)" = "ji_bd", 
                                 "Variation of information" = "voi"),
                  selected = "moc"),
      
      
      
      
      conditionalPanel(
        condition="input.demo == false",
        tags$hr(),
        
        fileInput("allDomainFiles", label=("Files with TAD coordinates")),
        tags$hr()
      ),
      
      
      
      
      conditionalPanel(
        condition="input.demo == false",
        
                    
          HTML('<u>Input format:</u>'),
                    
        
          checkboxInput('fileHeader', 'File with header', TRUE),
          radioButtons('fieldSep', 'Separator',
                       c(Comma=',',
                         Semicolon=';',
                         Tab='\t'),
                       ','),

        useShinyjs(),                                                  # Include shinyjs in the UI
        actionButton("reset_button", "Reset your choices")
      )
    ),
    mainPanel(
      tabsetPanel(
        ###### TAB WITH HELP #####
        tabPanel("Help",
                 htmlOutput("help")
        ),
        ###### TAB WITH THE OUTPUT #####
        
        # tabPanel("Your analysis",
        #          ##### Analysis summary
        #          h2("Your analysis:"),
        #          htmlOutput("your_analysis"),
        #          tags$br(),
        #          
        #          ##### Consensus tree
        #          h4("Heatmap:"),
        #          HTML('<table>
        #               <tr>
        #               <td rowspan="2">'),
        #          htmlOutput("partition_similarity"),
        #          HTML(   '</td>
        #                  <td>
        #                  <style>button#plotSimilarityBtn{width:200px;height:30px;padding-top:5px;border:0px;}</style>
        #                  <button id="plotSimilarityBtn" style="background-color:#2780E3" type="button" icon="fa fa-tree"
        #                  class="btn btn-default action-button"><i class="fa fa-tree"></i> Plot similarity heatmap</button>'),
        #          
        #          bsModal(id="plotPartitionSimilarity", title="Parition similarity", trigger="plotSimilarityBtn", size = "large",
        #                  plotOutput("consensusTree"), downloadLink("outputPlot", "Download")),
        #          HTML(    '</td>
        #                   </tr>
        #                   <tr>
        #                   <td>
        #                   <style>a#outputDwld{width:200px;height:30px;padding-top:5px;}</style>'),
        #          downloadButton('outputDwld', 'Download')
        #          ),
        
        
        
        ###### TAB WITH DATA OVERVIEW ##### # => add statitsitcs (# of domains, size of domains, etc.)
        tabPanel("Your data",
                 conditionalPanel(
                   h4("Descriptive stats. of current TAD lists")),
                 
                 tableOutput("showDescStat")
      
        ),
        
        ##### TAB WITH THE HEAD OF LOADED FILES
        
        tabPanel("Loaded data",
                 conditionalPanel(
                   condition="input.demo==false",
                   h4("First lines of the loaded files")),
                 conditionalPanel(
                   condition="input.demo==true",
                   h4("First lines of the demo files")),
                 tableOutput("showDataHead")
        ),
        
        tabPanel("Current parameters",
                 tableOutput("showCurrentOptions")
        ),
        
        
        
        
        ###### TAB WITH CONTACT #####
        tabPanel("Contact",
                 htmlOutput("contact")
            )
                                  ) # close tabsetPanel
                 )
    
                 ), # close mainPanel

  # ADD THE FOOTER WITH PICTURES FOR THE LOGO AND THE LICENSE TEXT
HTML('<hr><footer>
            <a href="http://www.unil.ch" target="_blank">
              <img src="images/logo_UNIL.png" alt="UNIL logo" width="5%" 
                  style="vertical-align:bottom;">
              </a><div style="text-align:center;margin-top:-21px">
              <a align=center href="https://opensource.org/licenses/AGPL-3.0" style="font-size:70%;"  target="_blank"> 
              Open Source AGPL v3 License </a></div>

<div style="text-align:right;margin-top:-21px">
        <a href="http://ciriellolab.org" target="_blank">

              <img src="images/logo_CSO.png" alt="CSO logo" width="5%" 
                  style="vertical-align:bottom;"></a><div style="text-align:center;margin-top:-21px">
</div>

              </p>
              </footer>')
        )  # close fluidPage
  ) # close shinyUI