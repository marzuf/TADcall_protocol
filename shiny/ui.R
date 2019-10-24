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
      
      HTML('<u>Input format:</u>'),
      
      conditionalPanel(
        condition="input.demo == false",
        
        
        conditionalPanel(
          condition="input.demo == false",
          checkboxInput('headP', 'File with header', TRUE)
          
        ),
        
        conditionalPanel(
          condition="input.demo == false",
          radioButtons('sepP', 'Separator',
                       c(Comma=',',
                         Semicolon=';',
                         Tab='\t'),
                       ',')
        ),
        
        
        selectInput("var", 
                    label = "Choose a comparison metric",
                    choices = list("Measure of Concordance (MoC)" = "moc", 
                                   "Jaccard Index (bins)" = "ji_bin",
                                   "Jaccard Index (boundaries)" = "ji_bd", 
                                   "Variation of information" = "voi"),
                    selected = "moc"),
        
        
        
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
        
        tabPanel("Your analysis",
                 ##### Analysis summary
                 h2("Your analysis:"),
                 htmlOutput("your_analysis"),
                 tags$br(),
                 
                 ##### Consensus tree
                 h4("Consensus tree:"),
                 HTML('<table>
                      <tr>
                      <td rowspan="2">'),
                 htmlOutput("newick_consensus"),
                 HTML(   '</td>
                         <td>
                         <style>button#plotConsensusBtn{width:200px;height:30px;padding-top:5px;border:0px;}</style>
                         <button id="plotConsensusBtn" style="background-color:#2780E3" type="button" icon="fa fa-tree"
                         class="btn btn-default action-button"><i class="fa fa-tree"></i> Plot the tree</button>'),
                 
                 bsModal("plotConsensusTree", "Consensus tree", "plotConsensusBtn", size = "large",
                         plotOutput("consensusTree"), downloadLink("outputPlot", "Download")),
                 
                 HTML(    '</td>
                          </tr>
                          <tr>
                          <td>
                          <style>a#outputDwld{width:200px;height:30px;padding-top:5px;}</style>'),
                 downloadButton('outputDwld', 'Download'),
                 
                 HTML(      '</td>
                            </tr>
                            </table>
                            
                            <script>
                            hide = 0;
                            function changeText(id) {
                            if(hide==0){
                            id.innerHTML = "Hide ancient and recent phylogenies";
                            hide = 1;
                            }else {
                            id.innerHTML="Show ancient and recent phylogenies";
                            hide = 0;
                            }
                            }
                            </script>
                            <!--Button to display also ancient and recent polymorphism phylogenies -->
                            </br><style>button#AR{width:300px;height:25px;padding-top:4px; font-size:80%;border:0px;color:black;}</style>
                            <button id="AR" type="button" class="btn btn-default action-button" style="background-color:#E3E3E3" onclick="changeText(this)" value="show">
                            <i class="glyphicon glyphicon-tree-deciduous"> </i>
                            Show ancient and recent phylogenies</button></br></br>'),
                 
                 ##### If click: ancient and recent trees
                 conditionalPanel("hide==1",
                                  h4("Recent polymorphism phylogeny:"),
                                  HTML('<table>
                                       <tr>
                                       <td rowspan="2">'),
                                  htmlOutput("newick_recent"),
                                  HTML(    '</td>
                                           <td>
                                           <style>button#plotRecentBtn{width:200px;height:30px;padding-top:5px;border:0px;}</style>
                                           <button id="plotRecentBtn" style="background-color:#2780E3" type="button" icon="fa fa-tree"
                                           class="btn btn-default action-button"><i class="fa fa-tree"></i> Plot the tree</button>'),
                                  
                                  bsModal("plotRecentTree", "Recent polymorphism phylogeny", "plotRecentBtn", size = "large",
                                          plotOutput("recentTree"), downloadLink("outputPlot_recent", "Download")),
                                  
                                  HTML(    '</td>
                                           </tr>
                                           <tr>
                                           <td>
                                           <style>a#outputDwld_recent{width:200px;height:30px;padding-top:5px;}</style>'),
                                  downloadButton('outputDwld_recent', 'Download'),
                                  
                                  HTML(      '</td>
                                             </tr>
                                             </table>')
                                  
                                  ),
                 conditionalPanel("hide==1",
                                  h4("Ancient polymorphism phylogeny:"),
                                  #htmlOutput("newick_ancient")
                                  HTML('<table>
                                               <tr>
                                                   <td rowspan="2">'),
                                  htmlOutput("newick_ancient"),
                                  HTML(    '</td>
                                                   <td>
                                                    <style>button#plotAncientBtn{width:200px;height:30px;padding-top:5px;border:0px;}</style>
                                      <button id="plotAncientBtn" style="background-color:#2780E3" type="button" icon="fa fa-tree"
                                       class="btn btn-default action-button"><i class="fa fa-tree"></i> Plot the tree</button>'),
                                  bsModal("plotAncientTree", "Ancient polymorphism phylogeny", "plotAncientBtn", size = "large",
                                          plotOutput("ancientTree"),  downloadLink("outputPlot_ancient", "Download")),
                                  HTML(    '</td>
                                                 </tr>
                                                     <tr>
                                                         <td>
                                                         <style>a#outputDwld_ancient{width:200px;height:30px;padding-top:5px;}</style>'),
                                  downloadButton('outputDwld_ancient', 'Download'),
                                  HTML(      '</td>
                                                      </tr>
                                                        </table>')                             
                 )
                 
                                  ),
        ###### TAB WITH DATA OVERVIEW #####
        tabPanel("Your data",
                 conditionalPanel(
                   condition="input.haploData == true || input.demo==true",
                   h4("First lines of your recent polymorphism data")),
                 
                 conditionalPanel(
                   condition="input.haploData == false & input.demo==false",
                   h4("First lines of your polymorphism and haplogroup data")),
                 
                 tableOutput("showDataTableP"),
                 
                 conditionalPanel(
                   condition="input.haploData == true || input.demo==true",
                   h4("First lines of your haplogroup data"),
                   tableOutput("showDataTableH"))
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