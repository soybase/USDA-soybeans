library(shiny)

# Define UI for page that allows selection of genetic lines with corresponding facets
shinyUI(fluidPage(
  titlePanel("Soybean Copy Number Variation",
              tags$head(tags$style(type="text/css", "label.checkbox { display: inline-table; width: 17%; margin-right:2%;} label.radio { display: inline-table; width: 18%; margin-right:2%;}")) # format checkbox inputs in 3 columns with some padding.
              ),
  
  # Sidebar with checkbox inputs for varieties and chromosomes, plus a reset button.
  fluidRow( 
    column(width=3, 
           wellPanel(
    
    
    tagList(
      tags$table(width='100%',
                 tags$tr(
                   tags$td( # Reset Button
                     conditionalPanel(condition="input.plottype!='Methodology'",
                                      actionButton("resetButton", "Clear Selections")), 
                     style="text-align:left;"
                     ),
                   tags$td( # Download Button for GlymaID Table
                     conditionalPanel(condition="input.plottype=='CNV List'",
                                      downloadButton("DataFrameDownload", 
                                                     label="Download", class=NULL)),
                     style="text-align:right;"
                     ),
                   tags$td( # Download Button for Search by Location Table
                     conditionalPanel(condition="input.plottype=='Search CNVs by Location'",
                                      downloadButton("DataFrameDownload2", 
                                                     label="Download", class=NULL)),
                     style="text-align:right;"),
                   tags$td( # Generations Slider for Genealogy
                     div(conditionalPanel(condition="input.plottype=='Genealogy'", sliderInput("gens", "Number of Generations to Show", value=3, min=1, max=10))),
                     style="float:right;margin-right:10px;"                     
                   )
                 )
      )
    ),br(),
    
    # Chromosome Slider for CNV Search
    conditionalPanel(condition="input.plottype=='Search CNVs by Location'", 
                     uiOutput("ChrSlider")
                     ),
    
    # Variety List
    conditionalPanel(condition="input.plottype!='Genealogy' & 
                                input.plottype!='Phenotype Data' & 
                                input.plottype!='Search CNVs by Location' &
                                input.plottype!='Methodology'", 
                     uiOutput("VarietyList"), 
                     helpText("Type the name of the variety or click on the text box for a list of options")
                     ),  
    
    # Chromosome Lists
    conditionalPanel(condition="input.plottype!='Genealogy' & 
                                input.plottype!='Phenotype Data' & 
                                input.plottype!='Search CNVs by Location' &
                                input.plottype!='Methodology'",
                     uiOutput("ChromosomeList"), 
                     helpText("Type the chromosome number or click on the text box for a list of options")
                     ),
    
    # CNV Type Lists
    conditionalPanel(condition="input.plottype!='Genealogy' & 
                                input.plottype!='Phenotype Data' & 
                                input.plottype!='Methodology'",
                     uiOutput("CNVTypeList"), 
                     helpText("Select the types of CNV regions you would like to display")
    ),
    
    conditionalPanel(condition="input.plottype=='Search CNVs by Location'", 
                     uiOutput("ChrSelect2")
                     ),
      
    
    # Help text for Methodology
    conditionalPanel(condition="input.plottype=='Methodology'", 
                     helpText("This section details the data collection and analysis  process used to generate the results shown in this application.")
    ),

    # Help text for GlymaID list
    conditionalPanel(condition="input.plottype=='CNV List'", 
                     helpText("Use the checkboxes to filter the table output.", 
                              "Use the Search field on the right to filter by Glyma ID.")
                     ),
  
    # Animint plots
    conditionalPanel(condition="input.plottype=='Phenotype Data'", 
                     helpText("Click on a data point to the right to see genealogy 
                              and field trial data."), 
                     tagList(
                       tags$table(width='100%',
                         tags$tr(
                           tags$td(
                             tags$button(onClick="changeAll([\'yield\', 
                                                             \'protein\', 
                                                             \'oil\']);", 
                                         class="btn btn-info", 'Show/Hide')
                             )),
                         tags$tr(
                           tags$td('Yield, Protein, Oil by Year')
                           ),
                         tags$tr(tags$td(br())),
                         tags$tr(
                           tags$td(
                             tags$button(onClick="changeAll([\'maturity\', 
                                                             \'lodging\', 
                                                             \'seedsize\', 
                                                             \'seedquality\', 
                                                             \'plcount\']);", 
                                         class="btn btn-info", 'Show/Hide')
                           )),
                         tags$tr(
                           tags$td('More Field Trial Data by Year', br(), tags$small('Maturity, Lodging, Seeds'))
                         ),
                         tags$tr(),
                         tags$tr(tags$td(br())),
                         tags$tr(),
                         tags$tr(
                           tags$td(
                             tags$button(onClick="changeAll([\'yieldOil\', 
                                                             \'yieldProtein\', 
                                                             \'proteinOil\']);", 
                                         class="btn btn-info", 'Show/Hide')
                           )),
                       tags$tr(
                           tags$td('Yield, Protein, Oil Pairwise Plots')
                         )
                       )
                     )
                     ),

    # Variety List for Genealogy (separate from other Variety Lists)
    conditionalPanel(condition="input.plottype=='Genealogy'", 
                     helpText("Type the name of the variety or click on the text box for a list of options"), 
                     uiOutput("VarietyList2"))
  )),
  
  # Show a plot of the relevant chromosomes and varieties
  column(width=9,
    tabsetPanel(id="plottype",
      tabPanel("CNV Location", plotOutput("ChromosomePlot", width="100%"), helpText("Green lines indicate a significant CNV at that location")), 
      tabPanel("Copy Number", plotOutput("CopyNumberPlot", width="100%"), helpText("Open circles indicate a significant CNV at that location; darker blue lines indicate higher copy number, lighter lines indicate lower copy number.")),
      tabPanel("Search CNVs by Location", uiOutput("CNVList")),
      tabPanel("CNV List", uiOutput("GlymaTable")),
      tabPanel("Phenotype Data", uiOutput("YieldInfo")),
      tabPanel("Genealogy", plotOutput("FamilyTree", width="100%", height=600)),
      tabPanel("Methodology", includeHTML("Documentation.html"))
#       tabPanel("Methodology", uiOutput("Methods"))
      )
#     ,textOutput("DebugText")
  ))
))