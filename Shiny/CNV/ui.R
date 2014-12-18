library(shiny)
# Load Chromosome names and varieties
load("ShinyStart.rda")
# Sort varieties and chromosomes
varieties <- sort(unique(varieties))
seqnames <- sort(chrs <- unique(seqnames))

cnvtypes <- c("gene", "CDS", "mRNA", "exon")

yield.btn <- tags$button(
  onClick="changeAll([\'yield\', \'protein\', \'oil\']);", 
  class="btn btn-info", 'Show/Hide')

seedinfo.btn <- tags$button(
  onClick="changeAll([\'maturity\', \'lodging\', \'seedsize\', 
                      \'seedquality\', \'plcount\']);", 
  class="btn btn-info", 'Show/Hide')

pairwise.btn <- tags$button(
  onClick="changeAll([\'yieldOil\', \'yieldProtein\', \'proteinOil\']);", 
  class="btn btn-info", 'Show/Hide')

header <- {
  wellPanel(
    class="well-sm",
    fluidRow(
      # Overview Tab
      conditionalPanel(condition="input.plottype=='Overview'",
                       div(helpText("Use the first two plots to select a chromosome and one or more varieties."), 
                       display="inline-table", align="center")),
      
      # Phenotype Tab
      conditionalPanel(condition="input.plottype=='Phenotype Data'",
                       column(3, helpText("Click on data points in the plot to see field trial data.")),
                       column(3, h6("Yield, Protein, and Oil by Year"), yield.btn),
                       column(3, h6("More Field Trial Data by Year"), 
                              tags$small('Maturity, Lodging, Seeds'), 
                              seedinfo.btn),
                       column(3, h6("Yield, Protein, Oil Pairwise Plots"), pairwise.btn)                     
      ),
      
      # Methodology Tab
      conditionalPanel(condition="input.plottype=='Methodology'",
                       div(helpText("This section details the data collection and analysis process used to generate the results shown in this application."), 
                           display="inline-table", align="center")),
      
      # Other Tabs
      conditionalPanel(
        condition="input.plottype!='Methodology' & 
                     input.plottype!='Overview' & 
                     input.plottype!='Phenotype Data'",
        # Reset (and Download) Buttons
        column(3,
               conditionalPanel(
                 condition="input.plottype!='CNV List' & 
                              input.plottype!='Search CNVs by Location'",
                 tagList(
                   tags$table(width='100%',
                              tags$tr(
                                tags$td(style="text-align:center;",
                                        actionButton("resetButton", "Clear Selections")
                                )
                              )
                   )
                 )
               ),
               conditionalPanel(
                 condition="!(input.plottype!='CNV List' & 
                                input.plottype!='Search CNVs by Location')",
                 tagList(
                   tags$table(width='100%',
                              tags$tr(
                                tags$td(style="text-align:left;",
                                        actionButton("resetButton", "Clear Selections")),
                                tags$td(style="text-align:right;",
                                        conditionalPanel(
                                          condition="input.plottype=='CNV List'",
                                          downloadButton("DataFrameDownload", 
                                                         label="Download", class=NULL)
                                        ),
                                        conditionalPanel(
                                          condition="input.plottype=='Search CNVs by Location'",
                                          downloadButton("DataFrameDownload2", 
                                                         label="Download", class=NULL)
                                        )
                                )
                              )
                   )
                 )
               )
        ),
        
        # Chromosome (ish) Information
        column(3,
               # Chromosome Checkboxes
               conditionalPanel(condition="input.plottype!='Search CNVs by Location' & 
                                             input.plottype!='Genealogy'", 
                                helpText("Type the chromosome number or click on the 
                                            text box for options"),
                                div(
                                  selectizeInput("chromosomes", "Choose Chromosomes", 
                                                 seqnames, NULL, multiple=TRUE), 
                                  display="inline-table", align="center")
               ),
               # Chromosome Radio Buttons
               conditionalPanel(condition="input.plottype=='Search CNVs by Location'", 
                                div(
                                  selectizeInput("locationChrs", "Choose Chromosome", 
                                                 seqnames, NULL, multiple=FALSE), 
                                  display="inline-table", align="center"
                                )
               ),
               # Generations Slider
               conditionalPanel(condition="input.plottype=='Genealogy'", 
                                div(
                                  sliderInput("gens", "Number of Generations to Show", 
                                              value=3, min=1, max=10)), 
                                display="inline-table", align="center")
        ),
        # Variety (ish) Information
        column(3, 
               # Variety Help Text
               conditionalPanel(condition="input.plottype!='Search CNVs by Location'", 
                                helpText("Type the name of the variety or click on the text box for options")),
               
               # Variety Selectize Input
               conditionalPanel(condition="input.plottype!='Genealogy' & 
                                             input.plottype!='Search CNVs by Location'", 
                                div(
                                  selectizeInput("varieties", "Choose Varieties", varieties, NULL, multiple=TRUE), 
                                  display="inline-table", align="center")),
               
               # Variety Selectize Input (Genealogy Tab)
               conditionalPanel(condition="input.plottype=='Genealogy'",  
                                div(
                                  selectizeInput("genvarieties", "Choose Varieties", varieties, NULL, multiple=TRUE), 
                                  display="inline-table", align="center")),
               
               # Chromosome Slider (Search CNVs by Location Tab)
               conditionalPanel(condition="input.plottype=='Search CNVs by Location'", 
                                div(uiOutput("ChrSlider"), 
                                    display="inline-table", align="center"))
        ),
        column(3,
               conditionalPanel(condition="input.plottype!='Genealogy'",
                                helpText("Select the types of CNV regions you would like to display"),
                                div(
                                  selectizeInput("featuretypes", "Choose Features (optional)", 
                                                 cnvtypes, cnvtypes, multiple=TRUE), 
                                  display="inline-table", align="center"))
        )
      )
    )
  )
}


# Define UI for page that allows selection of genetic lines with corresponding facets
navbarPage(
  title="Copy Number Variation", 
  tabPanel("Overview", uiOutput("overview")),
  tabPanel("CNV Location", 
           plotOutput("ChromosomePlot", width="100%"),
           br(), 
           helpText("Green lines indicate a significant CNV at that location")), 
  tabPanel("Copy Number", 
           plotOutput("CopyNumberPlot", width="100%"), 
           br(),
           helpText("Open circles indicate a significant CNV at that location; 
                       darker blue lines indicate higher copy number, 
                       lighter lines indicate lower copy number.")),
  tabPanel("Search CNVs by Location", uiOutput("CNVList")),
  tabPanel("CNV List", uiOutput("GlymaTable")),
  tabPanel("Phenotype Data", uiOutput("YieldInfo")),
  tabPanel("Genealogy", plotOutput("FamilyTree", width="100%", height=600)),
  tabPanel("Methodology", includeHTML("Documentation.html")),
#   tags$head(
#     tags$link(href="bootstrap/css/bootstrap.min.css", rel="stylesheet"),
#     tags$link(href="bootstrap/css/themes/cerulean/bootstrap.min.css", rel="stylesheet"),
#     tags$script(src="https://ajax.googleapis.com/ajax/libs/jquery/1.11.1/jquery.min.js"),
#     tags$script(src="bootstrap/js/bootstrap.min.js", type="text/javascript")
#   ,
  header=header,
  id="plottype", 
  inverse=T,
  collapsable=F,
  windowTitle="Soybean CNVs"
)
