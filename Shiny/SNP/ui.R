library(shiny)

head.scripts <- 
  tags$head(tags$link(href="libs/bootstrap3-3.2.0/css/bootstrap.min.css", rel="stylesheet"),
            tags$link(href="libs/bootstrap3-3.2.0/css/themes/cerulean/bootstrap.min.css", rel="stylesheet"),
            tags$link(href="libs/highlightjs-8.2/highlight/idea.css", rel="stylesheet"),
            tags$script(src="libs/highlightjs-8.2/highlight.pack.js", type="text/javascript"),
            tags$link(href="libs/MagnificPopup-0.9.9/magnific-popup.css", rel="stylesheet"),
            tags$script(src="libs/MagnificPopup-0.9.9/magnific-popup.js", type="text/javascript"),
            tags$link(href="libs/knitrBootstrap-0.0.1/css/knitrBootstrap.css", rel="stylesheet" ),
            tags$script(src="libs/knitrBootstrap-0.0.1/js/knitrBootstrap.js", type="text/javascript"),
            tags$style(type="text/css", "input#bases { height:30pt;}"),
            tags$style(type="text/css", "input#glymaID { height:30pt;}"),
            tags$style(type="text/css", "input#glymaID2 { height:30pt;}"),
            tags$style(type="text/css", "input#glymaID3 { height:30pt;}"),
            tags$style(type="text/css", "input#chrStart { height:30pt;}"))

load("ShinyStart.rda")

AggSNPBrowser <- function(){
  tabPanel("Aggregate SNP Browser", 
           # Sidebar with checkbox inputs for varieties and chromosomes, 
           # plus a reset button.
           fluidRow(
             column(width=3, 
                    wellPanel(
                      textInput("glymaID", #done
                                "Locate Position by Glyma ID",
                                value=""
                      ),
                      helpText("Ex: 01g004400 will search Chr 01 for IDs containing 004400"),
                      br(),
                      numericInput("bases",  # done
                                   "# downstream SNPs (up to 50)", 
                                   value=20, min=5, max=50, step=5
                      )),
                    wellPanel(
                      helpText("Manually choose chromosome and location:"),
                      selectInput("locationChrs", "Choose Chromosome of Interest",  # done
                                  choices=unique(seqnames), multiple=FALSE, selectize=T),
                      helpText("Enter a numeric start point on the chromosome"),
                      textInput("chrStart",  # done
                                "Start point", value=0)
                      )
                    ),
             column(width=9, 
                    plotOutput("AggregatePlot", width="100%", height="400px"),
                    tagList(
                      tags$table(
                        tags$tr(
                          tags$td(
                            actionButton("up", 
                                         "Move upstream", 
                                         icon=icon('angle-double-left', 
                                                   class="2x pull-left")
                            )
                          ),
                          tags$td(rowspan=7, width='75%', display='none', ''),
                          tags$td(align="right",
                            actionButton("down", 
                                         "Move downstream", 
                                         icon=icon('angle-double-right', 
                                                   class="2x pull-right")
                                         )
                            )
                          )
                        )
                      ),
                    br(),
                    h3("Matching Glyma IDs"),
                    dataTableOutput("glymaTable")
                    )
             )
  )
}

VarSNPBrowser <- function(){
  tabPanel("Variety-Level SNP Browser", 
           # Sidebar with checkbox inputs for varieties and chromosomes, 
           # plus a reset button.
           fluidRow(
             column(width=3, 
                    wellPanel(
                      textInput("glymaID2",  # done
                                "Locate Position by Glyma ID",
                                value=""
                      ),
                      helpText("Ex: 01g004400 will search Chr 01 for IDs containing 004400"),
                      selectizeInput("varieties", "Choose up to 10 Cultivars of Interest", # done
                                         choices=unique(varieties), multiple=TRUE, options=list(maxItems=10))
                      ),
                    h3("Matching Glyma IDs"),
                    dataTableOutput("glymaTable2")
           ),
           column(width=9, 
                  plotOutput("VarietySnpPlot", width="100%", height="800px"),
                  br(),
                  h3("Matching SNPs"),
                  dataTableOutput("snpTable")
           )
           )
           )
}

SNPDensity <- function(){
  tabPanel("SNP Density", 
           fluidRow(
             column(width=2, 
                    wellPanel(
                      selectizeInput("densityChrs", "Choose Chromosome(s) of Interest", 
                                     choices=unique(seqnames), multiple=TRUE, 
                                     options=list(maxItems=5)), # done
                      selectizeInput("densityVars", "Choose Varieties of Interest", 
                                     choices=unique(varieties), multiple=TRUE) # done
                      )
                    ),
             column(width=10, 
                    plotOutput("DensityPlot", width="100%", height="800px"), 
                    helpText("Todo: Fill in sections with color key showing QTL category?")
                    )
           )
  )
}

SNPSummary <- function(){
  tabPanel("SNP Counts by GlymaID", 
           wellPanel(
             fluidRow(
               column(width=2),
               column(width=3, 
                      selectizeInput("glymaChrs", "Filter GlymaIDs by Chromosome(s)", # done
                                     choices=unique(seqnames), multiple=TRUE, options=list(maxItems=5)),
                      helpText("The table below will show GlymaIDs with SNPs")
                      ), 
               column(width=2),
               column(width=3,
                      textInput("glymaID3", # done
                                "View SNP sites within a GlymaID",
                                value=""),
                      helpText("Input a glyma ID such as 01g000200")
                      ),
               column(width=2)
             )
           ),
           fluidRow(
             column(width=6, wellPanel(
               helpText("Use the table below to see what GlymaIDs have SNPs on a specific chromosome."),
               dataTableOutput("glymaSummary")
               )),
             column(width=6, wellPanel(
               helpText("Use the table below to see SNPs at specific locations within a GlymaID."),
               dataTableOutput("positionSummary")
             ))
           )
  )
}

SNPKinship <- function(){
  tabPanel("Kinship via SNPs", fluidRow(uiOutput("KinshipSNP")))
}

methodology <- function(){
  tabPanel("Methodology", 
           fluidRow(
             includeHTML("Documentation.html"),
             tagList(head.scripts)
             )
           )
}

# Define UI for page that allows selection of genetic lines with corresponding facets
shinyUI(navbarPage(title="Soybean SNPs", 
                   header=tags$script(
                     'Shiny.addCustomMessageHandler(
                      \'setTab\',
                      function(data) {
                        var nav_ref = \'li a:contains(\\"\' + data + \'\\")\';
                        $(nav_ref).tab(\'show\');
                      });'
                   ),
                   AggSNPBrowser(),
                   SNPSummary(),
                   VarSNPBrowser(),
                   SNPDensity(), 
                   SNPKinship(),
                   methodology(),
                   inverse=TRUE
                   )
        )
