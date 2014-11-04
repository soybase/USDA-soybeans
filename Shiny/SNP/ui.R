library(shiny)

head.scripts <- 
  tags$head(
#             tags$script(src="libs/jquery-1.11.0/jquery.min.js", type="text/javascript"),
            tags$link(href="libs/bootstrap3-3.2.0/css/bootstrap.min.css", rel="stylesheet"),
            tags$link(href="libs/bootstrap3-3.2.0/css/themes/cerulean/bootstrap.min.css", rel="stylesheet"),
#             tags$script(src="libs/bootstrap3-3.2.0/js/bootstrap.min.js", type="text/javascript"),
            tags$link(href="libs/highlightjs-8.2/highlight/idea.css", rel="stylesheet"),
            tags$script(src="libs/highlightjs-8.2/highlight.pack.js", type="text/javascript"),
            tags$link(href="libs/MagnificPopup-0.9.9/magnific-popup.css", rel="stylesheet"),
            tags$script(src="libs/MagnificPopup-0.9.9/magnific-popup.js", type="text/javascript"),
            tags$link(href="libs/knitrBootstrap-0.0.1/css/knitrBootstrap.css", rel="stylesheet" ),
            tags$script(src="libs/knitrBootstrap-0.0.1/js/knitrBootstrap.js", type="text/javascript"),
            tags$style(type="text/css", "input#bases { height:30pt;}"),
            tags$style(type="text/css", "input#glymaID { height:20pt;}"),
            tags$style(type="text/css", "input#glymaID2 { height:20pt;}"),
            tags$style(type="text/css", "input#chrStart { height:20pt;}"))

load("ShinyStart.rda")

AggSNPBrowser <- function(){
  tabPanel("Aggregate SNP Browser", 
           # Sidebar with checkbox inputs for varieties and chromosomes, 
           # plus a reset button.
           fluidRow(
             column(width=3, 
                    wellPanel(
                      textInput("glymaID", 
                                "Locate Position by Glyma ID",
                                value=""
                      ),
                      helpText("Ex: 01g004400 will search Chr 01 for IDs containing 004400"),
                      br(),
                      numericInput("bases",
                                   "# downstream SNPs (up to 50)", 
                                   value=20, min=5, max=50, step=5
                      )),
                    wellPanel(
                      helpText("Manually choose chromosome and location:"),
                      selectInput("locationChrs", "Choose Chromosome of Interest", 
                                  choices=unique(seqnames), multiple=FALSE, selectize=T),
                      helpText("Enter a numeric start point on the chromosome"),
                      textInput("chrStart", 
                                "Start point", value=0)
                      )
                    ),
             column(width=9, 
                    plotOutput("AggregatePlot", width="100%"),
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
                          tags$td(rowspan=6, width='70%', display='none', ''),
                          tags$td(
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
                      textInput("glymaID2", 
                                "Locate Position by Glyma ID",
                                value=""
                      ),
                      helpText("Ex: 01g004400 will search Chr 01 for IDs containing 004400"),
                      selectizeInput("varieties", "Choose up to 10 Cultivars of Interest", 
                                         choices=unique(varieties), multiple=TRUE, options=list(maxItems=10))
                      ),
                    h3("Matching Glyma IDs"),
                    dataTableOutput("glymaTable2")
           ),
           column(width=9, 
                  plotOutput("VarietySnpPlot", width="100%"),
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
                      selectizeInput("densityChrs", "Choose Chromosome(s) of Interest", choices=unique(seqnames), multiple=TRUE, options=list(maxItems=5)),
                      selectizeInput("densityVars", "Choose Varieties of Interest", choices=unique(varieties), multiple=TRUE)
                      )
                    ),
             column(width=10, 
                    plotOutput("DensityPlot", width="100%", height="800px"), 
                    helpText("Todo: Fill in sections with color key showing QTL category?")
                    )
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
                   AggSNPBrowser(),
                   VarSNPBrowser(),
                   SNPDensity(), 
                   SNPKinship(),
                   methodology(),
                   inverse=TRUE
                   )
        )
