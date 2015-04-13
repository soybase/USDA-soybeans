library(shiny)

head.scripts <- 
  tags$head(
    tags$link(href="libs/bootstrap3-3.2.0/css/bootstrap.min.css", rel="stylesheet"),
    tags$link(href="libs/bootstrap3-3.2.0/css/themes/cerulean/bootstrap.min.css", rel="stylesheet"),
    tags$link(href="libs/highlightjs-8.2/highlight/idea.css", rel="stylesheet"),
    tags$script(src="libs/highlightjs-8.2/highlight.pack.js", type="text/javascript"),
    tags$link(href="libs/MagnificPopup-0.9.9/magnific-popup.css", rel="stylesheet"),
    tags$script(src="libs/MagnificPopup-0.9.9/magnific-popup.js", type="text/javascript"),
    tags$link(href="libs/knitrBootstrap-0.0.1/css/knitrBootstrap.css", rel="stylesheet" ),
    tags$script(src="libs/knitrBootstrap-0.0.1/js/knitrBootstrap.js", type="text/javascript"),
    tags$link(href="shiny.css", rel="stylesheet")
  )

load("ShinyStart.rda")

headerDef <- function(){
  tagList(
    head.scripts,
    tags$script(
      'Shiny.addCustomMessageHandler(
                      \'setTab\',
                      function(data) {
                        var nav_ref = \'li a:contains(\\"\' + data + \'\\")\';
                        $(nav_ref).tab(\'show\');
                      });'
    ),
    conditionalPanel(
      condition="!(input.tabname=='Methodology' | input.tabname=='Kinship via SNPs')",
      fluidRow(
        column(
          width=5, offset=1,
          wellPanel(
            fluidRow(
              conditionalPanel(
                condition="input.tabname=='Aggregate SNP Browser' | input.tabname=='Variety-Level SNP Browser'",
                column(
                  width=4,
                  h4("GlymaID navigation"),
                  helpText("Ex: 01g004700 will search Chr 01 for IDs containing 004700")
                ),
                column(
                  width=4,
                  textInput("glymaID", "Locate Position by Glyma ID", value="")
                )
              ),
              conditionalPanel(
                condition="input.tabname=='Aggregate SNP Browser'",
                column(
                  width=4,
                  numericInput("bases", "# downstream SNPs (up to 50)", 
                               value=20, min=5, max=50, step=5)
                )
              ),
              conditionalPanel(
                condition="input.tabname=='Variety-Level SNP Browser'",
                column(
                  width=4,
                  helpText("The table on the left will show matching GlymaIDs for partial entries.")
                )
              ),
              conditionalPanel(
                condition="input.tabname=='SNP Counts by GlymaID'",
                column(
                  width=6,
                  helpText("Select one or more chromosomes. The first table will show the number of SNP sites for each glymaID on the chromosome.")
                )
              ),
              conditionalPanel(
                condition="input.tabname=='SNP Density'",
                column(width=3,
                       helpText("Select one or more chromosomes to see the distribution of SNPs."))
              ),
              conditionalPanel(
                condition="input.tabname=='SNP Counts by GlymaID' | input.tabname=='SNP Density'",
                column(
                  width=4,
                  selectizeInput("glymaChrs", "Show Chromosome(s)", 
                                 choices=unique(seqnames), multiple=TRUE, options=list(maxItems=5))
                )
              )
            )
          )
        ),
        column(
          width=5,
          wellPanel(
            fluidRow(
              conditionalPanel(
                condition="input.tabname=='Aggregate SNP Browser'",
                column(
                  width=3,
                  h4("Manual navigation"),
                  helpText("Choose Chromosome and position")
                ),
                column(
                  width=3,
                    selectInput("locationChrs", "Choose Chromosome",  
                                choices=unique(seqnames), multiple=FALSE, selectize=T)
                ),
                column(
                  width=2,
                  helpText("Enter a numeric start point")
                ),
                column(
                  width=4,
                  textInput("chrStart", "Start point", value=0)
                )
              ),
              conditionalPanel(
                condition="input.tabname=='SNP Counts by GlymaID'",
                column(
                  width=2,
                  helpText("Input a glyma ID such as 01g004700")
                ),
                column(
                  width=4,
                  textInput("glymaID3", "View SNP sites within a GlymaID", value="")
                ),
                column(
                  width=6,
                  helpText("For a selected glymaID, the second table shows each SNP location, and the number of varieties with SNPs at that position.")
                )
              ),
              conditionalPanel(
                condition="input.tabname=='Variety-Level SNP Browser' | input.tabname=='SNP Density'",
                column(
                  width=4,
                  offset=4,
                  selectizeInput("varieties", "Cultivars of Interest (up to 10)", 
                                 choices=unique(varieties), multiple=TRUE, options=list(maxItems=10))
                )
              )
            )
          )
        )
      )
    )
  )
}

AggSNPBrowser <- function(){
  tabPanel("Aggregate SNP Browser", 
           # Sidebar with checkbox inputs for varieties and chromosomes, 
           # plus a reset button.
           fluidRow(
             column(width=3, 
                    wellPanel(
                      h3("Matching Glyma IDs"),
                      dataTableOutput("glymaTable")
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
                    )
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
  tabPanel(
    "SNP Density", 
    fluidRow(plotOutput("DensityPlot", width="100%", height="800px"))
  )
}

SNPSummary <- function(){
  tabPanel("SNP Counts by GlymaID", 
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
shinyUI(
  navbarPage(
    title="Soybean SNPs", 
    header=headerDef(),
    AggSNPBrowser(),
    SNPSummary(),
    VarSNPBrowser(),
    SNPDensity(), 
    SNPKinship(),
    methodology(),
    inverse=TRUE,
    id="tabname"
  )
)
