library(shiny)
library(knitr) # required to create a simple table of output for popover variety list
library(stringr)
library(dplyr)

head.scripts <- 
  tags$head(
    singleton(tags$link(href="shiny.css", rel="stylesheet")),
    singleton(tags$link(href="libs/bootstrap-3.3.4/css/themes/cerulean/bootstrap.min.css", rel="stylesheet")),
    # singleton(tags$link(href="libs/bootstrap-3.3.4/css/bootstrap.min.css", rel="stylesheet")),
    # singleton(tags$script(src="libs/bootstrap-3.3.4/js/bootstrap.min.js", type="text/javascript")),
    
    # Needed for display of the methodology section tab
    # singleton(includeHTML("www/knitrBootstrapIncludes.html")),
    # singleton(tags$script(src="libs/toc.js", type="text/javascript")),
    singleton(tags$script(src="libs/refs.js", type="text/javascript")),
    singleton(tags$link(href="libs/knitrBootstrap.css", type="text/javascript")),
    
    singleton(tags$link(href="libs/highlightjs-8.2/highlight/idea.css", rel="stylesheet")),
    singleton(tags$script(src="libs/highlightjs-8.2/highlight.pack.js", type="text/javascript")),
    singleton(tags$link(href="libs/MagnificPopup-0.9.9/magnific-popup.css", rel="stylesheet")),
    singleton(tags$script(src="libs/MagnificPopup-0.9.9/magnific-popup.js", type="text/javascript")),
    singleton(tags$link(href="libs/knitrBootstrap-0.0.1/css/knitrBootstrap.css", rel="stylesheet" )),
    singleton(tags$script(src="libs/knitrBootstrap-0.0.1/js/knitrBootstrap.js", type="text/javascript"))
  )

# Load list of varieties and chromosomes
load("ShinyStart.rda")


# Create variety list popover HTML
#-------------------------------------------------------------------------------
# HTML for list of varieties, to be placed in the "data-content" position in the popover definition
ncols <- 4
blanks <- ceiling(length(varieties)/ncols)*ncols-length(varieties)
varietyList <- 
  kable(
    # Create a matrix of varieties, with 3 varieties to a row, where each column is filled in order
    # Append blanks to varieties list so that the vector fits into this matrix with no missing values
    matrix(
      c(varieties, rep("", blanks)), 
      ncol=ncols, byrow=F
    ), 
    format="html",
    row.names=F,
    col.names=rep("", ncols),
    align="l",
    padding=2) %>% 
  # Replace the table header of blank column names with an empty string
  str_replace(pattern="<thead>.*?</thead>", "") %>%
  # Replace '"' with '\'' to prevent breaking the string HTML
  str_replace_all(pattern="\"", replacement="\\\\'")
varietyListButton <- 
  # button definition which produces the popover when clicked
  tagList(
    tags$div(
      br(),
      tags$button(
        type = "button",
        class = "btn btn-info btn-small", # this will create a blue button (btn-info) with our css
        "data-toggle" = "popover",
        title = "Sequenced Cultivars", # title of popover
        "data-placement" = "bottom",
        "data-content" = "insertHTMLhere",
        "data-trigger" = "click",
        "data-html" = TRUE,
        "data-viewport" = list(selector="body > div.container-fluid > div.row", padding=0),
        "List of Cultivars" # button text
      ) %>%
        as.character %>%
        str_replace("insertHTMLhere", as.character(varietyList)) %>%
        HTML()
    )
  )
  
#-------------------------------------------------------------------------------


headerDef <- function(){
  tagList(
    conditionalPanel(
      # Load head.scripts after tabs have been loaded to prevent file not found errors
      condition="!input.tabname==''",
      head.scripts,    
      tags$script(
        'Shiny.addCustomMessageHandler(
                      \'setTab\',
                      function(data) {
                        var nav_ref = \'li a:contains(\\"\' + data + \'\\")\';
                        $(nav_ref).tab(\'show\');
                      });
         Shiny.addCustomMessageHandler(
                      \'javascript\',
                      function(data) {
                        eval(data);
                      });'
      )
    ),
    conditionalPanel(
      # If using reactivity (i.e. not relatedness and SNPs or methodology) include
      # a header for inputs
      condition="!(input.tabname=='Methodology' | input.tabname=='Relatedness and SNPs')",
      fluidRow(
        column(
          width=5, offset=1,
          # Left wellpanel: 
          # Position information (chromosome or GlymaID)
          wellPanel(
            fluidRow(
              conditionalPanel(
                condition="input.tabname=='Aggregated SNPs' | input.tabname=='Cultivar-Level SNPs' | input.tabname=='Inheritance of SNPs'",
                column(
                  width=4,
                  h4("GlymaID navigation"),
                  helpText("Ex: 01g004700 will search Chr 01 for IDs containing 004700")
                ),
                column(
                  width=4,
                  textInput("glymaID", "Locate Position by Glyma ID", value="01g004700")
                )
              ),
              conditionalPanel(
                condition="input.tabname=='Aggregated SNPs'",
                column(
                  width=4,
                  numericInput("bases", "# downstream SNPs (up to 50)", 
                               value=20, min=5, max=50, step=5)
                )
              ),
              conditionalPanel(
                condition="input.tabname=='Cultivar-Level SNPs'",
                column(
                  width=4,
                  helpText("The table on the left will show matching GlymaIDs for partial entries.")
                )
              ),
              conditionalPanel(
                condition="input.tabname=='find a GlymaID'",
                column(
                  width=6,
                  helpText("Select one or more chromosomes. The first table will show the number of SNP sites for each glymaID on the chromosome.")
                )
              ),
              conditionalPanel(
                condition="input.tabname=='Overview: SNP Locations'",
                column(width=3,
                       helpText("Select one or more chromosomes to see the distribution of SNPs."))
              ),
              conditionalPanel(
                condition="input.tabname=='find a GlymaID' | input.tabname=='Overview: SNP Locations'",
                column(
                  width=4,
                  selectizeInput("glymaChrs", "Show Chromosome(s)", 
                                 choices=unique(seqnames), selected="Chr01", multiple=TRUE, options=list(maxItems=5))
                )
              )
            )
          )
        ),
        column(
          width=5,
          # Right wellPanel
          # Variety information 
          # or manual navigation for aggregate browser
          wellPanel(
            fluidRow(
              conditionalPanel(
                condition="input.tabname=='Aggregated SNPs'",
                column(
                  width=3,
                  h4("Manual navigation"),
                  helpText("Choose Chromosome and position")
                ),
                column(
                  width=3,
                    selectInput("locationChrs", "Choose Chromosome",  
                                choices=unique(seqnames), selected="Chr01", 
                                multiple=FALSE, selectize=T)
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
                condition="input.tabname=='find a GlymaID'",
                column(
                  width=2,
                  helpText("Input a glyma ID such as 01g004700")
                ),
                column(
                  width=4,
                  textInput("glymaID3", "View SNP sites within a GlymaID", value="01g004700")
                ),
                column(
                  width=6,
                  helpText("For a selected glymaID, the second table shows each SNP location, and the number of varieties with SNPs at that position.")
                )
              ),
              conditionalPanel(
                condition="input.tabname=='Cultivar-Level SNPs' | input.tabname=='Overview: SNP Locations'",
                column(
                  width=4,
                  offset=2,
                  selectizeInput("varieties", "Cultivars of Interest (up to 10)", 
                                 choices=unique(varieties), multiple=TRUE,
                                 options=list(maxItems=10))
                )
              ),
              conditionalPanel(
                condition="input.tabname=='Cultivar-Level SNPs' | input.tabname=='Overview: SNP Locations'",
                column(
                  width=6,
                  varietyListButton,
                  tags$script("$('body button[data-toggle=\"popover\"').popover('toggle').popover('toggle')")
                )
              ),
              conditionalPanel(
                condition="input.tabname=='Inheritance of SNPs'",
                column(
                  width=3,
                  offset=1,
                  selectizeInput("variety0", "Cultivar of Interest", 
                                 choices=unique(varieties), multiple=FALSE, 
                                 selected="A.K.")
                ),
                column(
                  width=1,
                  varietyListButton,
                  tags$script("$('body button[data-toggle=\"popover\"').popover('toggle').popover('toggle')")
                ),
                column(
                  offset=1,
                  width=3,
                  checkboxInput("ancestors", "Search Ancestors?", value=TRUE),
                  checkboxInput("descendants", "Search Descendants?", value=TRUE)
                ),
                column(
                  width=2,
                  numericInput("gens", "# Generations", value=2, min=1, max=4)
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
  tabPanel("Aggregated SNPs", 
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
                        style="border:0px solid;background-color:#ffffff;",
                        tags$tr(
                          style="border:0px solid;background-color:#ffffff;",
                          tags$td(
                            style="border:0px solid;background-color:#ffffff;",
                            actionButton("up", 
                                         "Move upstream", 
                                         icon=icon('backward', class="glyphicon", lib="glyphicon"),
                                         class="btn-info"
                            )
                          ),
                          tags$td(rowspan=7, width='75%', display='none',
                                  style="border:0px solid;background-color:#ffffff;", ''),
                          tags$td(align="right",
                                  style="border:0px solid;background-color:#ffffff;",
                                  actionButton("down", 
                                               span("Move downstream", 
                                               icon('forward', class="glyphicon", lib="glyphicon")),
                                               class="btn-info"
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
  tabPanel("Cultivar-Level SNPs", 
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

GenealogySNPBrowser <- function(){
  tabPanel("Inheritance of SNPs", 
           # Sidebar with checkbox inputs for varieties and chromosomes, 
           # plus a reset button.
           fluidRow(
             column(width=3, 
                    h3("Genealogy"),
                    plotOutput("GenealogyTree", width="100%", height="600px"),
                    helpText("SNP data is available for cultivars shown in black.")
             ),
             column(width=9, 
                    plotOutput("GenealogySnpPlot", width="100%", height="800px")
             )
           )
  )
}

SNPDensity <- function(){
  tabPanel(
    "Overview: SNP Locations", 
    fluidRow(plotOutput("DensityPlot", width="100%", height="800px"))
  )
}

SNPSummary <- function(){
  tabPanel("find a GlymaID", 
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
  tabPanel("Relatedness and SNPs", fluidRow(uiOutput("KinshipSNP")))
}

methodology <- function(){
  text <- includeHTML("Documentation.html")
  header <- str_extract(text, "<head>(.*?)</head>")
  header <- str_replace_all(header, "</?head>", "") %>%
    str_replace_all(fixed("libs/bootstrap3-3.2.0"), "libs/bootstrap-3.3.4")
  text <- str_replace(text, "<head>(.*?)</head>", "")
  
  tabPanel("Methodology", 
           fluidRow(
             HTML(text)
           )
  )
}

# Define UI for page that allows selection of genetic lines with corresponding facets
shinyUI(
  navbarPage(
    title="Soybean SNPs", 
    header=headerDef(),
    SNPDensity(), 
    navbarMenu(
      "Browse SNPs Visually",
      AggSNPBrowser(),
      VarSNPBrowser(),
      GenealogySNPBrowser()
    ),
    navbarMenu(
      "Locate SNPs",
      SNPSummary()
    ),
    SNPKinship(),
    methodology(),
    inverse=TRUE,
    id="tabname"
  )
)
