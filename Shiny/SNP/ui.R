library(shiny)
library(knitr) # required to create a simple table of output for popover variety list
library(stringr)
library(dplyr)

# Shiny unexported function
`%AND%` <- function (x, y) {
  if (!is.null(x) && !is.na(x)) 
    if (!is.null(y) && !is.na(y)) 
      return(y)
  return(NULL)
}

#------- Head Scripts ------

head.scripts <- 
  tags$head(
#     singleton(tags$script(type="text/javascript", charset="utf8", src="https://code.jquery.com/ui/1.10.3/jquery-ui.min.js")),
#     singleton(tags$script(type="text/javascript", charset="utf8", src="http://cdn.datatables.net/1.10.3/js/jquery.dataTables.js")),
#     singleton(tags$link(rel="stylesheet", type="text/css", href="http://cdn.datatables.net/1.10.3/css/jquery.dataTables.css")),
#     singleton(tags$link(href="libs/bootstrap-3.3.4/css/themes/cerulean/bootstrap.min.css", rel="stylesheet")),
#     singleton(tags$link(href="libs/bootstrap-3.3.4/css/bootstrap.min.css", rel="stylesheet")),
    singleton(tags$script("$('.dropdown-toggle').dropdown()", type="text/javascript")),
    singleton(tags$script(src="libs/bootstrap-3.3.4/js/bootstrap.min.js", type="text/javascript")),
    
    singleton(tags$link(href="shiny.css", rel="stylesheet"))
    
    # Needed for display of the methodology section tab
    # singleton(includeHTML("www/knitrBootstrapIncludes.html"))#,
    # singleton(tags$script(src="libs/toc.js", type="text/javascript")),
    # singleton(tags$script(src="libs/refs.js", type="text/javascript"))#,
    # singleton(tags$link(href="libs/knitrBootstrap.css", type="text/javascript")),
    
    # singleton(tags$link(href="libs/highlightjs-8.2/highlight/idea.css", rel="stylesheet")),
    # singleton(tags$script(src="libs/highlightjs-8.2/highlight.pack.js", type="text/javascript")),
    # singleton(tags$link(href="libs/MagnificPopup-0.9.9/magnific-popup.css", rel="stylesheet")),
    # singleton(tags$script(src="libs/MagnificPopup-0.9.9/magnific-popup.js", type="text/javascript")),
    # singleton(tags$link(href="libs/knitrBootstrap-0.0.1/css/knitrBootstrap.css", rel="stylesheet" )),
    # singleton(tags$script(src="libs/knitrBootstrap-0.0.1/js/knitrBootstrap.js", type="text/javascript"))
  )

#--------- Load Data ------

# Load list of varieties and chromosomes
load("ShinyStart.rda")
#-----

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

helpButton <- function(label="?", popTitle = "Title", text="Help text", style=""){
  # button definition which produces a help popover when clicked
  tagList(
    tags$button(
      type = "button",
      class = "btn btn-xs", # this will create a button 
      "data-toggle" = "popover",
      title = popTitle, # title of popover
      "data-placement" = "right",
      "data-content" = text,
      "data-trigger" = "click",
      "data-html" = TRUE,
      "data-viewport" = list(selector="body .container-fluid .row .well", padding=0),
      style = style,
      label # button text
    ) 
  )
}

#-------------------------------------------------------------------------------


headerDef <- function(){
  tagList(
#-------
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
      ) # End custom message handler
    ), # End conditional panel for script loading
#-------
    conditionalPanel(
      # If using reactivity (i.e. not relatedness and SNPs or methodology) include
      # a header for inputs
      condition="!(input.tabname=='Methodology' | input.tabname=='Relatedness and SNPs')",
      fluidRow(
        style="margin-left:0px;margin-right:0px",
        column(
          #-------
          width=6, 
          # Left wellpanel: 
          # Position information (chromosome or GlymaID)
          wellPanel(
            fluidRow(
              style="margin-left:0px;margin-right:0px",
              conditionalPanel(
                condition="input.tabname=='Aggregated SNPs' | input.tabname=='Cultivar-Level SNPs' | input.tabname=='Inheritance of SNPs'",
                column(
                  width=4,
                  style="padding:0px;",
                  h4("GlymaID navigation"),
                  helpText("Ex: 01g000900 searches Chr 01 for ID 000900")
                ),
                column(
                  width=4,
                  div(class = "form-group shiny-input-container", 
                      "Navigate by GlymaID" %AND% 
                        tags$label("Navigate by GlymaID", `for` = "glymaID"), 
                      helpButton(
                        label = "?", 
                        popTitle = "Navigate by GlymaID", 
                        text = "<p>Enter a GlymaID to navigate to that region of a specific chromosome. </p> The table on the left will show matching GlymaIDs for partial entries.",
                        style="float:right"),
                      div(
                        tags$input(
                          id = "glymaID", type = "text", 
                          class = "form-control", value = "01g000900", style="width:100%;"),
                        style="padding-right:.5em;overflow:hidden;"
                      )
                  )
                )
              ),
              conditionalPanel(
                condition="input.tabname=='Aggregated SNPs'",
                column(
                  width=3,
                  numericInput("bases", "# SNPs shown", 
                               value=20, min=5, max=50, step=5)
                )
              ),
              conditionalPanel(
                condition="input.tabname=='Find a GlymaID'",
                column(
                  width=6,
                  helpText("Select one or more chromosomes. The first table will show the number of SNP sites for each glymaID on the chromosome.")
                )
              ),
              conditionalPanel(
                condition="input.tabname=='Overview: SNP Locations'",
                column(width=6,
                       helpText("Select one or more chromosomes to see the distribution of SNPs."))
              ),
              conditionalPanel(
                condition="input.tabname=='Find a GlymaID' | input.tabname=='Overview: SNP Locations'",
                column(
                  width=6,
                  selectizeInput("glymaChrs", "Show Chromosome(s)", 
                                 choices=unique(seqnames), selected="Chr01", multiple=TRUE, options=list(maxItems=5))
                )
              )
            )
          )
        ), # End first column definition
        column(
          width=6,
          # Right wellPanel
          # Variety information 
          # or manual navigation for aggregate browser
          wellPanel(
            fluidRow(
              style="margin-left:0px;margin-right:0px",
              conditionalPanel(
                condition="input.tabname=='Aggregated SNPs'",
                column(
                  width=4,
                  style="padding:0px;",
                  h4("Manual navigation"),
                  helpText("Choose chromosome and position")
                ),
                column(
                  width=4,
                    selectInput("locationChrs", "Choose Chromosome",  
                                choices=unique(seqnames), selected="Chr01", 
                                multiple=FALSE, selectize=T)
                ),
                column(
                  width=4,
                  div(class = "form-group shiny-input-container", 
                      "Start point" %AND% 
                        span(
                          tags$label("Start point", `for` = "chrStart"), 
                          helpButton(
                            label="?", 
                            popTitle = "Chromosome Start Point", 
                            text="Enter a numeric start point on the Chromosome")
                        ),
                        tags$input(id = "chrStart", type = "text", 
                                   class = "form-control", value = 0)
                  )
                )
              ),
              conditionalPanel(
                condition="input.tabname=='Find a GlymaID'",
                column(
                  width=6,
                  div(class = "form-group shiny-input-container", 
                      "View SNPs within a GlymaID" %AND% 
                        span(
                          tags$label("View SNPs within a GlymaID", `for` = "glymaID3"), 
                          helpButton(
                            label = "?", 
                            popTitle = "Help", 
                            text = "Input a glyma ID such as 01g000900")
                        ),
                      tags$input(id = "glymaID3", type = "text", 
                                 class = "form-control", value = "01g000900")
                  )
                ),
                column(
                  width=6,
                  helpText("For a selected glymaID, the second table shows each SNP location, and the number of varieties with SNPs at that position.")
                )
              ),
              conditionalPanel(
                condition="input.tabname=='Cultivar-Level SNPs' | input.tabname=='Overview: SNP Locations'",
                column(
                  width=8,
                  selectizeInput("varieties", "Cultivars of Interest (up to 10)", 
                                 choices=unique(varieties), multiple=TRUE,
                                 options=list(maxItems=10))
                )
              ),
              conditionalPanel(
                condition="input.tabname=='Cultivar-Level SNPs' | input.tabname=='Overview: SNP Locations'",
                column(
                  width=4,
                  varietyListButton,
                  tags$script("$('body button[data-toggle=\"popover\"]').popover('toggle').popover('toggle');")
                )
              ),
              conditionalPanel(
                condition="input.tabname=='Inheritance of SNPs'",
                column(
                  width=3,
                  selectizeInput("variety0", "Show Cultivar", 
                                 choices=unique(varieties), multiple=FALSE, 
                                 selected="A3127")
                ),
                column(
                  width=3,
                  varietyListButton,
                  tags$script("$('body button[data-toggle=\"popover\"]').popover('toggle').popover('toggle');")
                ),
                column(
                  width=4,
                  checkboxInput("ancestors", "Show Ancestors", value=TRUE),
                  checkboxInput("descendants", "Show Descendants", value=TRUE)
                ),
                column(
                  width=2,
                  numericInput("gens", "Generations", value=2, min=1, max=4)
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
                    plotOutput("AggregatePlot", width="100%", height="500px"),
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
                    plotOutput("VarietySnpPlot", height="500px"),
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
                    plotOutput("GenealogyTree", width="100%", height="450px"),
                    helpText("Data is available for cultivars shown in black.")
             ),
             column(width=9, 
                    plotOutput("GenealogySnpPlot", width="100%", height="550px")
             )
           )
  )
}

SNPDensity <- function(){
  tabPanel(
    "Overview: SNP Locations", 
    fluidRow(plotOutput("DensityPlot", width="100%", height="550px"))
  )
}

SNPSummary <- function(){
  tabPanel("Find a GlymaID", 
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
    SNPSummary(),
    navbarMenu(
      "Browse SNPs Visually",
      AggSNPBrowser(),
      VarSNPBrowser(),
      GenealogySNPBrowser()
    ),
    SNPKinship(),
    methodology(),
    singleton(tags$script("$('.dropdown-toggle').dropdown()", type="text/javascript")),
    inverse=TRUE,
    id="tabname"
  )
)
