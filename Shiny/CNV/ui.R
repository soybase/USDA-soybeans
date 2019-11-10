library(shiny)
library(dplyr)
library(magrittr)
library(stringr)
library(knitr)

# Load Chromosome names and varieties
load("ShinyStart.rda")
# Sort varieties and chromosomes
varieties <- sort(unique(varieties))
seqnames <- sort(chrs <- unique(seqnames))

cnvtypes <- c("gene", "CDS", "mRNA", "exon")

yield.btn <- tags$button(
  onClick="changeAll([\'yield\', \'protein\', \'oil\']);", 
  class="btn btn-info btn-xs", 'Show/Hide')

seedinfo.btn <- tags$button(
  onClick="changeAll([\'maturity\', \'lodging\', \'seedsize\', 
                      \'seedquality\', \'plcount\']);", 
  class="btn btn-info btn-xs", 'Show/Hide')

pairwise.btn <- tags$button(
  onClick="changeAll([\'yieldOil\', \'yieldProtein\', \'proteinOil\']);", 
  class="btn btn-info btn-xs", 'Show/Hide')

#------- Head Scripts ------
head.scripts <- 
  tags$head(
    singleton(tags$head(tags$script(src = "https://d3js.org/d3.v3.min.js"))),
    singleton(tags$script(src="google-analytics.js", type="text/javascript")),
    singleton(tags$style(type="text/css", "a[data-value='Copy Number Variation'] {font-size:18px; line-height:20px; height:50px;}")),
    singleton(tags$link(href="shiny.css", rel="stylesheet"))
  )
#--------------------------------------------------------------------------------

# Create variety list popover HTML
#-------------------------------------------------------------------------------
# HTML for list of varieties, to be placed in the "data-content" position in the popover definition
ncols <- 4
blanks <- ceiling(length(varieties)/ncols)*ncols-length(varieties)
varietyList <- 
  kable(
    # Create a matrix of varieties, with 4 varieties to a row, where each column is filled in order
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
  gsub(pattern="<thead>.*?</thead>", replacement = "") %>%
  # Replace '"' with '\'' to prevent breaking the string HTML
  str_replace_all(pattern="\"", replacement="\\\\'")
varietyListButton <- 
  # button definition which produces the popover when clicked
  tagList(
    tags$div(
      br(),
      tags$button(
        type = "button",
        class = "btn btn-info btn-xs", # this will create a blue button (btn-info) with our css
        "data-toggle" = "popover",
        container = ".well",
        style="width:100%",
        title = "Sequenced Varieties", # title of popover
        "data-placement" = "bottom",
        "data-content" = "insertHTMLhere",
        "data-trigger" = "click",
        "data-html" = TRUE,
        "data-viewport" = list(selector="body > div.container-fluid > div.row", padding=0),
        "Variety List" # button text
      ) %>%
        as.character %>%
        str_replace("insertHTMLhere", as.character(varietyList)) %>%
        HTML()
    )
  )
#--------------------------------------------------------------------------------

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
#--------------------------------------------------------------------------------

header <- {
  wellPanel(
    class="well-sm",
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
    ),
    fluidRow(
      # Overview Tab
      conditionalPanel(condition="input.tabname=='Copy Number Variation'",
                       div(helpText("Use the first two plots to select a chromosome and one or more varieties."), 
                       display="inline-table", align="center")),
      
      # Phenotype Tab
      conditionalPanel(condition="input.tabname=='Field Trials'",
                       column(3, helpText("Click on data points in the plot to see field trial data.")),
                       column(3, h5("Yield, Protein, and Oil by Year"), yield.btn),
                       column(3, h5("More Field Trial Data by Year"), 
                              tags$small('Maturity, Lodging, Seeds'), 
                              seedinfo.btn),
                       column(3, h5("Yield, Protein, Oil Pairwise Plots"), pairwise.btn)                     
      ),
      
      # Methodology Tab
      conditionalPanel(condition="input.tabname=='Methodology'",
                       div(helpText("This section details the data collection and analysis process used to generate the results shown in this application."), 
                           display="inline-table", align="center")),
      
      # Other Tabs
      conditionalPanel(
        condition="input.tabname!='Methodology' & 
                     input.tabname!='Copy Number Variation' & 
                     input.tabname!='Field Trials'",
        # Reset (and Download) Buttons
        column(3,
               conditionalPanel(
                 condition="input.tabname!='by Variety' & 
                              input.tabname!='by Location'",
                 tagList(
                   tags$table(width='100%',
                              style="border:0px solid;",
                              tags$tr(
                                style="border:0px solid;",
                                tags$td(style="text-align:center;border:0px solid;",
                                        actionButton("resetButton", "Clear Selections")
                                )
                              )
                   )
                 )
               ),
               conditionalPanel(
                 condition="!(input.tabname!='by Variety' & 
                                input.tabname!='by Location')",
                 tagList(
                   tags$table(width='100%',
                              style="border:0px solid;",
                              tags$tr(
                                style="border:0px solid;",
                                tags$td(style="text-align:left;border:0px solid;",
                                        actionButton("resetButton", "Clear Selections")),
                                tags$td(style="text-align:right;border:0px solid;",
                                        conditionalPanel(
                                          condition="input.tabname=='by Variety'",
                                          downloadButton("DataFrameDownload", 
                                                         label="Download", class=NULL)
                                        ),
                                        conditionalPanel(
                                          condition="input.tabname=='by Location'",
                                          downloadButton("DataFrameDownload2", 
                                                         label="Download", class=NULL)
                                        )
                                )
                              ),
                              tags$tr(
                                style="border:0px solid;",
                                tags$td(
                                  style="border:0px solid;",
                                  colspan=2, br(), 
                                  helpText("Use the Search field above the table to filter by Glyma ID."))
                              )
                   )
                 )
               )
        ),
        
        # Chromosome (ish) Information
        column(3,
               # Chromosome Checkboxes
               conditionalPanel(condition="input.tabname!='by Location' & 
                                             input.tabname!='Genealogy'", 
                                helpText("Type the chromosome number or click on the 
                                            text box for options"),
                                div(
                                  selectizeInput("chromosomes", "Choose Chromosomes", 
                                                 seqnames, NULL, multiple=TRUE), 
                                  display="inline-table", align="center")
               ),
               # Chromosome Radio Buttons
               conditionalPanel(condition="input.tabname=='by Location'", 
                                div(
                                  selectizeInput("locationChrs", "Choose Chromosome", 
                                                 seqnames, NULL, multiple=FALSE), 
                                  display="inline-table", align="center"
                                )
               ),
               # Generations Slider
               conditionalPanel(condition="input.tabname=='Genealogy'", 
                                div(
                                  sliderInput("gens", "Number of Generations to Show", 
                                              value=3, min=1, max=10)), 
                                display="inline-table", align="center")
        ),
        # Variety (ish) Information
        column(3, 
               # Variety Help Text
               conditionalPanel(condition="input.tabname!='by Location'", 
                                helpText("Type the name of the variety or click on the text box for options")),
               
              
               # Variety Selectize Input
               conditionalPanel(condition="input.tabname!='Genealogy' & 
                                             input.tabname!='by Location'", 
                                fluidRow(
                                  column(
                                    8, 
                                    div(
                                      selectizeInput("varieties", 
                                                     "Choose Varieties", 
                                                     varieties, NULL, multiple=TRUE), 
                                      display="inline-table", align="center")
                                  ), 
                                  column(
                                    4,
                                    varietyListButton,
                                    tags$script("$('body button[data-toggle=\"popover\"]').popover('toggle').popover('toggle');")
                                    )
                                )),
               
               # Variety Selectize Input (Genealogy Tab)
               conditionalPanel(condition="input.tabname=='Genealogy'",
                                fluidRow(
                                  column(
                                    8, 
                                    div(
                                      selectizeInput("genvarieties", "Choose Varieties", varieties, NULL, multiple=TRUE), 
                                      display="inline-table", align="center")
                                  ), 
                                  column(
                                    4,
                                    varietyListButton,
                                    tags$script("$('body button[data-toggle=\"popover\"]').popover('toggle').popover('toggle');")
                                  )
                                )), 
               
               # Chromosome Slider (Search CNVs by Location Tab)
               conditionalPanel(condition="input.tabname=='by Location'", 
                                div(sliderInput("chrRange", "Range to search for CNVs", 
                                                min=0, max=60000000, value=c(0, 600000000), 
                                                step=10000, round=FALSE), 
                                    display="inline-table", align="center"))
        ),
        column(3,
               conditionalPanel(condition="input.tabname!='Genealogy'",
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

overview <- function(){
  tabPanel("Copy Number Variation", 
           fluidRow(
             column(5, offset=1, helpText("The first plot ('Whole-Genome CNV Density') shows normalized CNV counts for all varieties over all chromosomes. Click on a chromosome region to view more detailed information for that chromosome.")),
             column(5, helpText("The second plot ('Number of CNVs by Variety') shows CNV counts for each variety, sorted from most CNVs to fewest CNVs. Click on one or more varieties to see more detailed information about which regions of the chromosome contain the most CNVs."))
           ),
           uiOutput("overview"),
           fluidRow(
             column(5, offset=1, helpText("The third plot ('Density of CNVs') shows the overall CNV distribution across all varieties in black, with selected varieties (from the top-right plot) shown in color.")),
             column(5, helpText("The last plot ('Difference in Normalized CNV Count') shows, for each selected variety, which regions of the chromosome contain relatively more CNVs than the average variety, and which regions of the chromosome contain relatively fewer CNVs."))
           ))
}
copyLocation <- function(){
  tabPanel("Location", 
           plotOutput("ChromosomePlot", width="100%"),
           br(), 
           helpText("Green lines indicate a significant CNV at that location"))
}
copyNumber <- function(){
  tabPanel("Copy Number", 
           plotOutput("CopyNumberPlot", width="100%"), 
           br(),
           helpText("Open circles indicate a significant CNV at that location; 
                       darker blue lines indicate higher copy number, 
                       lighter lines indicate lower copy number."))
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
navbarPage(
  title="", 
  overview(),
  navbarMenu("Barcode Plots",
             copyLocation(), 
             copyNumber()),
## 2017-07-03 (weeks): this hasn't worked since at least Feb. 2016
##                     (potential JavaScript issue); disabling for now
# navbarMenu("Search CNVs",
#   tabPanel("by Location", uiOutput("CNVList")),
#   tabPanel("by Variety", uiOutput("GlymaTable"))
#   ),
  tabPanel("Field Trials", uiOutput("YieldInfo")),
  tabPanel("Genealogy", plotOutput("FamilyTree", width="100%", height=600)),
  methodology(),
  singleton(tags$script("$('.dropdown-toggle').dropdown()", type="text/javascript")),
  header=header,
  id="tabname", 
  inverse=T,
  collapsible=F,
  windowTitle="Soybean CNVs"
)
