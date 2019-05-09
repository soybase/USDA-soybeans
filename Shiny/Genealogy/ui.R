library(shiny)

# format checkbox inputs in 3 columns with some padding.
format.checkboxlist <- tags$head(tags$style(type="text/css", 
                                            "label.radio { display: inline-table; width: 25%; margin:2%;} 
                                            label.checkbox { display: inline-table; width: 45%; margin:2%;} 
                                            input[type='text']{ padding: 2px; height:30px;} 
                                            input[type='number'] { padding: 2px; height:35px;}")) 
format.selectize <- tags$head(tags$style(type="text/css",
                                         ".selectize-control {
                                            width: 80%;
                                            line-height: 14px; 
                                            }"))

# need only 1st two columns (child and parent)
tree <- read.csv("soybeanGenealogy.csv.gz", colClasses=c("character","character", rep(NULL,7)))
varieties <- unique(c(tree$child, tree$parent))
varieties <- varieties[!grepl(" x ", varieties)]

TreeBrowser <- function(){
  tabPanel("Genealogy Browser", 
           fluidRow(
             column(6,
                    format.selectize,
                    selectizeInput("variety", "Lines to Query", choices=varieties, multiple=TRUE, options=list(maxItems=6))),
             column(6,
                    sliderInput("gens", "# Generations to show: ", min=1, max=6, value = 2)
                    )
             ),
           fluidRow(
             column(12, 
                    plotOutput("FamilyTree", width='100%', height="600px")
                    )
             )
           )
  
}

TraitLinkage <- function(){
  tabPanel("Phenotype Similarity", NULL
           )
}

KevinBacon <- function(){
  tabPanel("Generational Distance", 
           fluidRow(
             column(4, 
                    format.selectize,
                    selectInput("var1", "Variety 1", choices=varieties, selected="Tokyo", multiple=FALSE, selectize=T)),
             column(4, 
                    format.selectize,
                    selectInput("var2", "Variety 2", choices=varieties, selected="Amsoy", multiple=FALSE, selectize=T)),
             column(4, 
                    helpText("Select two varieties to see if/how they are related"))
             ),
           fluidRow(
             column(8, plotOutput("KevinBaconDistance", width='100%', height="600px")),
             column(4, plotOutput("Var12", width='100%', height="600px"))
           )
  )
}

Phenotype <- function(){
  tabPanel("Phenotype Data Browser", NULL
           )
}

methodology <- function(){
  tabPanel("Methodology", fluidRow(includeHTML("Documentation.html")))
}

#---- Page UI ----#
shinyUI(navbarPage(title="Soybean Genealogy", 
                   TreeBrowser(),
                   KevinBacon(),
#                    TraitLinkage(),
#                    KevinBacon(), 
#                    Phenotype(),
                   methodology(),
                   inverse=TRUE
)
)
