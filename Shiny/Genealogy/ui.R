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


load("varieties.rda")

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
                    selectizeInput("var1", "Variety 1", choices=varieties, multiple=FALSE)),
             column(4, 
                    format.selectize,
                    selectizeInput("var2", "Variety 2", choices=varieties, multiple=FALSE)),
             column(4, 
                    helpText("Select two varieties to see if/how they are related"))
             ),
           column(12, 
                  plotOutput("KevinBaconDistance", width='100%', height="600px")
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
