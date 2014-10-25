library(shiny)
library(ggplot2)
library(RColorBrewer)
library(plyr)
library(reshape2)


# Data required for this app:
# created in generateDataObjects.R (run once to set up data for display)

load("tree.rda")
source("SelectGeneology.R")
library(phyViz)
treeGraph <- processTreeGraph(tree)

objlist <- ls()
if(!"kevinbacon"%in%objlist){
  load("kevinbacon.rda")
}

fix <- function(x){
  tolower(gsub("\\.", "", gsub(" ", "", gsub("-", "", x))))
}

# Define server logic required to generate and plot a subset of varieties with a subset of chromosomes
shinyServer(function(input, output, session) {

  output$FamilyTree <- renderPlot({
    objlist <- ls()
    if(!"tree"%in%objlist){
      load("tree.rda")
    }
    gen.vars2 <- input$variety
    vals <- list()
    vals$gen.vars <- gen.vars2
    
    if(length(input$variety)>0){
      temp2 <- ldply(gen.vars2, function(i){
        if(i %in% tree$child | i %in% tree$parent){
          temp <- rbind(node.to.data.frame(getancestors(i, maxgen=input$gens)), node.to.data.frame(getdescendants(i, maxgen=input$gens)))
          # Subset temp to get accurate number of generations according to the slider
#           temp <- subset(temp, abs(gen)<=input$gens)
          # get plot coords of subsetted data frame
          temp <- cbind(variety=i, plotcoords(temp))
          temp$label2 <- temp$label
        } else {
          temp <- data.frame(variety=i, label=i, root=NA, root.gen=0, gen=0, type=0, branch=0, x=0, y=0, xstart=0, ystart=0, xend=0, yend=0, branchx=0, branchy=0, id=0, par.id=0, size=6, label2=paste(i, "-", "Query not found in database", sep=""))
        }
        unique(temp[,-which(names(temp)%in%c("id", "par.id"))])
      })
      
      cols <- hcl(h=seq(0, 300, by=50), c=80, l=55, fixup=TRUE)
      temp2$color <- "#000000"
      if(length(gen.vars2)<=7){
        var.list <- which(temp2$label%in%gen.vars2)
        temp2$color[var.list] <- cols[as.numeric(factor(temp2$label[var.list]))]
      }
      vals$df <- temp2
      
      temp <- merge(data.frame(variety=input$variety,NewName=vals$gen.vars), ddply(vals$df, .(label), summarize, gen=mean(gen*c(-1,1)[(type=="descendant")+1])), by.x=2, by.y=1)
      vals$match <- temp[order(temp$gen, temp$variety),]
    } else {
      vals$df <- data.frame()
      vals$match <- data.frame()
      vals$gen.vars <- input$variety
    }
    
    if(nrow(vals$df)>0){
      plot <- ggplot(data=vals$df) + 
        geom_text(aes(x=x, y=y, label=label2, colour=color, size=size), vjust=-.25, hjust=.5) +
        geom_segment(data=vals$df, aes(x=xstart, xend=xend, y=ystart, yend=yend),inherit.aes=F) + 
        geom_segment(data=vals$df, aes(x=xend, xend=branchx, y=yend, yend=branchy),inherit.aes=F) +
        facet_wrap(~variety, scales="free", ncol=2) +
        scale_size_continuous(range=c(4,6),guide="none") +
        scale_colour_identity() +
        theme_bw() +
        theme(axis.title=element_blank(), 
              axis.text=element_blank(), 
              axis.ticks=element_blank()) + 
        scale_x_continuous(expand = c(.1, 1.075)) + 
        scale_y_continuous(expand = c(.1, 1.075))
    } else {
      plot <- ggplot() + 
        geom_text(aes(x=0, y=0, label="Please select varieties\n\n Note: It may take a minute to process the input")) +         
        theme_bw() + 
        theme(axis.text=element_blank(), 
              axis.ticks=element_blank(), 
              axis.title=element_blank())
    }
    print(plot)
  })


  output$KevinBaconDistance <- renderPlot({  
    path <- getPath(input$var1, input$var2, ig=treeGraph)
    print(generatePathPlot(path))
  })
  
#   output$nameList <- renderDataTable({
#     if(length(input$varieties)>0){
#       temp <- data.frame(variety=unique(tree$child), fixvar=fix(unique(tree$child)))
#       temp <- subset(temp, !grepl(" x ", variety) & grepl(fix(input$variety), fixvar, fixed=TRUE))
#       
#       temp[,1]
#     } else {
#       data.frame(Error="No Lines Selected", Hint="Select a line in the first text box")
#     }
#    
#   })

  
  output$Methods <- renderUI({
    tagList()
  })
  
})