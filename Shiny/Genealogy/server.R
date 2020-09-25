library(shiny)
library(ggplot2)
library(plyr)
library(reshape2)


# Data required for this app:
# created in generateDataObjects.R (run once to set up data for display)

source("SelectGeneology.R") # defines tree
library(phyViz)
treeGraph <- treeToIG(tree, isDirected=FALSE)

fix <- function(x){
  tolower(gsub("\\.", "", gsub(" ", "", gsub("-", "", x))))
}

plotFamilyTree <- function(variety, gens, ncol=2){
  gen.vars2 <- variety
  vals <- list()
  vals$gen.vars <- gen.vars2
  
  if(length(variety)>0){
    temp2 <- ldply(gen.vars2, function(i){
      if(i %in% tree$child | i %in% tree$parent){
        temp <- rbind(node.to.data.frame(getancestors(i, maxgen=gens)), node.to.data.frame(getdescendants(i, maxgen=gens)))
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
    
    temp <- merge(data.frame(variety=variety,NewName=vals$gen.vars), ddply(vals$df, .(label), summarize, gen=mean(gen*c(-1,1)[(type=="descendant")+1])), by.x=2, by.y=1)
    vals$match <- temp[order(temp$gen, temp$variety),]
  } else {
    vals$df <- data.frame()
    vals$match <- data.frame()
    vals$gen.vars <- variety
  }
  
  if(nrow(vals$df)>0){
    plot <- ggplot(data=vals$df) + 
      geom_text(aes(x=x, y=y, label=label2, colour=color, size=size), vjust=-.25, hjust=.5) +
      geom_segment(data=vals$df, aes(x=xstart, xend=xend, y=ystart, yend=yend),inherit.aes=F) + 
      geom_segment(data=vals$df, aes(x=xend, xend=branchx, y=yend, yend=branchy),inherit.aes=F) +
      facet_wrap(~variety, scales="free", ncol=ncol) +
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
  return(plot)
}

# Define server logic required to generate and plot a subset of varieties with a subset of chromosomes
shinyServer(function(input, output, session) {

  output$FamilyTree <- renderPlot({
    print(plotFamilyTree(variety=input$variety, gens=input$gens))
  })

  output$Var12 <- renderPlot({
    print(plotFamilyTree(variety=c(input$var1, input$var2), gens=1, ncol=1))
  })

  output$KevinBaconDistance <- renderPlot({  
    if(input$var1==input$var2){
      tmp <- data.frame(facet="Path Between Varieties", label="The two varieties are the same", x=0, y=0)
      plotdf <- qplot(data=tmp, label=label, x=x, y=y, geom="text") + 
        theme(axis.text = element_blank(), axis.ticks = element_blank(), 
              axis.title = element_blank(), legend.position = "none", 
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background=element_blank()) + 
        facet_wrap(~facet)
    } else {
      path <- getPath(input$var1, input$var2, ig=treeGraph, tree=tree)
      if(length(path)!=2){
        tmp <- data.frame(facet="Path Between Varieties", label="No path between the two varieties", x=0, y=0)
        plotdf <- qplot(data=tmp, label=label, x=x, y=y, geom="text") + 
          theme(axis.text = element_blank(), axis.ticks = element_blank(), 
                axis.title = element_blank(), legend.position = "none", 
                panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.background=element_blank()) + 
          facet_wrap(~facet)
      } else {
        dframe <- phyViz:::buildPathDF(path)
        dframe$facet <- "Path Between Varieties"
        dframe <- transform(dframe, 
                            w = strwidth(label, "inches") + 0.25, 
                            h = strheight(label, "inches") + 0.25)
        textFrame <- data.frame(x = dframe$x, y = dframe$y, label = dframe$label)
        textFrame <- transform(textFrame, 
                               w = strwidth(label, "inches") + 0.25, 
                               h = strheight(label, "inches") + 0.25)
        textFrame$facet <- "Path Between Varieties"
        plotdf <- ggplot(data = dframe, aes(x = x, y = y)) + 
                    geom_rect(data = textFrame, aes(xmin = x - w * 1.5 - 1, ymin = y - 0.1,
                                                    xmax = x + w * 1.5 + 1, ymax = y + 0.1), 
                              fill = "grey80") + 
                    geom_segment(aes(x = x - w * 1.5 - 1, 
                                     y = y - 0.1, 
                                     xend = x + w * 1.5 + 1, 
                                     yend = y - 0.1)) + 
                    geom_segment(aes(x = x - w * 1.5 - 1, 
                                     y = y + 0.1, 
                                     xend = x + w * 1.5 + 1, 
                                     yend = y + 0.1)) + 
                    geom_segment(aes(x = xstart, y = ystart, 
                                     xend = xend, yend = yend)) + 
                    geom_text(data = textFrame, aes(x = x, y = y, label = label), size = 4) + 
                    xlab("Year") + 
                    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), 
                          axis.title.y = element_blank(), legend.position = "none", 
                          panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
                          panel.background=element_blank(),
                          panel.border=element_rect(colour="grey80", fill=NA)) +
                    facet_wrap(~facet)
      }
    }
    print(plotdf)
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
