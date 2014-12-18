library(shiny)
library(ggplot2)
library(RColorBrewer)
library(plyr)
library(reshape2)

# Data required for this app:
# created in generateDataObjects.R (run once to set up data for display)
# The creation of these objects is rather time-intensive on a multicore server
# so pre-processing makes a faster applet.
## res.df - Used to plot the cnv regions (not the copy numbers)
##          columns seqnames, start, end, Variety, CN, and width. 
##          Variety is used to match with segranges.df2, segments, segranges.df
##          seqnames contains strings Chr01, Chr02, etc.
##          res.df is the data frame output (using mold) of cnvs(res) 
##             where res is the cn.mops output object. 
## segments - Used to output the table of CNVs and glyma IDs 
##            columns seqnames, source, feature, start, end, score, strand, frame, ID
## segranges.df - Used to plot copy numbers
##                columns seqnames, start, end, width, strand, midpoint, Variety, count, CN
##                data frame output of mold(bamSegmentRanges) melted to have separate 
##                entries for each variety and gene. integer copy numbers are merged in 
##                from the command mold(integerCopyNumber(res)) where res is the output from 
##                cn.mops and then melted similarly to the segranges.df object.
## segranges.df2 - Used to plot read counts (the grey bars) in the chromosome plot
##                 columns seqnames, start, end, Variety, count.
##                 In order to reduce number of plotted objects, binned by 60K intervals
##                 (hence why it's a different object than segranges)

# Note: As of 1/21/14, may need to subtract 30 bp from the granges.start object, and add 30 bp from the granges.end object to account for input to cn.mops that should provide more stability.

load("tree.rda")
load("ChrPlot.rda")
load("GlymaIDs.rda")
load("ShinyStart.rda")
source("SelectGeneology.R")
glymacols <- c(11,1:4,17,8)


# Define server logic required to generate and plot a subset of varieties with a subset of chromosomes
shinyServer(function(input, output, session) {

  output$ChrSlider <- renderUI({
    if(length(input$locationChrs)>0){
      start <- min(floor(subset(chr.summary, seqnames%in%input$locationChrs)$start/10000)*10000)
      end <- max(ceiling(subset(chr.summary, seqnames%in%input$locationChrs)$end/10000)*10000)
    } else {
      start <- 0
      end <- 60000000
    }
    min <- 0
    max <- 60000000
    sliderInput("chrRange", "Range to search for CNVs", min=min, max=max, value=c(start, end), step=10000, round=FALSE)
  })
  
  # Reset input lists
  observe({
    if (input$resetButton == 0)
      return()
    
    isolate({
      # Reset each input using the appropriate update* function
      updateSelectizeInput(session = session, 
                               inputId = "varieties", 
                               label = "Choose Varieties", 
                               choices = sort(unique(varieties)), 
                               selected = NULL)
      updateSelectizeInput(session = session, 
                               inputId = "genvarieties", 
                               label = "Choose Varieties", 
                               choices = sort(unique(varieties)), 
                               selected = NULL)
      updateSelectizeInput(session = session, 
                               inputId = "chromosomes", 
                               label = "Choose Chromosomes", 
                               choices = sort(unique(seqnames)), 
                               selected = NULL)
      updateSelectizeInput(session = session, 
                           inputId = "featuretypes", 
                           label = "Choose Features (optional)", 
                           choices = c("gene", "CDS", "mRNA", "exon"), 
                           selected = c("gene", "CDS", "mRNA", "exon"))
      updateSliderInput(session = session, 
                        inputId = "chrRange", 
                        label = "Range to search for CNVs",
                        value = c(0, 60000000))
      updateRadioButtons(session = session, 
                         inputId = "locationChrs",
                         label = "Choose Chromosome of Interest",
                         choices = sort(unique(seqnames))
                         )
      
    })
  })
  
#   vals <- reactiveValues()
#   treedims <- reactiveValues()
  
  output$FamilyTree <- renderPlot({

    gen.vars2 <- input$genvarieties
    vals <- list()
    vals$gen.vars <- gen.vars2
    
    if(length(input$genvarieties)>0){
      temp2 <- ldply(gen.vars2, function(i){
        if(i %in% tree$child | i %in% tree$parent){
          temp <- rbind(node.to.data.frame(getancestors(i)), node.to.data.frame(getdescendants(i)))
          # Subset temp to get accurate number of generations according to the slider
          temp <- subset(temp, abs(gen)<=input$gens)
          # get plot coords of subsetted data frame
          temp <- cbind(variety=i, plotcoords(temp))
          temp$label2 <- temp$label
        } else {
          temp <- data.frame(variety=i, label=i, root=NA, root.gen=0, gen=0, type=0, branch=0, x=0, y=0, xstart=0, ystart=0, xend=0, yend=0, branchx=0, branchy=0, id=0, par.id=0, size=6, label2=paste(i, "-", "No Information Found", sep=""))
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
      
      temp <- merge(data.frame(variety=input$genvarieties,NewName=vals$gen.vars), ddply(vals$df, .(label), summarize, gen=mean(gen*c(-1,1)[(type=="descendant")+1])), by.x=2, by.y=1)
      vals$match <- temp[order(temp$gen, temp$variety),]
    } else {
      vals$df <- data.frame()
      vals$match <- data.frame()
      vals$gen.vars <- input$genvarieties
    }
    
    if(nrow(vals$df)>0){
      plot <- qplot(data=vals$df, x=x, y=y, label=label2, geom="text", vjust=-.25, hjust=.5, 
                    size=size, colour=color) +
                geom_segment(aes(x=xstart, y=ystart, xend=xend, yend=yend),inherit.aes=F) + 
                geom_segment(aes(x=xend, y=yend, xend=branchx, yend=branchy),inherit.aes=F) +
                facet_wrap(~variety, scales="free", ncol=2) +
                scale_size_continuous(range=c(3,6),guide="none") +
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
 
  # Expression that generates a plot of the chromosome. The expression
  # is wrapped in a call to renderPlot to indicate that:
  #
  #  1) It is "reactive" and therefore should be automatically 
  #     re-executed when inputs change
  #  2) Its output type is a plot 
  #
  output$ChromosomePlot <- renderPlot({

    if(length(input$chromosomes)>0 & length(input$varieties)>0){

      cnvs <- subset(glymaIDs, 
                       Chromosome%in%input$chromosomes & 
                       Variety%in%input$varieties & 
                       Feature%in%input$featuretypes)
      
      
      tmp <- data.frame(old.name=input$varieties, new.name=input$varieties)
      tmp$year <- sapply(tmp$new.name, function(i) min(tree$year[which(tree$child==i)]))
      tmp <- tmp[order(tmp$year),]
      
      data1 <- subset(segranges.df2, seqnames%in%input$chromosomes &
                        Variety%in%input$varieties)
      data1$Variety <- factor(data1$Variety, levels=tmp$old.name)
      
      
      plot <- ggplot() + 
        geom_rect(data=data1, 
                  aes(xmin=start, xmax=end, ymin=-.5, ymax=.5,
                      fill=count, colour=count)) +  
        scale_fill_gradient("# Reads", trans="log", 
                            low="#FFFFFF", high="#000000", 
                            breaks=c(1, 20, 400, 8000, 160000), 
                            na.value="#FFFFFF") + 
        scale_colour_gradient("# Reads", trans="log", 
                              low="#FFFFFF", high="#000000",
                              breaks=c(1, 20, 400, 8000, 160000),
                              na.value="#FFFFFF") + 
        facet_grid(Variety~seqnames, scales="free") + 
        theme_bw() + 
        theme(axis.text.y=element_blank(), 
              axis.ticks.y=element_blank(), 
              axis.title.y=element_blank())

      if(nrow(cnvs)>0){
        plot <- plot + 
          geom_rect(data=cnvs, 
                    aes(xmin=Cnv.Start, xmax=Cnv.End, ymin=-.6, ymax=.6),
                    fill="green", colour="green")
      }
    } else {
      plot <- ggplot() + 
        geom_text(aes(x=0, y=0, label="Please select varieties and chromosomes\n\n Note: It may take a minute to process the input")) +         
        theme_bw() + 
        theme(axis.text=element_blank(), 
              axis.ticks=element_blank(), 
              axis.title=element_blank())
    }
    
    print(plot)
  })
  
  output$CopyNumberPlot <- renderPlot({
#     set.df()
    pal = brewer.pal(9, "Blues")[-c(2:5)]
    
    tmp <- data.frame(old.name=input$varieties, new.name=input$varieties)
    tmp$year <- sapply(tmp$new.name, function(i) min(tree$year[which(tree$child==i)]))
    tmp <- tmp[order(tmp$year),]
    
    data1 <- subset(segranges.df, seqnames%in%input$chromosomes &
                      Variety%in%input$varieties)
    data1$Variety <- factor(data1$Variety, levels=tmp$old.name)
    
    if(length(input$chromosomes)>0 & length(input$varieties)>0){
      plot <- ggplot() + 
        geom_rect(data=subset(chr.summary, seqnames%in%input$chromosomes),
                  aes(xmin=start, xmax=end, ymin=-.5, ymax=.5, fill=CN, colour=CN)) +
        geom_rect(data=data1, 
                  aes(xmin=start, xmax=end, ymin=-.5, ymax=.5, fill=CN, colour=CN)) +
        scale_fill_gradientn(colours=pal,na.value = "grey50") +
        scale_colour_gradientn(colours=pal, na.value = "grey50") +
        facet_grid(Variety~seqnames, scales="free") + 
        theme_bw() + 
        theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.title.y=element_blank()) +
        xlab("")
      
      cnvs <- subset(glymaIDs, Chromosome%in%input$chromosomes & Variety%in%input$varieties & Feature%in%input$featuretypes)
      if(nrow(cnvs)>0){
        plot <-  plot + geom_point(data=subset(cnvs), aes(x=(Cnv.Start+Cnv.End)/2, y=-.55, shape="CNV")) + 
          scale_shape("CNV Location", solid=FALSE)
      }
      
    } else {
      plot <- ggplot() + 
        geom_text(aes(x=0, y=0, label="Please select varieties and chromosomes\n\n Note: It may take a minute to process the input")) +         
        theme_bw() + 
        theme(axis.text=element_blank(), 
              axis.ticks=element_blank(), 
              axis.title=element_blank())
    }
    
    print(plot)
  })
    
  output$glymaIDs <- renderDataTable({

#     
#     if(length(input$chromosomes)==0 & length(input$varieties)>0){
#       temp <- unique(subset(glymaIDs, Variety %in% input$varieties)[,glymacols])
#     } else if(length(input$chromosomes)>0 & length(input$varieties)==0){
#       temp <- unique(subset(glymaIDs, Chromosome %in% input$chromosomes)[,glymacols])
#     } else if(length(input$chromosomes)>0 & length(input$varieties)>0){
#       temp <- unique(subset(glymaIDs, Chromosome%in%input$chromosomes & Variety%in%input$varieties)[,glymacols])
#     } else { # neither varieties nor chromosomes selected
#       temp <- unique(subset(glymaIDs)[,glymacols])
#     }
    
    if(length(input$chromosomes)>0){
      temp <- unique(subset(glymaIDs, Chromosome %in% input$chromosomes)[,glymacols])
    } else {
      temp <- unique(glymaIDs)[,glymacols]
    }
    
    if(length(input$varieties)>0){
      temp <- unique(subset(temp, Variety%in%input$varieties))
    }
    
    if(length(input$featuretypes)>0){
      temp <- unique(subset(temp, Feature%in%input$featuretypes))
    }
    
    if(nrow(temp)>0){
      temp
    }else{
      data.frame()
    }
  })

  output$glymaIDs2 <- renderDataTable({
    if(length(input$locationChrs)==0 & length(input$chrRange)<=1){
      temp <- unique(subset(glymaIDs)[,glymacols])
    } else if(length(input$locationChrs)>0 & length(input$chrRange)<=1){
      temp <- unique(subset(glymaIDs, Chromosome %in% input$locationChrs)[,glymacols])
    } else if(length(input$locationChrs)>0 & length(input$chrRange)==2){
      temp <- unique(subset(glymaIDs, Chromosome %in% input$locationChrs & Cnv.Start>=input$chrRange[1] & Cnv.End <= input$chrRange[2])[,glymacols])
    } else { # neither start/end nor chromosomes selected
      temp <- unique(subset(glymaIDs)[,glymacols])
    }
    
    if(length(input$featuretypes)>0){
      temp <- unique(subset(temp, Feature%in%input$featuretypes))
    }
    
    if(nrow(temp)>0){
      temp
    }else{
      data.frame()
    }
  })

  output$DataFrameDownload <- downloadHandler(  
    filename=function(){
      chrs <- ifelse(length(input$chromosomes)>0, paste0("Chr", paste(gsub("Chr", "-", input$chromosomes), collapse="")), "AllChr")
      paste("CNVList-", chrs, '.csv', sep='')
      },
    content=function(con){    
      if(length(input$chromosomes)==0 & length(input$varieties)>0){
        temp <- unique(subset(glymaIDs, Variety %in% input$varieties)[,glymacols])
      } else if(length(input$chromosomes)>0 & length(input$varieties)==0){
        temp <- unique(subset(glymaIDs, Chromosome %in% input$chromosomes)[,glymacols])
      } else if(length(input$chromosomes)>0 & length(input$varieties)>0){
        temp <- unique(subset(glymaIDs, Chromosome%in%input$chromosomes & 
                                Variety%in%input$varieties)[,glymacols])
      } else { # neither varieties nor chromosomes selected
        temp <- unique(subset(glymaIDs)[,glymacols])
      }
      if(nrow(temp)>0) write.csv(temp, con) else write.csv(data.frame(), con)
    }
    )

  output$DataFrameDownload2 <- downloadHandler(  
    filename=function(){
      if(length(input$chrRange)>1){
        paste(input$locationChrs, "start=", input$chrRange[1], end=input$chrRange[2], '.csv', sep='')
      } else {
        paste(input$locationChrs, '.csv', sep='')
      }
    },
    content=function(con){      
      if(length(input$locationChrs)==0 & length(input$chrRange)<=1){
        temp <- unique(subset(glymaIDs)[,glymacols])
      } else if(length(input$locationChrs)>0 & length(input$chrRange)<=1){
        temp <- unique(subset(glymaIDs, Chromosome %in% input$locationChrs)[,glymacols])
      } else if(length(input$locationChrs)>0 & length(input$chrRange)==2){
        temp <- unique(subset(glymaIDs, Chromosome%in%input$locationChrs & Cnv.Start>=input$chrRange[1] & Cnv.End <= input$chrRange[2])[,glymacols])
      } else { # neither start/end nor chromosomes selected
        temp <- unique(subset(glymaIDs)[,glymacols])
      }
      
      if(nrow(temp)>0) write.csv(temp, con) else write.csv(data.frame(), con)
    }
  )

  output$GlymaTable <- renderUI({
    tagList(
      dataTableOutput("glymaIDs"),
      br(),
      helpText("Use the Search field (top-right) to filter by Glyma ID.")
      )
  })

  output$CNVList <- renderUI({
    dataTableOutput("glymaIDs2")
  })

  output$Methods <- renderUI({
    tagList(h3("Data Collection"), 
            p("Biology info goes here... plant methods, climate conditions, whatever other green/slimy stuff..."),
            h3("Processing"), 
            p("Andrew's stuff goes here... raw data to BAM file pipeline with parameters"),
            h3("Copy Number Identification"),
            h5("cn.mops"),
            p("Using BAM files from the previous step as input to ", a(href="http://www.bioinf.jku.at/software/cnmops/", target="_blank", "cn.mops."), "The program was executed separately on each geneomic feature (gene, exon, CDS, mRNA) to provide internal verification as well as reduce the problem to a more computationally manageable size. As suggested in the cn.mops manual, each region of the genome was extended by 30 bp on each side to aid in identification of CNV regions. Once these modifications were performed, the command used to run the analysis was:", tags$blockquote("exomecn.mops(bamSegmentRanges, parallel=32, I=c(0.025, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5), classes=paste0('CN', as.character(0:10)))"), "The results from the algorithm were back-transformed (regions were reduced by 30 bp on each side) and merged with annotation files."),
            h5("Additional Analysis"),
            p("In order to verify the reliability of cn.mops, cross-validation was performed on a subset of lines. This step ensures that the results from the algorithm are not unduly affected by one or two varieties. Results of cross validation confirm that most of the variants identified by cn.mops are present in 90% or more of the cross-validation runs."),
            h3("Genealogy"),
            p("Link sources for genealogy compilation here")
    )
  })

  output$YieldInfo <- renderUI({
    tagList(
      singleton(tags$head(tags$script(src = "fieldtrials/animint.js", type='text/javascript'))),
      singleton(tags$head(tags$script(src = "fieldtrials/vendor/d3.v3.js", type='text/javascript'))),
      singleton(tags$head(tags$link(rel = "stylesheet", type = "text/css", href="fieldtrials/styles.css"))),
      singleton(tags$body(tags$script("
                      changeHiding = function(id){
                        var plotEl = document.getElementById(id).parentNode.parentNode.parentNode;
                        if(plotEl.style.display=='inline-block'){
                          return(plotEl.style.display='none');
                        } else {
                          return(plotEl.style.display='inline-block');
                        }
                      }
                      changeAll = function(vec){
                        vec.forEach(function(d){ changeHiding(d)});
                      }
                                      "))),
      tags$div(id="plot"),
      tags$script("var plot = new animint('#plot','fieldtrials/plot.json');"),
      tags$script("var monitor = true; 
                   var doMonitor = function(){
                      monitor=document.getElementById('yieldOil')==null | document.getElementById('maturity')==null;
                      if(!monitor){
                        nonyear = document.getElementById('yieldOil').parentNode.parentNode.parentNode; 
                        extrainfo = document.getElementById('maturity').parentNode.parentNode.parentNode;
                        var br = document.createElement('br'); 
                        document.getElementById('plot').insertBefore(br, nonyear);
                        var br = document.createElement('br'); 
                        document.getElementById('plot').insertBefore(br, extrainfo);
                      } else {
                        setTimeout(doMonitor, 500);
                      }
                   }
                   setTimeout(doMonitor,500);
                   ")
    )
  })

  output$overview <- renderUI({
    tagList(
      singleton(tags$head(tags$script(src = "overview/animint.js", type='text/javascript'))),
      singleton(tags$head(tags$script(src = "overview/vendor/d3.v3.js", type='text/javascript'))),
      singleton(tags$head(tags$link(rel = "stylesheet", type = "text/css", href="overview/styles.css"))),
      tags$div(id="plot", align="center"),
      tags$script("var plot = new animint('#plot','overview/plot.json');"),
      tags$script("var monitor = true; 
                  var doMonitor = function(){
                  monitor=document.getElementById('chromosomeplot')==null | document.getElementById('overview')==null;
                  if(!monitor){
                  tmp = document.getElementById('chromosomeplot').parentNode.parentNode.parentNode; 
                  var br = document.createElement('br'); 
                  document.getElementById('plot').insertBefore(br, tmp);
                  } else {
                  setTimeout(doMonitor, 500);
                  }
                  }
                  setTimeout(doMonitor,500);
                  ")
      )
  })

})
