library(shiny)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(reshape2)
library(stringr)

load("ShinyStart.rda")
pal <- brewer.pal(8, "Paired")[c(1, 2, 7, 8)]

objlist <- ls()
if(!"snp.counts"%in%objlist){
  load("SNPCounts.RData")
  snp.counts <- snp.counts %>% arrange(Chromosome, Position) %>% group_by(Chromosome) 
}

if(!"GlymaIDList"%in%objlist){
  load("GlymaID.rda")
  GlymaIDList <- GlymaIDList %>% group_by(chrnum, numid)
}

if(!"snpList"%in%objlist){
  load("snpList.rda")
  snpList <- snpList %>% group_by(Chromosome, Variety)
}

if(!"snp.density"%in%objlist){
#   snp.density <- tbl(src_sqlite("SNPdb.sqlite"), "snpDensity") %>% group_by(Chromosome, Variety)
  load("SNPDensity.RData")
}

# Data required for this app:
## snp.df - data frame containing information from the info() and geno() commands applied
##          to the VCF file. 

# Define server logic required to generate and plot a subset of varieties with a subset of chromosomes
shinyServer(function(input, output, session) {

  
  # Return SNPs for selected varieties
  varsnps <- reactive({
    id2.cur <- id2()
    tmp <- filter(id2.cur, shown)
    
    chr <- paste0("Chr", substr(tmp$ID, 7, 8))
    position.min=min(tmp$start)
    position.max=max(tmp$end)
    tmplist <- as.data.frame(
                 filter(snpList, Chromosome == chr) %>%
                 filter(Position >= position.min & 
                        Position <= position.max))
    if(nrow(tmplist)>0 & length(input$varieties)>0){
      tmplist <- filter(tmplist, Variety%in%input$varieties)
    }
    tmplist
    })
  
  snps <- reactive(filter(snp.counts, 
                          Chromosome%in%c(input$locationChrs, "")) %>%
                   filter(Position>=min(as.numeric(input$chrStart), Inf)))
  
  snpDensity <- reactive(filter(snp.density, 
                                Chromosome%in%c(input$densityChrs, "")) %>%
                         filter(Variety%in%c(input$densityVars, "")))
  
  # search for glyma IDs
  id <- reactive({
    if(nchar(input$glymaID)>0){
      # strip off extra glyma stuff to get 01g000000
      tmp <- gsub("wm82a2v1", "", gsub("glyma", "", gsub(".", "", tolower(input$glymaID), fixed=TRUE)))
      # attempt to match strings nicely
      chr.str <- word(tmp, sep="[Gg]")
      pos.str <- word(tmp, sep="[Gg]", start=-1)
      if(nchar(chr.str)==2 & nchar(pos.str)<=6){
        row.idx <- chr.str==tolower(as.data.frame(select(GlymaIDList, chrnum))$chrnum) & grepl(pos.str, as.data.frame(select(GlymaIDList, numid))$numid)
        res <- as.data.frame(GlymaIDList)[row.idx,]
      } else{
        row.idx <- grepl(tmp, as.data.frame(select(GlymaIDList, searchstr))$searchstr)
        res <- as.data.frame(GlymaIDList)[row.idx,]
      }
      if(nrow(res)>0){
        res$shown <- FALSE
        res$shown[unique(res$numid)[1]==res$numid] <- TRUE
      }      
    } else {
      res <- data.frame(shown=NULL)
    }
    
    res
  })
  
  # search for glyma IDs (variety level)
  id2 <- reactive({
    if(nchar(input$glymaID2)>0){
      # strip off extra glyma stuff to get 01g000000
      tmp <- gsub("wm82a2v1", "", gsub("glyma", "", gsub(".", "", tolower(input$glymaID2), fixed=TRUE)))
      # attempt to match strings nicely
      chr.str <- word(tmp, sep="g")
      pos.str <- word(tmp, sep="g", start=-1)
      if(nchar(chr.str)==2 & nchar(pos.str)<=6){
        row.idx <- chr.str==tolower(as.data.frame(select(GlymaIDList, chrnum))$chrnum) & grepl(pos.str, as.data.frame(select(GlymaIDList, numid))$numid)
        res <- as.data.frame(GlymaIDList)[row.idx,]
      } else{
        row.idx <- grepl(tmp, as.data.frame(select(GlymaIDList, searchstr))$searchstr)
        res <- as.data.frame(GlymaIDList)[row.idx,]
      }
      if(nrow(res)>0){
        res$shown <- FALSE
        res$shown[unique(res$numid)[1]==res$numid] <- TRUE
      }    
    } else {
      res <- data.frame(shown=NULL)
    }
    res
  })
  
  # Update start and chromosome with glyma ID search
  observe({
    x <- id()
    # Check to see if there is more than one ID matching the search string; 
    # if so, update the start location and chromosome with corresponding values, otherwise, 
    # update the start location with the minimum matching location and most common chromosome. 
    if(nrow(x)==1){
      updateNumericInput(session, "chrStart", value=x$start)
      updateRadioButtons(session, "locationChrs", choices=seqnames, selected=x$seqnames)
    } else if(nrow(x)>1){
      z <- table(x$seqnames)
      updateNumericInput(session, "chrStart", value=min(x$start))
      updateRadioButtons(session, "locationChrs", choices=seqnames, selected=names(z)[which.max(z)])
    }
  })
#   
#   # Update start and chromosome with glyma ID search
#   observe({
#     x <- id2()
#     # Check to see if there is more than one ID matching the search string; 
#     # if so, update the start location and chromosome with corresponding values, otherwise, 
#     # update the start location with the minimum matching location and most common chromosome. 
#     if(nrow(x)==1){
#       updateNumericInput(session, "chrStart2", value=x$start)
#       updateRadioButtons(session, "locationChrs2", choices=seqnames, selected=x$seqnames)
#     } else if(nrow(x)>1){
#       z <- table(x$seqnames)
#       updateNumericInput(session, "chrStart2", value=min(x$start))
#       updateRadioButtons(session, "locationChrs2", choices=seqnames, selected=names(z)[which.max(z)])
#     }
#   })
#   
  glymacols <- c("seqnames", "start", "end", "ID")
  
  # Output a data table of glyma IDs matching the text box
  output$glymaTable <- renderDataTable({
    if(input$glymaID!=""){
      x <- id()
      
      if(nrow(x)>0) {
        y <- x[,glymacols]
        names(y) <- c("Chromosome", "Start", "End", "GlymaID")
        y
      } else {
        data.frame(Problem = "Query not found")
      }

    } else {
      data.frame(Hint = "Enter a partial Glyma ID to search")
    }
  }#, options=list(pageLength=25, 
#                   lengthMenu= "[ [10, 25, 50, -1], [10, 25, 50, 'All'] ]", 
#                  searchDelay=100,
#                   bsort=FALSE,
#                  sDom='<"top"i>rt<"bottom"p><"clear">')
)
  
  # Output a data table of glyma IDs matching the text box
  output$glymaTable2 <- renderDataTable({
    if(nchar(input$glymaID2)>0){
      x <- id2()
      if(nrow(x)>0){
        y <- x[,glymacols[-c(2:3)]]
        names(y) <- c("Chr", "GlymaID")
        y
      } else {
        data.frame(Problem = "Query not found")
      }
    } else {
      data.frame(Hint = "Enter a glymaID at the top")
    }
  }, searchDelay=250, 
  options=list(pageLength=25, 
               lengthMenu= "[ [10, 25, 50, -1], [10, 25, 50, 'All'] ]", 
               bsort=FALSE,
               sDom='<"top"i>rt<"bottom"p><"clear">')
)
  
  # Output a data table of displayed SNPs + Varieties
  output$snpTable <- renderDataTable({
    if(length(input$glymaID2)>0 & nrow(id2())>=1){
      x <- varsnps()

      if(nrow(x)>0) {
        x
      } else {
        data.frame(Result="No SNPs found for query", query=input$glymaID2)
      }
    } else {
      data.frame(Hint="Enter a glymaID in the left panel", Problem="No glymaIDs match your query")
    }
  }, searchDelay=250, 
  options=list(pageLength=25, 
               lengthMenu= "[ [10, 25, 50, -1], [10, 25, 50, 'All'] ]", 
               bsort=FALSE,
               sDom='<"top"i>rt<"bottom"p><"clear">')
  )
  
  # get all snps before the current start value, 
  # arranged so that the top row is the closest to the current start value
  rev.snps <- reactive(filter(snp.counts, Chromosome==input$locationChrs) %>% 
                       filter(Position < as.numeric(input$chrStart)) %>%
                       arrange(desc(Position)))
  
  # Scan up genome
  observe({
    # depend on input$up
    if(input$up==0){
      return()
    }
    
    isolate({
      tmp <- rev.snps()
      
      # If there are more than 4*bases rows, then get the position of the 4*base-th row
      if(nrow(tmp)>4*input$bases){
        updateTextInput(session, inputId = "chrStart",
                        value=tmp$Position[4*input$bases]-1)
      } else { 
        # otherwise, find the minimum position
        updateTextInput(session, inputId = "chrStart",
                        value=min(tmp$Position))
      }
    })
  })
  # Scan down genome
  observe({
    # depend on input$down
    if(input$down==0){
      return()
    }
    
    isolate({
      tmp <- snps()
      if(nrow(tmp)<4*input$bases){
        # if there aren't 4*bases rows, then keep current TextInput value
        return()
      } else {
        # otherwise, move downstream by input$bases snps.
        updateTextInput(session, inputId = "chrStart", value=tmp$Position[4*input$bases]-1)
      }
    })
  })
  
  
  # Expression that generates a plot of snps along the specified range of 
  # the chosen chromosome. 
  output$AggregatePlot <- renderPlot({
    if(length(input$chrStart)>0 & length(input$locationChrs)>0){
      if(input$glymaID!=""){
        matchingGlymas <- filter(id(), shown)$ID
        plot.title <- switch(min(length(matchingGlymas), 2)+1, input$glymaID, matchingGlymas, paste(matchingGlymas, collapse=", "))
      } else {
        plot.title=paste0(input$chrStart, " on ", input$locationChrs)
      }
      
      
      tmp <- snps()
      loc.chrs <- input$locationChrs
      loc.start <- input$chrStart
      
      if(nrow(tmp)>0){
        snps.sub <- tmp[1:min(nrow(tmp), 4*input$bases),]
        plot <- ggplot(data=snps.sub) + 
          geom_histogram(aes(x=factor(Position), 
                             fill=factor(Nucleotide, levels=c("A", "G", "T", "C")), 
                             weight=Count), 
                         position="stack") + 
          scale_fill_manual("Nucleotide", values=c("A"=pal[1], "G"=pal[2], 
                                                   "T"=pal[3], "C"=pal[4]),
                            drop=FALSE) + 
          xlab(paste("Position on ", gsub("Chr", "Chromosome ", loc.chrs))) + 
          ylab("Total number of SNPs in 10 lines") + 
          theme_bw() + 
          theme(axis.text.x=element_text(angle=90), legend.position="bottom") + 
          ggtitle(paste0("SNPs downstream of ", plot.title))
        
      } else{
        plot <- ggplot() + 
          geom_text(aes(x=0, y=0, label=paste("No SNPs found downstream of ", 
                                              as.numeric(loc.start), " on ", 
                                              gsub("Chr", "Chromosome ", loc.chrs)))) +         
          theme_bw() + 
          theme(axis.text=element_blank(), 
                axis.ticks=element_blank(), 
                axis.title=element_blank())
      }
    } else {
      plot <- ggplot() + 
        geom_text(aes(x=0, y=0, label="Please select chromosome and an index on that chromosome\n\n Note: It may take a minute to process the input")) +         
        theme_bw() + 
        theme(axis.text=element_blank(), 
              axis.ticks=element_blank(), 
              axis.title=element_blank())
    }
    
    print(plot)
  })
  
  # Expression that generates a plot of snps for a given gene/glymaID, 
  # facetted by variety
  output$VarietySnpPlot <- renderPlot({
    
    if(length(input$glymaID2)>0){
      id2.cur <- id2()

      if(nrow(id2.cur)>0){
        matchingGlymas <- filter(id2.cur, shown)$ID
                
        tmp <- varsnps()
        if(nrow(tmp)>0){
          if(length(input$varieties)>0){
            tmp$Variety <- factor(tmp$Variety, levels=input$varieties)
          } else {
            tmp$Variety <- factor(tmp$Variety, levels=varieties)
            tmp$Variety <- droplevels.factor(tmp$Variety)
            if(length(levels(tmp$Variety))>10){
              tmp <- filter(tmp, Variety %in% levels(tmp$Variety)[1:10])
              tmp$Variety <- droplevels.factor(tmp$Variety)
            }
          }
          
          plot.title=switch(min(length(matchingGlymas), 2)+1, input$glymaID2, matchingGlymas, paste(matchingGlymas, collapse=", "))
          
          plot <- ggplot(data=tmp) + 
            geom_histogram(aes(x=factor(Position), 
                               fill=factor(Alternate, levels=c("A", "G", "T", "C")), 
                               weight=Alt_Allele_Freq)) + 
            scale_fill_manual("Nucleotide", values=c("A"=pal[1], "G"=pal[2], 
                                                     "T"=pal[3], "C"=pal[4]),
                              drop=FALSE) + 
            facet_grid(.~Variety, drop=FALSE) +
            coord_flip() + 
            xlab("Position on Chromosome") + 
            ggtitle(paste0("SNPs for gene ", plot.title)) +
            theme_bw() + 
            theme(axis.text.x=element_text(angle=90)) + 
            scale_y_continuous("Alternate Allele Frequency", breaks=c(0, 1, 2), limits=c(0,2))
        } else{
          labeldf <- data.frame(x=0, y=0, label=paste0("No SNPs found for glymaID(s) containing\n", input$glymaID2))
          plot <- ggplot() + 
            geom_text(data=labeldf, aes(x=x, y=y, label=label)) +         
            theme_bw() + 
            theme(axis.text=element_blank(), 
                  axis.ticks=element_blank(), 
                  axis.title=element_blank())
        }
      } else {
        plot <- ggplot() + 
          geom_text(aes(x=0, y=0, label="Please enter a valid glymaID\n\n 
                        Note: It may take a minute to process the input")) +         
          theme_bw() + 
          theme(axis.text=element_blank(), 
                axis.ticks=element_blank(), 
                axis.title=element_blank())
      }
    } else {
      plot <- ggplot() + 
        geom_text(aes(x=0, y=0, label="Please enter a glymaID\n\n Note: It may take a minute to process the input")) +         
        theme_bw() + 
        theme(axis.text=element_blank(), 
              axis.ticks=element_blank(), 
              axis.title=element_blank())
    }
    
    print(plot)
  })
  
  # Expression that generates a plot of snp density along the chromosome
  output$DensityPlot <- renderPlot({
    if(length(input$densityChrs)>0 & length(input$densityVars)>0){
      tmp <- snpDensity()
      if(nrow(tmp)>0){
        plot <- ggplot(data=tmp) + 
          geom_line(aes(x=x, y=y)) + 
          xlab("Position on Chromosome") + 
          ylab("Density") + 
          theme_bw() + 
          ggtitle("SNP Density") + 
          theme(axis.text.x=element_text(angle=45, hjust=1)) + 
          facet_grid(Variety~Chromosome, scales="free_x")
      } else{
        labeldf <- data.frame(x=0, y=0, label=paste("No SNPs found for \n(", 
                                                    paste(input$densityVars, collapse=", "),
                                                    ")\n on \n Chromosomes (", 
                                                    paste(gsub("Chr", "", input$densityChrs), 
                                                          collapse=", "), ").", sep=""))
        plot <- ggplot() + 
          geom_text(data=labeldf, aes(x=x, y=y, label=label)) +         
          theme_bw() + 
          theme(axis.text=element_blank(), 
                axis.ticks=element_blank(), 
                axis.title=element_blank())
      }
    } else {
      plot <- ggplot() + 
        geom_text(aes(x=0, y=0, label="Please select chromosome(s) and varieties of interest\n\n Note: It may take a minute to process the input")) +         
        theme_bw() + 
        theme(axis.text=element_blank(), 
              axis.ticks=element_blank(), 
              axis.title=element_blank())
    }
    
    print(plot)
  })
  
  output$KinshipSNP <- renderUI({
    tagList(
      singleton(tags$head(tags$script(src = "animint/animint.js", type='text/javascript'))),
      singleton(tags$head(tags$script(src = "animint/vendor/d3.v3.js", type='text/javascript'))),
      singleton(tags$head(tags$link(rel = "stylesheet", type = "text/css", href="styles.css"))),
      tags$div(id="plot"),
      tags$script("var plot = new animint('#plot','animint/plot.json');")
      )
  })
})
