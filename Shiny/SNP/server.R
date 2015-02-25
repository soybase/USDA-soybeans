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

if(!"snpList.GlymaSummary"%in%objlist | !"snpList.PositionSummary"%in%objlist){
  load("GlymaSNPsummary.rda")
  names(snpList.GlymaSummary) <- gsub("Number.of.Varieties", "TotalVarietiesWithSNPs", names(snpList.GlymaSummary))
  names(snpList.GlymaSummary) <- gsub("Number.of.SNP.Sites", "NumberOfSNPs", names(snpList.GlymaSummary))
  snpList.GlymaSummary <- snpList.GlymaSummary %>% group_by(Chromosome, ID)
  snpList.PositionSummary <- snpList.VarietySummary %>% group_by(Chromosome, Position, ID)
  rm(snpList.VarietySummary)
}

if(!"snp.density"%in%objlist){
  load("SNPDensity.RData")
}

# Data required for this app:
## snp.df - data frame containing information from the info() and geno() commands applied
##          to the VCF file. 

# Define server logic required to generate and plot a subset of varieties with a subset of chromosomes
shinyServer(function(input, output, session) {

  observe({
    str <- parseQueryString(session$clientData$url_search)
    fix.vecs <- names(str)[grepl("c\\(.*\\)", str)]
    for(i in fix.vecs){
      # ensure variable names and values are preserved while removing programmatic stuff (injection prevention)
      values <- gsub(
        # Remove any character used for executing a function or storing values into variables
        "[\"\'()<->=\\*\\+/]?", "", 
        unlist( # Extract only valid R names/strings used in this app
          str_extract_all(str[[i]], "[\"\']{1}[[:alnum:]\\._ ]{1,}[\"\']{1}")
        ))
      if(length(values)>0){
        str[[i]] <- eval(parse(text=paste0("c(", paste(sprintf("'%s'", values), collapse=",", sep=""), ")")))
      } else {
        warning(sprintf("Possible injection attempt detected: Query String = %s", str[[i]]))
        str[[i]] <- NULL
        warning("Injection attempt removed from input successfully")
      }
    }
    
    # Reset each input using the appropriate update* function
    if("tabname"%in%names(str)){
      session$sendCustomMessage(type='setTab', str$tabname)
    }
    if("varieties"%in%names(str)){
      updateSelectizeInput(session = session, 
                           inputId = "varieties", 
                           label = "Choose up to 10 Cultivars of Interest", 
                           choices = sort(unique(varieties)), 
                           selected = str$varieties)
      updateSelectizeInput(session = session, 
                           "densityVars", "Choose Varieties of Interest", 
                           choices=sort(unique(varieties)), 
                           selected = str$varieties)
    }
    if("chromosomes"%in%names(str)){
      updateSelectInput(session = session, 
                        inputId = "locationChrs",
                        label = "Choose Chromosome of Interest", 
                        selected = str$chromosomes[1])
      updateSelectizeInput(session = session, 
                           "densityChrs", "Choose Chromosome(s) of Interest", 
                           choices=unique(seqnames), 
                           selected = str$chromosomes)
      updateSelectizeInput(session = session, 
                           "glymaChrs", "Filter GlymaIDs by Chromosome(s)", 
                           choices=unique(seqnames), 
                           selected = str$chromosomes)
    }
    if("glymaID"%in%names(str)){
      updateTextInput(session = session, 
                      "glymaID", 
                      "Locate Position by Glyma ID",
                      value=str$glymaID)
      updateTextInput(session = session, 
                      "glymaID2", 
                      "Locate Position by Glyma ID",
                      value=str$glymaID)
      updateTextInput(session = session, 
                      "glymaID3", 
                      "View SNP sites within a GlymaID",
                      value=str$glymaID)
    }
    if("chrStart"%in%names(str)){ 
      updateTextInput(session = session, 
                      "chrStart", 
                      "Start point", value=str$chrStart)
    }
    if("bases"%in%names(str)){
      if(as.numeric(str$bases)>=5 & as.numeric(str$bases)<=50){
        updateNumericInput(session = session, 
                           "bases",
                           "# downstream SNPs (up to 50)", 
                           value=as.numeric(str$bases), min=5, max=50, step=5)
      }
    }
  })
  
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
                   filter(Position>=min(as.numeric(input$chrStart), Inf)) %>% 
                   arrange(Chromosome, Position))
  
  # get all snps before the current start value, 
  # arranged so that the top row is the closest to the current start value
  rev.snps <- reactive(filter(snp.counts, Chromosome==input$locationChrs) %>% 
                         filter(Position < as.numeric(input$chrStart)) %>%
                         arrange(desc(Position)))
  
  #  filter density data reactively
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
  
  # search for glyma IDs to filter snpList.GlymaSummary
  id3 <- reactive({
    if(nchar(input$glymaID3)>0){
      # strip off extra glyma stuff to get 01g000000
      tmp <- gsub("wm82a2v1", "", gsub("glyma", "", gsub(".", "", tolower(input$glymaID3), fixed=TRUE)))
      # attempt to match strings nicely
      chr.str <- word(tmp, sep="[Gg]")
      pos.str <- word(tmp, sep="[Gg]", start=-1)
      if(nchar(chr.str)==2 & nchar(pos.str)<=6){
        row.idx <- chr.str==tolower(as.data.frame(select(GlymaIDList, chrnum))$chrnum) & grepl(pos.str, as.data.frame(select(GlymaIDList, numid))$numid)
        res <- as.data.frame(GlymaIDList)[row.idx,]
      } else if(length(input$glymaChrs)>0){
        row.idx <- grepl(tmp, as.data.frame(select(GlymaIDList, searchstr))$searchstr)
        res <- as.data.frame(GlymaIDList)[row.idx,] %>% filter(seqnames%in%input$glymaChrs)        
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
      updateSelectInput(session, "locationChrs", choices=seqnames, selected=x$seqnames)
    } else if(nrow(x)>1){
      z <- table(x$seqnames)
      updateNumericInput(session, "chrStart", value=min(x$start))
      updateSelectInput(session, "locationChrs", choices=seqnames, selected=names(z)[which.max(z)])
    }
  })

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
  })
  
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

  # Output a data table of tabulated SNPs + Varieties for each glymaID
  output$glymaSummary <- renderDataTable({
    if(length(input$glymaChrs)>0){
      x <- input$glymaChrs
    } else {
      x <- seqnames
    }
    
    if(nchar(input$glymaID3)>0){
      gid <- id3()$ID
      res <- snpList.GlymaSummary %>% filter(Chromosome%in%x) %>% 
        filter(str_detect(ID, gid))
      if(nrow(res)>0){
        res[,1:4]
      } else {
        data.frame(ID="No Matching GlymaIDs found", Chromosome="?")
      }
    } else {
      filter(snpList.GlymaSummary, Chromosome%in%x)[,1:4]
    }
    
  }, searchDelay=250)
  
  # Output a data table of tabulated SNPs + Varieties for each glymaID
  output$positionSummary <- renderDataTable({
    if(length(input$glymaChrs)>0){
      x <- input$glymaChrs
    } else {
      x <- seqnames
    }
    
    if(nchar(input$glymaID3)>0){
      gid <- id3()$ID
      res <- snpList.PositionSummary %>% filter(Chromosome%in%x) %>% 
        filter(str_detect(ID, gid)) %>% arrange(ID, Chromosome, Position)
      if(nrow(res)>0){
        res[,1:4]
      } else {
        data.frame()
      }
    } else {
      res <- filter(snpList.PositionSummary, Chromosome%in%x)%>% select(1:4) %>% arrange(ID, Chromosome, Position)
      res[1:min(50000, nrow(res)),]
    }
  }, searchDelay=250)
  
#   output$varietySummary <- renderDataTable({
#     if(length(input$glymaChrs)>0){
#       x <- input$glymaChrs
#     } else {
#       x <- seqnames
#     }
#     
#     if(length(input$glymaPosition)>0){
#       if(nchar(input$glymaID3)>0){
#         gid <- id3()$ID
#         res <- GlymaIDSNPs %>% filter(Chromosome%in%x) %>% 
#           filter(str_detect(ID, gid)) %>% 
#           filter(Position==as.numeric(input$glymaPosition)) %>%
#           arrange(ID, Chromosome, Position, Variety)
#       } else {
#         res <- GlymaIDSNPs %>% filter(Chromosome%in%x) %>% 
#           filter(Position==as.numeric(input$glymaPosition)) %>%
#           arrange(ID, Chromosome, Position, Variety)
#       }
#     } else {
#       res <- data.frame(problem = "Please Enter a GlymaID or Position", Variety = "No Varieties Found")
#     }
#     
#     res   
#     
#   }, searchDelay=250)
#   
  # Scan up genome
  observe({
    # depend on input$up
    if(input$up==0){
      return()
    }
    
    isolate({
      tmp <- rev.snps()
      
      # If there are more than 2*bases rows, then get the position of the 2*base-th row
      if(nrow(tmp)>2*input$bases){
        updateTextInput(session, inputId = "chrStart",
                        value=tmp$Position[2*input$bases]-1)
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
      if(nrow(tmp)<2*input$bases){
        # if there aren't 2*bases rows, then keep current TextInput value
        return()
      } else {
        # otherwise, move downstream by input$bases snps.
        updateTextInput(session, inputId = "chrStart", value=tmp$Position[2*input$bases]-1)
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
        pos.idx <- unique(tmp$Position)
        pos.idx <- sort(as.numeric(pos.idx))[1:min(length(pos.idx), input$bases)]
        snps.sub <- filter(tmp, Position%in%pos.idx)
        plot <- ggplot(data=snps.sub) + 
          geom_histogram(aes(x=factor(Position), 
                             fill=factor(Nucleotide, levels=c("A", "G", "T", "C")), 
                             weight=Count), 
                         position="stack") + 
          scale_fill_manual("Nucleotide", values=c("A"=pal[1], "G"=pal[2], 
                                                   "T"=pal[3], "C"=pal[4]),
                            drop=FALSE) + 
          xlab(paste("Position on ", gsub("Chr", "Chromosome ", loc.chrs))) + 
          ylab(paste0("Total number of SNPs in ", length(varieties)," lines")) + 
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
