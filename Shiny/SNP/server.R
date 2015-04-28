# Libraries required for the app
#-------------------------------------------------------------------------------
library(shiny)
library(ggplot2) # required for plots throughout the app
library(ggenealogy) # required for family tree display
library(dplyr) # required for table operations throughout the app
library(reshape2) 
library(stringr) # required for string operations throughout the app

# End libraries
#-------------------------------------------------------------------------------

# Supplemental libraries required to modify the app
#-------------------------------------------------------------------------------
# library(RColorBrewer) # required for color scheme generation 

# End supplemental libraries
#-------------------------------------------------------------------------------


# Initial datasets
#-------------------------------------------------------------------------------
# Data required for this app:
# - initial lists of varieties and chromosomes (varieties, seqnames)
# - GlymaIDList:
#     GlymaIDList is a data table of all GlymaIDs and positions. 
#     Columns: seqnames, source, feature, gene, start, end, score, strand, frame, 
#              group (list of attributes), ID (glymaID, with .Wm82.a2.v1 appended),
#              name (glymaID without thet appended text), Parent (for nested IDs), 
#              chrnum (raw number), searchstr (lowercase string of the form 00g000000),
#              and link (HTML link to soybase)
# - snpList:
#     snpList is a data table of all SNPs for each variety
#     Columns: Chromosome, Position, Reference, Alternate, Allele.Freq, Variety, Alt.Allele.Freq, 
#              Genotype_Probability, Gene_State (0|1), Alt_Allele_Count
# - Several other summary variables are also required but can be pre-computed from previously mentioned variables
#   - snpList.GlymaSummary:
#       snpList.GlymaSummary is a data table which summarizes SNPs by GlymaID
#       Columns: Chromosome, ID, TotalVarietiesWithSNPs, NumberOfSnps, chrnum, searchstr
#   - snpList.PositionSummary:
#       snpList.PositionSummary is a data table summarizing SNPs by position and GlymaID
#       Columns: Chromosome, Position, ID, Number.of.Varieties, chrnum, searchstr
#   - snp.counts:
#       snp.counts is a data frame with the aggregate (population total) 
#       number of snps at each position with each allele
#       Columns: Chromosome, Position, Reference, Nucleotide, Count
#   - snp.density:
#       snp.density is a data table containing the density of SNPs for each chromosome and variety.
#       Columns: Chromosome, Variety, x, y

load("ShinyStart.rda")

objlist <- ls()
if(!"snp.counts"%in%objlist){
  load("SNPCounts.RData")
  snp.counts <- snp.counts %>% arrange(Chromosome, Position) %>% group_by(Chromosome) 
}

if(!"GlymaIDList"%in%objlist){
  load("GlymaID.rda")
  idx <- names(GlymaIDList) %in% c("seqnames", "feature", "start", "end", "ID", "numid", "chrnum", "searchstr", "link")
  GlymaIDList <- GlymaIDList[,idx] %>% group_by(chrnum, numid) 
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

# End of dataset initialization
#-------------------------------------------------------------------------------

# Initial variables and functions
#-------------------------------------------------------------------------------
# Genealogy data and plotting function
source("plotFamilyTree.R")

# Options for DataTables
tableoptions <- list(searchDelay=250, pageLength=10, sDom='lti<br>rp')

# Palette for ATGC data
# pal <- brewer.pal(8, "Paired")[c(1, 2, 7, 8)]
pal <- c("#A6CEE3", "#1F78B4", "#FDBF6F", "#FF7F00") # output from brewer.pal line

# Function to fix variety names for links
fixVarieties <- function(x){
  y <- gsub("PI ?", "PI ", x)
  y <- gsub("LG ?", "LG ", y)
  gsub(".", "", gsub(" ", "%20", y, fixed=T), fixed=T)
}

# End initial variable/function defs
#-------------------------------------------------------------------------------

# Empty plots with relevant messages
#-------------------------------------------------------------------------------
# Plot to show with no GlymaID
emptyGlymaPlot <- ggplot() + 
  geom_text(aes(x=0, y=0, label="Please enter a glymaID\n\n Note: It may take a minute to process the input")) +         
  theme_bw() + 
  theme(axis.text=element_blank(), 
        axis.ticks=element_blank(), 
        axis.title=element_blank())

# Plot to show with invalid GlymaID
invalidGlymaPlot <- ggplot() + 
  geom_text(aes(x=0, y=0, label="Please enter a valid glymaID\n\n Note: It may take a minute to process the input")) +         
  theme_bw() + 
  theme(axis.text=element_blank(), 
        axis.ticks=element_blank(), 
        axis.title=element_blank())

# End empty plots
#-------------------------------------------------------------------------------

# Define server logic required to generate and plot a subset of varieties with a subset of chromosomes
shinyServer(function(input, output, session) {

  # Code to deal with query string
  observe({
    withProgress(
      message = "Parsing Query", 
      min=0, max=1, value=.05,
      {
        str <- parseQueryString(session$clientData$url_search)
        fix.vecs <- names(str)[grepl("c\\(.*\\)", str)]
        # ensure variable names and values are preserved while removing programmatic stuff (injection prevention)
        for(i in fix.vecs){
          
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
                               label = "Cultivars of Interest (up to 10)", 
                               choices = sort(unique(varieties)), 
                               selected = str$varieties)
        }
        if("chromosomes"%in%names(str)){
          updateSelectInput(session = session, 
                            inputId = "locationChrs",
                            label = "Choose Chromosome of Interest", 
                            selected = str$chromosomes[1])
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
      }
    )# End progress message
  })
  
  # Return SNPs for selected varieties
  varsnps <- reactive({
    withProgress({
      id.cur <- id()
      tmp <- filter(id.cur, shown)
      
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
    }, message="Retrieving SNPs", value=0)
    })
  
  getRelatives <- reactive({
    validate(
      need(!is.null(input$variety0) & input$variety0 != "", 
           "Select at least one variety!"),
      need(input$gens, 
           "Select the number of generations to search!")
    )
    withProgress({
      relatives <- data.frame(label=input$variety0, gen=0, stringsAsFactors=F)
      if(input$ancestors){
        anc <- try(getAncestors(input$variety0, tree, input$gens), silent=T)
        if(!is.character(anc)){
          relatives <- rbind(relatives, anc)
          relatives$gen <- -1*relatives$gen
        }
      }
      if(input$descendants){
        desc <- try(getDescendants(input$variety0, tree, input$gens), silent=T)
        if(!is.character(desc)){
          relatives <- rbind(relatives, desc)
        }
      }
      
      names(relatives) <- c("Variety", "generation")
      return(relatives %>% arrange(generation) %>% as.data.frame())
    }, message = "Retrieving Relatives", value=0)
  })
  
  # Return SNPs for varieties related to selected variety
  varsnps.genealogy <- reactive({
    withProgress(
      {
        id.cur <- id()
        tmp <- filter(id.cur, shown)
        
        chr <- paste0("Chr", substr(tmp$ID, 7, 8))
        position.min=min(tmp$start)
        position.max=max(tmp$end)
        tmplist <- as.data.frame(
          filter(snpList, Chromosome == chr) %>%
            filter(Position >= position.min & Position <= position.max))
        
        if(nrow(tmplist)>0){
          relatives <- getRelatives() %>% filter(Variety%in%varieties)
          tmplist <- filter(tmplist, Variety%in%relatives$Variety) %>% left_join(relatives)
        }
        tmplist
      }, message="Retrieving SNPs", value=0)
  })
  
  # filter snps by chromosome and input location
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
                                Chromosome%in%c(input$glymaChrs, "")) %>%
                         filter(Variety%in%c(input$varieties, "")))
  
  # Match glyma text to ID
  match <- function(glymaIDstring){
    # strip off extra glyma stuff to get 01g000000
    tmp <- gsub("wm82a2v1", "", gsub("glyma", "", gsub(".", "", tolower(glymaIDstring), fixed=TRUE)))
    # attempt to match strings nicely
    chr.str <- word(tmp, sep="[Gg]")
    pos.str <- word(tmp, sep="[Gg]", start=-1)
    if(nchar(chr.str)==2 & nchar(pos.str)<=6){
      possible.match <- with(GlymaIDList, str_detect(chrnum, chr.str) & str_detect(numid, pos.str))
    } else{
      possible.match <- with(GlymaIDList, str_detect(searchstr, tmp))
    }
    which(possible.match)
  }
  
  # search for glyma IDs
  id <- reactive({
    withProgress({
      validate(
        need(!is.null(input$glymaID), "Input part of a glymaID to begin")
      )
      if(nchar(input$glymaID)>0){
        res <- data.frame(GlymaIDList[match(input$glymaID),])
        if(nrow(res)>0){
          res$shown <- FALSE
          res$shown[res$numid[1]==res$numid & res$chrnum[1]==res$chrnum] <- TRUE
        }      
      } else {
        res <- data.frame(seqnames=NULL, start=NULL, end=NULL, ID=NULL, numid=NULL, chrnum=NULL, searchstr=NULL, link=NULL, shown=NULL)
      }
      res
    }, message="Finding Glyma IDs", value=0)
  })
  
  # search for glyma IDs to filter snpList.GlymaSummary
  id3 <- reactive({
    withProgress({
      validate(
        need(!is.null(input$glymaID3), "Input part of a glymaID to begin")
      )
      
      if(nchar(input$glymaID3)>0){
        res <- data.frame(GlymaIDList[match(input$glymaID3),])
      } else {
        res <- data.frame(seqnames=NULL, start=NULL, end=NULL, ID=NULL, numid=NULL, chrnum=NULL, searchstr=NULL, link=NULL, shown=NULL)
      }
      
      if(length(input$glymaChrs)>0){
        if(nrow(res)>0){
          res <- res[res$seqnames%in%input$glymaChrs,]
        } else {
          res <- GlymaIDList[GlymaIDList$seqnames%in%input$glymaChrs,]
        }
      } 
      
      if(nrow(res)>0){
        res$shown <- FALSE
        res$shown[res$numid[1]==res$numid & res$chrnum[1]==res$chrnum] <- TRUE
      } else {
        res <- data.frame(seqnames=NULL, start=NULL, end=NULL, ID=NULL, numid=NULL, chrnum=NULL, searchstr=NULL, link=NULL, shown=NULL)
      }
      res
    }, message="Finding GlymaIDs", value=0)
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

  glymacols <- c("seqnames", "link", "start", "end")
  
  glymaTableOutput <- reactive({
    withProgress(
      {
        validate(
          need(nchar(input$glymaID)>0, "Input a glymaID or partial glymaID to start.")
        )
        
        x <- id()
        
        validate(
          need(nrow(x)>0, "Query not found")
        )
        
        # idx <- if(input$tabname!="Aggregated SNPs") which(names(x)%in%glymacols[-c(3:4)]) else which(names(x)%in%glymacols)
        idx <- which(names(x)%in%glymacols)
        
        y <- x[,idx]
        names(y) <- c("Chr", names(x)[2], "Start", "End", names(x)[5:8], "GlymaID")[idx]
        y <- y  %>% transform(Chr = str_replace(Chr, "Chr0?", ""), 
                              GlymaID=str_replace(GlymaID, ">\\s?Glyma\\.", ">"))
        showCols <- if(input$tabname!="Aggregated SNPs") c("GlymaID", "Start", "End") else c("GlymaID", "Start", "End")
        y[,showCols]
      }, message="Making Glyma Table", value=.5)
  })
  
  # Output a data table of glyma IDs matching the text box
  output$glymaTable <- renderDataTable({
    glymaTableOutput()
  }, escape=FALSE, options=tableoptions)
  
  # Output a data table of glyma IDs matching the text box
  output$glymaTable2 <- renderDataTable({
    glymaTableOutput()
  }, escape=FALSE, options=tableoptions)
  
  # Output a data table of displayed SNPs + Varieties
  output$snpTable <- renderDataTable({
    if(length(input$glymaID)>0 & nrow(id())>=1){
      x <- varsnps()

      if(nrow(x)>0) {
        x$Variety <- sprintf("<a href='http://www.ars-grin.gov/cgi-bin/npgs/html/acc_list_post.pl?lopi=&hipi=&lono=&hino=&plantid=%s&pedigree=&taxon=%s&family=&cname=&country=&state=&site=%s&acimpt=%s&uniform=%s&recent=anytime&pyears=1&received=&records=100' target='_blank'>%s</a>", fixVarieties(x$Variety), "Glycine%20max", "ALL%20-%20All%20Repositories", "Any%20Status", "Any%20Status", x$Variety)
        x <- x[,c("Chromosome", "Position", "Reference", "Alternate", "Variety", "Alt_Allele_Count", "Alt_Allele_Freq")]
        names(x)[3:4] <- c("Ref", "Alt")
        names(x)[6:7] <- c("Alt Allele Count", "Alt Allele Freq")
        x
      } else {
        data.frame(Result="No SNPs found for query", query=input$glymaID)
      }
    } else {
      data.frame(Hint="Enter a glymaID in the left panel", Problem="No glymaIDs match your query")
    }
  }, escape=FALSE, options=tableoptions)

  # Output a data table of SNPs tabulated by Chromosome and GlymaID
  output$glymaSummary <- renderDataTable({
    if(length(input$glymaChrs)>0){
      x <- input$glymaChrs
    } else {
      x <- seqnames
    }
    
    if(nchar(input$glymaID3)>0){
      gid <- id3()$ID
      matches <- rowSums(sapply(gid, str_detect, string=snpList.GlymaSummary$ID))>0
      res <- unique(bind_rows(snpList.GlymaSummary[matches,],snpList.GlymaSummary[snpList.GlymaSummary$Chromosome%in%x,] ))
    } else {
      res <- snpList.GlymaSummary[with(snpList.GlymaSummary, Chromosome%in%x),]
    }
    
    if(nrow(res)>0){
      res2 <- as.data.frame(res[,c("Chromosome", "link", "TotalVarietiesWithSNPs", "NumberOfSNPs")])
      names(res2)[2:4] <- c("ID", "Total Varieties with SNPs", "Number of SNP sites")
      res2
    } else {
      data.frame(ID="No Matching GlymaIDs found", Chromosome="?")
    }
  }, escape=FALSE, options=tableoptions)
  
  # Output a data table of tabulated SNPs and # Varieties for each glymaID
  output$positionSummary <- renderDataTable({
    if(length(input$glymaChrs)>0){
      x <- input$glymaChrs
    } else {
      x <- seqnames
    }
    
    if(nchar(input$glymaID3)>0){
      gid <- id3()$ID
      res <- snpList.PositionSummary[with(snpList.PositionSummary, Chromosome%in%x & str_detect(ID, gid[1])),]
    } else {
      res <- snpList.PositionSummary[with(snpList.PositionSummary, Chromosome%in%x),]
      res <- res[1:min(5000, nrow(res)),]
    }
    
    if(nrow(res)>0){
      res2 <- as.data.frame(res[,c("Chromosome", "Position", "link", "Number.of.Varieties")])
      names(res2)[3:4] <- c("ID", "Number of Varieties")
      res2 %>% arrange(Chromosome, Position)
    } else {
      data.frame(ID="No Matching GlymaIDs found", Chromosome="?")
    }
  }, escape=FALSE, options=tableoptions)

  # Scan up genome
  observe({
    if(input$tabname=='Aggregate SNP Browser'){
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
    }
  })
  # Scan down genome
  observe({
    if(input$tabname=='Aggregate SNP Browser'){
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
    }
  })
  
  
  # Expression that generates a plot of snps along the specified range of 
  # the chosen chromosome. 
  output$AggregatePlot <- renderPlot({
    if(length(input$chrStart)>0 & length(input$locationChrs)>0){
      withProgress(
        message="Working", value=0, 
        {
          if(input$glymaID!=""){
            matchingGlymas <- filter(id(), shown)$ID
            plot.title <- switch(min(length(matchingGlymas), 2)+1, input$glymaID, matchingGlymas, paste(matchingGlymas, collapse=", "))
          } else {
            plot.title=paste0(input$chrStart, " on ", input$locationChrs)
          }
          
          
          tmp <- snps()
          loc.chrs <- input$locationChrs
          loc.start <- input$chrStart
          
          if(nrow(tmp)>0){ # Plot displayed SNPs
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
            
          } else { # End branch for displaying SNPs, start "no SNPs"
            plot <- ggplot() + 
              geom_text(aes(x=0, y=0, label=paste("No SNPs found downstream of ", 
                                                  as.numeric(loc.start), " on ", 
                                                  gsub("Chr", "Chromosome ", loc.chrs)))) +         
              theme_bw() + 
              theme(axis.text=element_blank(), 
                    axis.ticks=element_blank(), 
                    axis.title=element_blank())
          } # End "No SNPs" plot option
        }) # End progress indicator
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
    
    if(length(input$glymaID)>0){
      id.cur <- id()

      if(nrow(id.cur)>0){
        matchingGlymas <- filter(id.cur, shown)$ID
                
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
          
          plot.title=switch(min(length(matchingGlymas), 2)+1, input$glymaID, matchingGlymas, paste(matchingGlymas, collapse=", "))
          
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
          labeldf <- data.frame(x=0, y=0, label=paste0("No SNPs found for glymaID(s) containing\n", input$glymaID))
          plot <- ggplot() + 
            geom_text(data=labeldf, aes(x=x, y=y, label=label)) +         
            theme_bw() + 
            theme(axis.text=element_blank(), 
                  axis.ticks=element_blank(), 
                  axis.title=element_blank())
        }
      } else {
        plot <- invalidGlymaPlot
      }
    } else {
      plot <- emptyGlymaPlot
    }
    
    print(plot)
  })
  
  # Expression that generates a plot of snps for a given gene/glymaID and a single
  # input variety (and its parents or children)
  output$GenealogySnpPlot <- renderPlot({
    
    if(length(input$glymaID)>0){
      id.cur <- id()
      
      if(nrow(id.cur)>0){
        matchingGlymas <- filter(id.cur, shown)$ID
        
        rels <- getRelatives()
        tmp <- varsnps.genealogy()
        tmp <- filter(tmp, Variety%in%varieties)
        if(nrow(tmp)>0){
          tmp$Variety <- factor(tmp$Variety, levels=unique(tmp$Variety[order(tmp$generation, tmp$Variety, decreasing = F)]))
          
          title.phrase <- ifelse(
            input$ancestors, 
            ifelse(input$descendants, 
                   'relatives of ', # True, True
                   'ancestors of '),  # True, False
            ifelse(input$descendants, 
                   'descendants of ', # False, True
                   '') # False, False
            )
          
          plot.title=switch( 
            min(length(matchingGlymas), 2)+1, 
            # default value
            paste0(title.phrase, "cultivar ", input$variety0, " on gene ", input$glymaID), 
            # if only one glymaID, use that as the title
            paste0(title.phrase, "cultivar ", input$variety0, " on gene ", matchingGlymas),
            # Otherwise, paste the ids together
            paste0(title.phrase, "cultivar ", input$variety0, " on genes ", 
                  paste(matchingGlymas, collapse=", "))
          )
          
          tmp$Label <- paste("Gen ", tmp$generation, "\n", tmp$Variety, sep="")
          # Varieties + generations for which there is data (even if no SNPs identified)
          data.vars <- paste("Gen ", rels$generation, "\n", rels$Variety, sep="")[rels$Variety%in%varieties]
          tmp$Label <- factor(tmp$Label, levels=data.vars)
          
          plot <- ggplot(data=tmp) + 
            geom_histogram(aes(x=factor(Position), 
                               fill=factor(Alternate, levels=c("A", "G", "T", "C")), 
                               weight=Alt_Allele_Freq)) + 
            scale_fill_manual("Nucleotide", values=c("A"=pal[1], "G"=pal[2], 
                                                     "T"=pal[3], "C"=pal[4]),
                              drop=FALSE) + 
            facet_grid(.~ Label, drop=FALSE) +
            coord_flip() + 
            xlab("Position on Chromosome") + 
            ggtitle(paste0("SNPs for ", plot.title)) +
            theme_bw() + 
            theme(axis.text.x=element_text(angle=90)) + 
            scale_y_continuous("Alternate Allele Frequency", breaks=c(0, 1, 2), limits=c(0,2))
        } else{
          labeldf <- data.frame(x=0, y=0, label=paste0("No SNPs found for glymaID(s) containing\n", input$glymaID, "\nwithin ", input$gens, " generations of ", input$variety0))
          plot <- ggplot() + 
            geom_text(data=labeldf, aes(x=x, y=y, label=label)) +         
            theme_bw() + 
            theme(axis.text=element_blank(), 
                  axis.ticks=element_blank(), 
                  axis.title=element_blank())
        }
      } else {
        plot <- invalidGlymaPlot
      }
    } else {
      plot <- emptyGlymaPlot
    }
    
    print(plot)
  })
  
  # Expression that generates a family tree of selected varieties
  output$GenealogyTree <- renderPlot({
    plotFamilyTree(input$variety0, input$gens, ncol=1) 
   })
  
  # Expression that generates a plot of snp density along the chromosome
  output$DensityPlot <- renderPlot({
    if(length(input$glymaChrs)>0 & length(input$varieties)>0){
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
                                                    paste(input$varieties, collapse=", "),
                                                    ")\n on \n Chromosomes (", 
                                                    paste(gsub("Chr", "", input$glymaChrs), 
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
