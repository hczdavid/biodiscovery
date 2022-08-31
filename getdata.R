# -------------------------------------------------------------------#
#
#      Bio-Discovery R Shiny Tool support code: load data to app
#
#      Author:  Caizhi David Huang
#      Project: 2022 Summer Intern Project
# -------------------------------------------------------------------#



library(data.table)
  
getdata <- function(techtype, tissuetype, inputprobe, trt, ifmerge){
  
  
  path <- "Rapp Data"
  
  
  allfile    <-  list.files(path, full.names = T)
  selectfile <- allfile[grepl(c(techtype), allfile) &  grepl(paste0("_",c(tissuetype)), allfile, ignore.case = T)]
  
  
  # load the log FC data
  logFC           <- fread(file = selectfile[1])
  logFC           <- as.data.frame(logFC)
  rownames(logFC) <- logFC[,1]
  logFC           <- logFC[,-1]
  
  # load the log gene expression data
  logGE           <- fread(file = selectfile[2])
  logGE           <- as.data.frame(logGE)
  
  rownames(logGE) <- logGE[,1]
  logGE           <- logGE[,-1]
  
  # load the thickness difference and thickness data
  thickdiff <- read.csv(file = selectfile[3])[,-1]
  thickness <- read.csv(file = selectfile[4])[,-1]
  
  
  thickness$Study <- sapply(strsplit(thickness$Study, split = "_"), `[`, 1) # remove p1 p2
  thickness$name  <- paste(thickness$Study, thickness$SubID, thickness$Trt, sep = "_")
  thickdiff       <- thickdiff %>% left_join(thickness[,c(3,4,6)], by = c("name" = "name")) # match the sample ID and colname
  
  
  # get all or subset of probes
  if(inputprobe != "all"){
    logFC <- logFC[inputprobe, ]
    logGE <- logGE[inputprobe,]
  }
  
  # get the data by trt 
  if(!ifmerge & length(trt) > 1){
    
    mydata <- list()
    for(i in 1:length(trt)){
      
      # get the retinol in the first two study
      thickdiff_trt <- thickdiff[thickdiff$Trt %in% trt[i], ]
      
      # get the sample ID for the baseline
      # thickdiff_trt$namev <- paste(thickdiff_trt$Study.x, thickdiff_trt$SubID, "Vehicle", sep = "_")
      
      thickdiff_trt$namev <- ifelse(thickdiff_trt$Study %in% c("GSS2842", "GSS2843"), paste(thickdiff_trt$Study, thickdiff_trt$SubID, "SC800", sep = "_"),
                                    paste(thickdiff_trt$Study, thickdiff_trt$SubID, "Vehicle", sep = "_"))
      
      
      thickdiff_trt       <- thickdiff_trt %>% left_join(thickness[,c(3,6)], by = c("namev" = "name"))
      thickdiff_trt       <- thickdiff_trt[thickdiff_trt$group != "Gray",]
      thickdiff_trt$name <- gsub("[/+]", ".", thickdiff_trt$name)
      
      # get the gene expression at baseline and FC
      logGE_baseline_trt <- logGE[, thickdiff_trt$Sample.y]
      logFC_trt          <- logFC[, thickdiff_trt$name]
      
      
      colnames(thickdiff_trt)[7] <- "SampleID_trt"
      colnames(thickdiff_trt)[10] <- "SampleID_vehicle"
      
      mydata[[i]] <- list(foldchange = logFC_trt, baseline = logGE_baseline_trt, metadata = thickdiff_trt)
      
    }
    names(mydata) <- trt
    return(mydata)
    
  }else{
    # get the retinol in the first two study
    thickdiff_trt <- thickdiff[thickdiff$Trt %in% trt, ]
    
    # get the sample ID for the baseline
    # thickdiff_trt$namev <- paste(thickdiff_trt$Study.x, thickdiff_trt$SubID, "Vehicle", sep = "_")
    
    thickdiff_trt$namev <- ifelse(thickdiff_trt$Study %in% c("GSS2842", "GSS2843"), paste(thickdiff_trt$Study, thickdiff_trt$SubID, "SC800", sep = "_"),
                                  paste(thickdiff_trt$Study, thickdiff_trt$SubID, "Vehicle", sep = "_"))
    
    
    thickdiff_trt       <- thickdiff_trt %>% left_join(thickness[,c(3,6)], by = c("namev" = "name"))
    thickdiff_trt       <- thickdiff_trt[thickdiff_trt$group != "Gray",]
    #thickdiff_trt$name <- gsub("[/+]", ".", thickdiff_trt$name)
    
    # get the gene expression at baseline and FC
    logGE_baseline_trt <- logGE[, thickdiff_trt$Sample.y]
    logFC_trt          <- logFC[, thickdiff_trt$name]
    
    
    colnames(thickdiff_trt)[7] <- "SampleID_trt"
    colnames(thickdiff_trt)[10] <- "SampleID_vehicle"
    
    mydata <- list(list(foldchange = logFC_trt, baseline = logGE_baseline_trt, metadata = thickdiff_trt))
    names(mydata) <- paste0(trt, collapse = ".")
    return(mydata)
    
  }
  
}
