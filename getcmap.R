# -------------------------------------------------------------------#
#
#      Bio-Discovery R Shiny Tool support code: material connection
#
#      Author:  Caizhi David Huang
#
# -------------------------------------------------------------------#

get_GSEA <- function(up, down, ref){
  
  GSEA.score <- function (x, n){
    x <- sort(x)
    l <- length(x)
    x1 <- (1:l)/l-x/n
    a <- max(x1)
    b <- 1/l - min(x1)
    score <- ifelse(b<a, a, -b)
    score
  }
  
  nprobe <- nrow(ref)
  
  upref <- ref[up,,drop = F]
  doref <- ref[down,, drop = F]
  
  upscore <- apply(upref, 2, GSEA.score, nprobe)
  doscore <- apply(doref, 2, GSEA.score, nprobe)
  
  toscore <- upscore-doscore
  toscore[(upscore*doscore) > 0 ] <- 0
  
  mys <- as.data.frame(cbind(up_score = upscore, down_score = doscore, total_score = toscore))
  mys[order(mys[,3], decreasing = TRUE) ,]
}
get_meta <- function(myss, chem, datalib){
  
  myss$chem <- chem[rownames(myss),]
  myss$ChipID <- rownames(myss)
  myss <- myss %>% left_join(datalib)
  naindex <- which(is.na(myss$CellLine))
  
  if(length(naindex) > 0){
    
    
    for(i in 1:length(naindex)){
      
      batchid <- stringr::str_sub(myss$ChipID[naindex[i]], start = 6,end = 8) %>% as.numeric()
      myss$CellLine[naindex[i]] <- unique(datalib[datalib$Batch == batchid,]$CellLine)
      myss$Che.Name[naindex[i]] <- "DMSO"
      myss$Batch[naindex[i]] <- batchid
      
    }
  }
  myss
}

# SIMsip
SimSIP     <- function(r1,r2){
  hF<-1/((r1)^0.5)
  hG<-1/((r2)^0.5)
  Texs<-hF%*%hG
  return(Texs)
}
get_SimSip <- function(query, ref){
  
  query <- query[rownames(ref)]
  
  query <- rank(abs(query)*-1)
  
  myscore <- apply(ref, 2, SimSIP, query)
  myscore <- sort(myscore, decreasing = T) %>% as.data.frame()
  colnames(myscore) <- "score"
  myscore
}


getcmap <- function(meanFC, cutoff = 500){
  
  
  load("Rapp Data/cmap_database.Rdata")
  
  allrank_ref1_nodmso <- allrank_ref1_nodmso[names(meanFC),]
  
  probename <- names(meanFC)
  
  
  corder    <- rank(meanFC*-1)
  upprobe   <- probename[which(corder <= cutoff)]
  downprobe <- probename[which(corder > nrow(allrank_ref1_nodmso) - cutoff)]
  
  mys    <- get_GSEA(upprobe, downprobe, allrank_ref1_nodmso)
  mys2   <- get_meta(mys, chem_name_nodmso, data_lib)

  
  SIM_score  <-  get_SimSip(meanFC, allrank_ref1_nodmso)
  SIM_score2 <-  get_meta(SIM_score, chem_name_nodmso, data_lib)

  list(gsa = mys2[, c(5:7, 1:3)], sim = SIM_score2[,c(3,4,5, 1)])
}