# -------------------------------------------------------------------#
#
#      Bio-Discovery R Shiny Tool support code: differential abundance
#
#      Author:  Caizhi David Huang
#
# -------------------------------------------------------------------#




performTest <- function(FC, metadata, method){
  
  
  pheno <- metadata$group
  
  
  # Fold change test
  FC_test <- apply(FC, 1, function(x){
    
    withres <- x[pheno == "Reponse"]
    nores   <- x[pheno == "No reponse"]
    
    logfc <- mean(withres) - mean(nores)
    
    if(method == "t test"){
      wp    <- t.test(withres, nores)$p.value
    }else{
      wp    <- wilcox.test(withres, nores)$p.value
    }
    
    c(effect_FC = logfc, p_FC = wp)
    
  })
  FC_test <- t(FC_test) %>% as.data.frame() 
  
  
  # treatment test
  TRT_test <- apply(FC, 1, function(x){
    
    logfc <- mean(x)
    
    
    if(method == "t test"){
      wp    <- t.test(x)$p.value
    }else{
      wp    <- wilcox.test(x)$p.value
    }
    
    
    c(effect_TRT = logfc, p_TRT = wp)
    
  })
  
  TRT_test <- t(TRT_test) %>% as.data.frame() 
  
  
  all_test <- cbind(probe = rownames(FC_test), FC_test, TRT_test) %>% arrange(p_FC, p_TRT)
  
  
  # add the gene symbol
  annofile  <- data.table::fread("Rapp Data/probe2gene.csv")
  all_test  <- all_test %>% left_join(annofile, by = c("probe" = "Probe Set ID"))
  all_test
}






























