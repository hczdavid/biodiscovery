

# -----------------------------------------------------------------#
#
#      Bio-Discovery R Shiny Tool support code : Supervised learning
#
#      Author:  Caizhi David Huang
#
# -----------------------------------------------------------------#



library(randomForest)
library(pROC)
library(foreach)
library(doParallel)
library(caret)



getprediction <- function(cdata_FC, pheno, testtable, topnumber){
  
  
  # control      <- trainControl(method="repeatedcv", number=5, repeats=5,search = "random")
  # hyper.model1 <- train(cdata_FC,as.factor(pheno),method="rf",ntree = 1000,tuneLength=15,trControl=control)
  # bestmtry <- as.matrix(hyper.model$bestTune)[1]
  
  if(is.null(testtable)){
    probe_FC <- colnames(cdata_FC)
    #probe_BL <- colnames(cdata_BL)
  }else{
    
    probe_FC  <- testtable[, c("probe", "p_FC")] %>% arrange(p_FC) %>% pull(probe) %>% `[`(1:topnumber)
    #probe_BL  <- testtable[, c("probe", "p_BL")] %>% arrange(p_BL) %>% pull(probe) %>% `[`(1:topnumber)

  }
  
  
  registerDoParallel(detectCores()-1);detectCores()
  
  result <- foreach(i = 1:length(pheno), .inorder = FALSE)%dopar%{
    
    library(randomForest)
    set.seed(i)
    
    
    
    rf.train.FC   <- randomForest(cdata_FC[-i, probe_FC], y = as.factor(pheno[-i]))
    rf.pred.FC <- predict(rf.train.FC, cdata_FC[i, probe_FC] , type="prob")[,2]
    
    # rf.train.BL   <- randomForest(cdata_BL[-i, probe_BL], y = as.factor(pheno[-i]))
    # rf.pred.BL <- predict(rf.train.BL, cdata_BL[i, probe_BL] , type="prob")[,2]
    # 
    return(c(rf.pred.FC))
    
  }
  stopImplicitCluster()
  
  
  rf.pred.FC  <- sapply(result, `[`, 1)
  #rf.pred.BL  <- sapply(result, `[`, 2)

  list(rf.pred.FC)
}


getroc <- function(rf.pred.FC, pheno){
  
  roc_FC  <- roc(pheno, rf.pred.FC)
  #roc_BL  <- roc(pheno, rf.pred.BL)
  
  auc_FC    <- round(roc_FC$auc, 3)
  #auc_BL    <- round(roc_BL$auc, 3)
  
  ggroc(list(FC = roc_FC),aes=c("color"), legacy.axes = TRUE, size = 2)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          panel.border = element_rect(fill = NA,size = 1), axis.title = element_text(size=18),
          legend.position = c(0.9,0.15),legend.title = element_blank(),axis.text = element_text(size=15),
          legend.text = element_text(size=15), legend.key.size = unit(1, "cm"))+
    scale_color_manual(values = scales::hue_pal()(1), label = c(paste0("BL (AUC: ", auc_FC, ")")))

}


getpredtable <- function(rf.pred.FC, meta){
  
  
  meta$pred_FC <- ifelse(rf.pred.FC > 0.5, "Responder", "Non-responder")
  #meta$pred_BL <- ifelse(rf.pred.BL > 0.5, "Response", "No response")
  meta$group <- ifelse(meta$group == "Reponse", "Responder", "Non-responder")
  meta
}