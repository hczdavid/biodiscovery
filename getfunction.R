# ------------------------------------------------------------------#
#
#      Bio-Discovery R Shiny Tool Support code: functional connection
#
#      Author:  Caizhi David Huang
#
# -------------------------------------------------------------------#



getfunctiongene <- function(all_test, cut_FC){

  selectprobe_FC <-all_test$probe[abs(all_test[,2]) > cut_FC[1] & all_test[,3] < cut_FC[2]]

  selectgene_FC  <- AnnotationDbi::select(hgu219.db, keys=selectprobe_FC,  columns=c("SYMBOL"), keytype="PROBEID")[,2] %>% unique()

  selectgene_FC

}


# getfunctiongene <- function(all_test, cut_FC, cut_TRT){
#   
#   selectprobe_FC <-all_test$probe[abs(all_test$effect_FC) > cut_FC[1] & all_test$p_FC < cut_FC[2]]
#   selectprobe_TRT <-all_test$probe[abs(all_test$effect_TRT) > cut_TRT[1] & all_test$p_TRT < cut_TRT[2]]
#   
#   selectgene_FC  <- AnnotationDbi::select(hgu219.db, keys=selectprobe_FC,  columns=c("SYMBOL"), keytype="PROBEID")[,2] %>% unique()
#   selectgene_TRT <- AnnotationDbi::select(hgu219.db, keys=selectprobe_TRT, columns=c("SYMBOL"), keytype="PROBEID")[,2] %>% unique()
#   
#   list(selectgene_FC, selectgene_TRT)
#   
# }


getfunction <- function(selectgene_FC, title){

  
  ego_FC  <- enrichGO(gene= selectgene_FC,  OrgDb =  get("org.Hs.eg.db"), ont= "BP",pAdjustMethod = "BH", keyType = "SYMBOL") 

  ego_res_FC   <- ego_FC@result

  ego_FC@result$Description <- ifelse(nchar(ego_FC@result$Description) < 40, ego_FC@result$Description, substr(ego_FC@result$Description, 1,40))
  
  pFC  <- clusterProfiler::dotplot(ego_FC, showCategory = 10, font.size=14, title = title)

  list(plot = pFC, table = ego_res_FC)
  
}


# 
# getfunction <- function(selectgene_FC, selectgene_TRT){
#   
#   
#   ego_FC  <- enrichGO(gene= selectgene_FC,  OrgDb =  get("org.Hs.eg.db"), ont= "BP",pAdjustMethod = "BH", keyType = "SYMBOL") 
#   ego_TRT <- enrichGO(gene= selectgene_TRT, OrgDb =  get("org.Hs.eg.db"), ont= "BP",pAdjustMethod = "BH", keyType = "SYMBOL") 
# 
#   ego_res_FC   <- ego_FC@result
#   ego_res_TRT  <- ego_TRT@result
# 
#   ego_FC@result$Description <- ifelse(nchar(ego_FC@result$Description) < 40, ego_FC@result$Description, substr(ego_FC@result$Description, 1,40))
#   ego_TRT@result$Description <- ifelse(nchar(ego_TRT@result$Description) < 40, ego_TRT@result$Description, substr(ego_TRT@result$Description, 1,40))
# 
#   
#   pFC  <- clusterProfiler::dotplot(ego_FC, showCategory = 10, font.size=14, title = "Based on fold change")
#   pTRT <- clusterProfiler::dotplot(ego_TRT, showCategory = 10, font.size=14, title = "Based on treatment effect")
# 
#   list(plot = list(FC = pFC, TRT = pTRT), table = list(FC = ego_res_FC, TRT = ego_res_TRT))
#   
# }

