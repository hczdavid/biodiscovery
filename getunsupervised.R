# -----------------------------------------------------------------#
#
#      Bio-Discovery R Shiny Tool support code : Supervised learning
#
#      Author:  Caizhi David Huang
#
# -----------------------------------------------------------------#




library(tsne)
library(snifter)
library(umap)
library(scales)
library(plotly)



getunsupervised <- function(data, meta, title){
  

  pco_fc  <- prcomp(t(data), center = TRUE, scale. = TRUE )
  pco_fc  <- summary(pco_fc)
  
  perplexity <- ncol(data)/3 - 1
  tsne_fc <- fitsne(t(data), perplexity = perplexity) %>% round(2)
  colnames(tsne_fc) <- c("t-SNE 1", "t-SNE 2")
  
  umap_fc    <- umap(t(data))
  umap_fc    <- round(umap_fc$layout,2)
  colnames(umap_fc) <- c("Umap 1", "Umap 2")
  
  
  plotdata1 <- pco_fc$x[,1:2] %>% round(2) %>%  as.data.frame() 
  plotdata1 <- cbind(plotdata1, meta, tsne_fc, umap_fc)
  
  plotdata1$thickdiff <- round(plotdata1$thickdiff, 2)

  ngroup <- unique(meta$group) %>% length()
  
  fig_pca <- plot_ly(
    plotdata1, x = ~PC1, y = ~PC2, color = ~group, colors = c(hue_pal()(ngroup)),
    text = ~paste("ThickDiff:", thickdiff, '<br>Name:', name, "<br>Group:", group),
    marker = list(size = 10)
  ) %>% layout(title = title)

  
  fig_tsne <- plot_ly(
    plotdata1, x = ~`t-SNE 1`, y = ~`t-SNE 2`, color = ~group, colors = c(hue_pal()(ngroup)),
    text = ~paste("ThickDiff:", thickdiff, '<br>Name:', name, "<br>Group:", group),
    marker = list(size = 10)
  ) %>% layout(title = title)
  
  fig_umap <- plot_ly(
    plotdata1, x = ~`Umap 1`, y = ~`Umap 2`, color = ~group, colors = c(hue_pal()(ngroup)),
    text = ~paste("ThickDiff:", thickdiff, '<br>Name:', name, "<br>Group:", group),
    marker = list(size = 10)
  ) %>% layout(title = title)
  
  
  
  
  list(pca = fig_pca, tsne = fig_tsne, umap = fig_umap)
}





gethclust <- function(data, meta, title){
  
  dist1 <- dist(t(data), method = "euclidean", diag = T, upper = T, p = 2)
  
  
  metadata <- cbind(pheno = meta$group, sample = colnames(data))
  
  colLab<<-function(n){
    if(is.leaf(n)){
      
      #I take the current attributes
      a=attributes(n)
      
      #I deduce the line in the original data, and so the treatment and the specie.
      ligne=base::match(attributes(n)$label,metadata[,2])
      treatment=metadata[ligne,1];
      if(treatment=="Reponse"){col_treatment=scales::hue_pal()(2)[2]};if(treatment=="No reponse"){col_treatment=scales::hue_pal()(2)[1]}
      
      #Modification of leaf attribute
      attr(n,"nodePar")<-c(a$nodePar,list(cex=1.5,lab.cex=1,pch=20,col=col_treatment,lab.col=col_treatment,lab.font=1,lab.cex=1))
    }
    return(n)
  }
  
  
  dend  <- hclust(dist1, method = "ward.D") %>% as.dendrogram(hang = 0.1)
  dL    <- dendrapply(dend, colLab)
  plot(dL , main = paste("Hierarchical Clustering", title))
  legend("topright", 
         legend = c("No response" , "Response"), 
         col = scales::hue_pal()(2), 
         pch = c(20,20), bty = "n",  pt.cex = 1.5, cex = 1.2 , 
         text.col = "black", horiz = FALSE, inset = c(0.05, 0.01))
  #p <- recordPlot()
}
