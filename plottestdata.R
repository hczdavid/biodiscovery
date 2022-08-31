# ------------------------------------------------------------------------#
#
#      Bio-Discovery R Shiny Tool support code: plot the testing box plot
#
#      Author:  Caizhi David Huang
#
# ------------------------------------------------------------------------#





plottestdata <- function(FC, metadata, all_test, s){
  
  i <- s
  rownames(all_test) <- all_test$probe
  pn <- all_test$probe
  
  
  pheno <- metadata$group
  
  genename <- sapply(strsplit(all_test$`Gene Symbol`,split = "///"), `[`,1)
  
  plotdata <- data.frame(FC = t(FC[pn[i],])[,1], pheno = pheno)
  plotdata$pheno <- ifelse(plotdata$pheno == "Reponse", "With Response", "No Response")
  plotdata$FCname <- genename[s]


  p_FC_lable <- ifelse(all_test$p_FC[i] < 1e-4, format(all_test$p_FC[i], digits = 2,nsmall = 2, scientific = TRUE), round(all_test$p_FC[i],3))
  p_TRT_lable <- ifelse(all_test$p_TRT[i] < 1e-4, format(all_test$p_TRT[i], digits = 2,nsmall = 2, scientific = TRUE), round(all_test$p_TRT[i],3))
  
  pFC <- ggplot(plotdata, aes(x = pheno, y = FC)) + 
    geom_boxplot(width = 0.4, outlier.shape = NA) +
    geom_jitter(width = 0.2) + 
    geom_hline(yintercept = 0, linetype = "dashed", col = "red") + 
    geom_hline(yintercept = mean(plotdata$FC), linetype = "dashed", col = "blue") + 
    annotate(geom = 'text', label =  paste0(" p_FC (boxplot) = ", p_FC_lable), x = -Inf, y = Inf, hjust = 0, vjust = 1.5, size = 5) + 
    annotate(geom = 'text', label =  paste0(" p_TRT (dashedline) = ", p_TRT_lable), x = -Inf, y = Inf, hjust = 0, vjust = 3.2, size = 5) + 
    xlab("Group") + ylab("Fold-Change") + 
    facet_wrap(~FCname) + 
    theme_bw()+ 
    theme(axis.title = element_text(size = 15), axis.text = element_text(size = 12), strip.text = element_text(size = 15))
  
  
  pFC
}



# VOLCANO plot

plotvolcano <- function(data, fc_cut, p_cut){
  
  cols   <- c("up" = "#ffad73", "down" = "#26b3ff", "ns" = "grey") 
  sizes  <- c("up" = 2, "down" = 2, "ns" = 1.5) 
  alphas <- c("up" = 1, "down" = 1, "ns" = 0.5)
  
  colnames(data)[2:3] <- c("effect_FC", "p_FC")
  
  myda <- data %>%  mutate(gene_type = ifelse(effect_FC >= fc_cut & p_FC < p_cut, "up",
                              ifelse(effect_FC <= -fc_cut & p_FC < p_cut, "down", "ns")))
  
  myp <-  ggplot(myda, aes(x = effect_FC,
                           y = -log10(p_FC),
                           fill = gene_type,    
                           size = gene_type,
                           alpha = gene_type,
                           col = gene_type)) + 
    geom_point(shape = 21) + 
    geom_hline(yintercept = -log10(0.05),
               linetype = "dashed") + 
    geom_vline(xintercept = c(-fc_cut, fc_cut),
               linetype = "dashed") +
    scale_color_manual(values = cols) + 
    scale_fill_manual(values = cols) + # Modify point colour
    scale_size_manual(values = sizes) + # Modify point size
    scale_alpha_manual(values = alphas) +# Modify point transparency
    theme_bw() + xlab("log2 (effect size)") + ylab("-log10(p-value)")  +
    theme(axis.title = element_text(size = 20), title = element_text(size = 15), axis.text = element_text(size = 15), legend.position = "none")
  
  myp
  
  
}
