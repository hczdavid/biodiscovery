# -----------------------------------------------------------------#
#
#      Bio-Discovery R Shiny Tool 
#
#      Author:  Caizhi David Huang
#
# -----------------------------------------------------------------#


library(BiocManager)
options(repos = BiocManager::repositories())


library(shiny)
library(DT)
library(tidyverse)
library(data.table)
require(openxlsx)
library(ggpubr)
library(clusterProfiler)
library(hgu219.db)
library(shinyjs)
library(randomForest)
library(pROC)
library(foreach)
library(doParallel)
library(tsne)
library(snifter)
library(umap)
library(scales)
library(plotly)



source("getdata.R")
source("performTest.R")
source("plottestdata.R")
source("getfunction.R")
source("getcmap.R")
source("getunsupervised.R")
source("getprediction.R")


ui <- fluidPage(

  # Title
  titlePanel(title=div(strong("Bio-Discovery Tool"), img(src="bios.jpg", height = 50), img(src="dms.png", height = 45))),
  
    # Use the main and side bar layout
    sidebarLayout(
      
        # --------- Start of side bar --------- #
      
        sidebarPanel(
          
            fileInput("infile", label = "Please upload new data:"),
            a(href="exampledata.csv", "(Example data)", download=NA, target="_blank"),
            
            tags$hr(),
            
            h3(strong("ACS Datasets")),
            
            tags$h3(),
            
            selectInput("trt", label = "Select one or multiple treatment(s):", choices = c(
              "Retinol", "Promatrixyl" ,"EGF+RF","Peptide","No Trt" ,"RP/SC800","EGF+RP+RF"), multiple = TRUE),
            checkboxInput("ifmerge", label = "merge multiple treatments", value = TRUE),
            
            tags$hr(),
            
           
            selectInput("tissuetype", "Choose tissue type:", choices = c("Epidermis", "Dermis"),selected = "Epidermis", selectize = T),
            selectInput("techtype", "Choose technical type:",choices = c("LCM", "ATC"), selected  =  c("LCM")),
           
            
            tags$h1(),
            
            textAreaInput("probe", label = "Input the gene or probe list (optional):"),
            checkboxInput("ifallprobes", label = "use all genes/probes", value = TRUE),
            
            tags$hr(),
            
            
            div("  click", downloadLink('downloadrawdata', 'here'), " to download raw data"),
            
            tags$h3(),
            span(strong("Perform analysis at different tabs!"), style = "color:blue"),
            
            width = 2 # 2 out of 12 
       
        ),
        
        # ************************************************************************ #
        
        # --------- Start of main panel --------- #
        
        #  1. Instruction
        #  2. Differential Abundance
        #  3. Functional Connection
        #  4. Material Connection
        #  5. Underlying Pattern (unsupervised learning)
        #  6. Response Prediction (supervised learning)
        
        mainPanel(
          
            # --------- Start of panel sets (6 panels) --------- #
    
            tabsetPanel(
              
                #  1. Instruction panel
                tabPanel( h4(strong(span("Instruction", style = "color:red"))),
                          
                          tags$h3(),

                          h2(strong("Welcome to use this Bio-Discovery Tool (BDT)!")),
                          tags$h3(),
                          h3(strong("BDT provides the analysis strategy (workflow) to get insights of multiple anti-aging clinical screening (ACS) studies.")),
                          tags$h3(),
                          img(src="scope.PNG", height = 550) #
                ),
                
                #  1.5. Data panel
                tabPanel( h4(strong(span("Data", style = "color:red"))),
                          
                          tags$h3(),
                          
                          h2(strong("Example Data: Anti-aging Clinical Screening (ACS) Datasets")),
                          tags$h3(),
                          span(h3(strong("Paired Gene Expression Data")), style = "color:green"),
                          
                          #h3(strong("Paired Gene Expression Data")),
                          img(src="genedata.PNG", height = 550), #genedata.png
                          tags$h3(),
                          span(h3(strong("Clinical Phenotype Data")), style = "color:green"),
                          
                          #h3(strong("Clinical Phenotype Data")),
                          img(src="phenotype.PNG", height = 420) #
                           
                ),
                
                #  2. Differential Abundance panel
                tabPanel(h4(strong("Differential Abundance")),
                         tags$h3(),
                         h2(strong("Explore the Gene Behavior between Treatment and Skin Epidermis Thickness Change")),
                         
                         tags$h3(),
                         img(src="DA1.PNG", height = 120),
                         br(),
                         h4("Note:"),
                         h4("- p-values is un-adjusted."),
                         h4("- select one row to display the boxplot of the fold change."),
                         h4("- click gene name to corresponding gene card webpage."),
                         
                         selectInput("da_method", "Choose method:", choices = c("t test", "Wilcoxon Rank Sum"),selected = "t test"),
                         
                         span(h4(strong("Current methods: paired/unpaired t-test, Wilcoxon rank-sum test")), style = "color:red"),
                         
                        actionButton( "run","Run Analysis",  icon = icon("bar-chart-o"), class = "btn-success"),
                        tags$h6(),
                        div("click", downloadLink('downloadDatatest', 'here'), " to download result table and ", downloadLink('downloadboxplot', 'here'), " for boxplot."),
                        tags$hr(),
                        splitLayout(  DTOutput("alldata",),
                                      plotOutput("alldataplot", height = "520px", width = "90%"), cellWidths = c("60%", "40%"))
                        ),

                #  3. Functional Connection panel
                tabPanel(h4(strong("Functional Connection")),

                         tags$h3(),
                         h2(strong("Understand the Mechanism of How the Gene Expression Affect Skin Behavior")),
                         
                         span(h4("Optional: adjust the thresholds to select the significant genes."), style = "color:green"),
                        
                         span(h4(strong("Current methods: over-representation analysis")), style = "color:red"),
                         
                         tags$h6(),
                         div(actionButton( "runfun","Run Analysis", icon = icon("bar-chart-o"), class = "btn-info"),
                             span("(please run differential abundance first!)"), style = "color:red"),
                         
                         tags$h6(),
                         div("click", downloadLink('downloadDataFC', 'here'), " to download results"),

                         tags$h2(),

                         tabsetPanel(
                           
                           tabPanel("Treatment Effect",
                                    h3("Functions related to treatment effect"),
                                    tags$h2(),
                                    splitLayout(cellWidths = c("12%", "12%"),
                                                sliderInput("effetcut_trt","Set the effect size cutoff:",0, 1, step = 0.05, value = 0.3,width = '100%'),
                                                sliderInput("pcut_trt","Set the p value cutoff:",0, 0.1, step = 0.001, value = 0.05,width = '100%')),
                                    
                                    
                                    tabsetPanel(
                                      
                                      tabPanel("Selected genes", 
                                               
                                               splitLayout(  
                                                 DTOutput("trt_table_gene"), 
                                                 plotOutput("trt_volcano", height = "500px", width = "90%"), cellWidths = c("30%", "40%"))
                                               ),
                                      
                                      
                                      tabPanel("function result",
                                               splitLayout(  
                                                 DTOutput("trt_table"), 
                                                 plotOutput("trt_plot", height = "500px", width = "90%"), cellWidths = c("60%", "40%")))
                                      
                                    )),
                           
                           tabPanel("Fold Change",
                                    h3("Functions related to fold change"),

                                    tags$h2(),

                                    splitLayout(cellWidths = c("12%", "12%"),
                                                sliderInput("effetcut_fc","Set the effect size cutoff:",0, 1, step = 0.05, value = 0.3,width = '100%'),
                                                sliderInput("pcut_fc","Set the p value cutoff:",0, 0.1, step = 0.001, value = 0.05,width = '100%')),

                                    
                                    tabsetPanel(
                                      
                                      tabPanel("Selected genes", 
                                               splitLayout(  
                                                 DTOutput("fc_table_gene"), 
                                                 plotOutput("fc_volcano", height = "500px", width = "90%"), cellWidths = c("30%", "40%"))
                                               
                                      ),
                                      
                                      
                                      tabPanel("function result",
                                               splitLayout(  
                                                 DTOutput("fc_table"),
                                                 plotOutput("fc_plot", height = "500px", width = "90%"), cellWidths = c("60%", "40%")))
                                      
                                    ))
                  
                                    
                                   
                           
                           
                         ),

                ),
                
                #  4. Material Connection panel
                tabPanel(h4(strong("Material Connection")),
                         tags$h3(),
                         h2(strong("Search Similar Chemicals in Cell Line Database")),
                         
                         tags$h3(),
                         img(src="cmap.PNG", height = 300),

                         span(h4("Cell line chemical effect database"), style = "color:green"),
                         span(h4(" - 1506 chemicals"), style = "color:green"),
                         span(h4(" - 2 cell lines (tert keratinocytes and BJ fibroblasts)"), style = "color:green"),
                         tags$h3(),
                         
                         sliderInput("cmapgene","Set number of probes in Cmap:",300, 2000, step = 100, value = 500,width = '15%'),
                         
                         span(h4(strong("Current methods: CMap, SimSip")), style = "color:red"),
                         
                         actionButton( "runcmap","Run Analysis", icon = icon("bar-chart-o"), class = "btn-primary"),
                         tags$h6(),
                         div("click", downloadLink('downloadDatamc', 'here'), " to download results"),

                         DTOutput("camptable", width = "60%")
                ),
                
                #  5. Underlying Pattern (unsupervised learning) panel
                tabPanel(h4(strong("Underlying Pattern")),
                         tags$h3(),
                         h2(strong("Discover Unknow Data Pattern")),
                         tags$h3(),
                  
                         selectInput("datatype", "Choose data type:", choices = c("Original data", "Rank"),selected = "Original data"),
                         
                         span(h4(strong("Current methods: PCA, t-SNE, Umap, hierarchical clustering")), style = "color:red"),
                         
                         actionButton( "rununsuper","Run Analysis", icon = icon("bar-chart-o"), class = "btn-warning"),
                         tags$h6(),
                         div("click", downloadLink('downloadDataup', 'here'), " to download results"),
                         
                         
                         tabsetPanel(
                           tabPanel("PCA",
                                    tags$h2(),
                                    
                                    splitLayout(  plotlyOutput("pcaplot_fc", height = "500px", width = "95%"),
                                                  plotOutput("none1"),
                                                  plotlyOutput("pcaplot_bl", height = "500px", width = "95%"), cellWidths = c("40%", "10%", "40%"))),
                           
                           tabPanel("t-SNE",
                                    tags$h2(),
                                    
                                    splitLayout(  plotlyOutput("tsneplot_fc", height = "500px", width = "95%"),
                                                  plotOutput("none2"),
                                                  plotlyOutput("tsneplot_bl", height = "500px", width = "95%"), cellWidths = c("40%", "10%", "40%"))),
                           tabPanel("Umap",
                                    tags$h2(),
                                    
                                    splitLayout(  plotlyOutput("umapplot_fc", height = "500px", width = "95%"),
                                                  plotOutput("none3"),
                                                  plotlyOutput("umapplot_bl", height = "500px", width = "95%"), cellWidths = c("40%", "10%", "40%"))),
                           tabPanel("hclust",
                                    tags$h2(),
                                    
                                    splitLayout(  plotOutput("hcplot_fc", height = "500px", width = "95%"),
                                                  plotOutput("none4"),
                                                  plotOutput("hcplot_bl", height = "500px", width = "95%"), cellWidths = c("40%", "10%", "40%")))
                           
                           
                         )
                         
                         
                ),
                
                #  6. Response Prediction (supervised learning) panel
                tabPanel(h4(strong("Prediction Capability")),
                         
                         tags$h3(),
                         h2(strong("Predict the Phenotypical Result using Genetic Information")),
                         tags$h3(),
                         
                         
                         selectInput("datatype_rp", "Choose data type:", choices = c("Original data", "Rank"),selected = "Original data"),
                         sliderInput("numbergene","Use top genes from DA analysis (optional):",1000, 49293, step = 1000, value = 3000, width = '15%'),
                         
                         tags$h3(),
                         span(h4(strong("Current methods: random forest")), style = "color:red"),
                         
                         actionButton( "runpred","Run Analysis", icon = icon("bar-chart-o"), class = "btn-warning"),
                         tags$h6(),
                         div("click", downloadLink('downloadDatapred', 'here'), " to download results"),
                         
                         splitLayout(cellWidths = c("50%", "50%"),
                           DTOutput("predtable", width = "80%"),
                           plotOutput("roc", height = "600px", width = "100%"))
                         
                         
                )
                
            ), 
            
            width = 10 # 10 out of 12
    
        )
    )
)



# ************************************************************************ #

server <- function(input, output, session) {

  
  # ---------------------------------------- get the input new data if uploaded ---------------------------- #
  
  inputdata <- reactive({
    inFile <- input$infile
    if (is.null(inFile)) return(NULL)
    indata <- data.table::fread(inFile$datapath)
    colnames(indata) <- c("probe", "fold-change", "p-value")
    annofile  <- data.table::fread("Rapp Data/probe2gene.csv")
    indata  <- indata %>% inner_join(annofile, by = c("probe" = "Probe Set ID"))
  })
  
  # output$alldata <- renderDataTable({
  #   
  #   req(input$infile)
  #   
  #   inputdata <- inputdata()
  # 
  #   alllink <- sapply(strsplit(inputdata$`Gene Symbol`,split = "///"), `[`,1)
  #   alllinks <- paste0("https://www.genecards.org/cgi-bin/carddisp.pl?gene=", alllink)
  #   inputdata <- inputdata %>%  mutate(`Gene Symbol` = paste0("<a href='", alllinks,"' target='_blank'>", alllink,"</a>"))
  #   
  #   DT::datatable(inputdata, escape = FALSE, rownames = FALSE,
  #                 options = list(columnDefs = list(list(className = 'dt-center', targets = "_all")))) %>% 
  #     formatSignif(columns = c("p-value"), digits = 3)
  # })
  
  
  
  
  
  # ----------------------------------------- get the input variables from the side bar input -------------------------------------------- #
  techtype <- reactive({
    input$techtype
  })
  
  tissuetype <- reactive({
   input$tissuetype
  })
  
  trt <- reactive({
    input$trt
  })
  
  ifmerge <- reactive({
    input$ifmerge
  })
  
  # change the trt input based on the techtype
  observeEvent(input$techtype, {
    
    if(input$techtype == "ATC"){
      updateSelectInput(session, "trt", choices = c( "Retinol", "Promatrixyl"))
    }
    if(input$techtype == "LCM"){
      updateSelectInput(session, "trt", choices = c(
        "Retinol", "Promatrixyl" ,"EGF+RF","Peptide","No Trt" ,"RP/SC800","EGF+RP+RF"))
    }
  })

  # get the probe list
  inputprobe <- eventReactive(c(input$ifallprobes, input$probe) , {

    if(input$ifallprobes){
     "all"
    }else{
      unlist(strsplit(x = input$probe, split = '[\r\n ]' ))
    }
  })
  
  
  # load the data for the analysis
  mydata <- reactive({
    inputprobe <- eventReactive(c(input$ifallprobes, input$probe) , {
      
      if(input$ifallprobes){
        "all"
      }else{
        unlist(strsplit(x = input$probe, split = '[\r\n ]' ))
      }
    })
    getdata(input$techtype, input$tissuetype, inputprobe(), input$trt, input$ifmerge)
  })
  
  
  output$downloadrawdata <- downloadHandler(
    
    
    filename =  function() {
      if(length(mydata()) > 1){
        paste('gene chip data', "_",Sys.Date(), "_",paste(names(mydata()), collapse = "_"),"_unmerged.xlsx", sep='')
      }else{
        paste('gene chip data', "_",Sys.Date(), "_",names(mydata()), ".xlsx", sep='')
      }
    },
    
    
    content = function(con) {
      shiny::withProgress(
        message = paste0("Downloading..."," less than 1 minute"),
        value = 0,
        {
          shiny::incProgress(1/20)
          Sys.sleep(2)
          shiny::incProgress(3/20)
          Sys.sleep(2)
          shiny::incProgress(5/20)
          Sys.sleep(2)
          shiny::incProgress(8/20)
          
          dataname = names(mydata())
          www <- list()
          for(i in 1:length(dataname)){
            www[[i]] = paste(dataname[i], names(mydata()[[i]]), sep = "_") #
          }
          
          outdata <- Reduce(append, mydata())
          
          newname <- unlist(www)
          names(outdata) <- newname
          
          write.xlsx(outdata, con, rowNames =TRUE)
          
        }
      )
    }
    
    
  )
  
  
  
  # ------------------------------------- #
  #         Start to do analyses  
  # ------------------------------------- #
  
  # ---------------------------------------- Panel 2. Differential abundance ---------------------------------------- #
  
  alltest <- eventReactive(input$run, {

    shiny::withProgress(
      message = paste0("Give me one second..."),
      value = 0,
      {
        shiny::incProgress(1/20)
        Sys.sleep(2)
        shiny::incProgress(5/20)
        Sys.sleep(2)
        shiny::incProgress(10/20)
        Sys.sleep(2)
        shiny::incProgress(15/20)

        testres <- performTest(mydata()[[1]][[1]], mydata()[[1]][[3]], method = input$da_method)
        testres[, c(2,4)] <- round(testres[,c(2,4)], 3)
        testres

      }
    )

  })
  
  
  # testing table
  output$alldata <- renderDataTable({
    
    
    if(!is.null(inputdata())){
        req(input$infile)
        
        inputdata <- inputdata()
        
        alllink <- sapply(strsplit(inputdata$`Gene Symbol`,split = "///"), `[`,1)
        alllinks <- paste0("https://www.genecards.org/cgi-bin/carddisp.pl?gene=", alllink)
        inputdata <- inputdata %>%  mutate(`Gene Symbol` = paste0("<a href='", alllinks,"' target='_blank'>", alllink,"</a>"))
        
        DT::datatable(inputdata, escape = FALSE, rownames = FALSE,
                      options = list(columnDefs = list(list(className = 'dt-center', targets = "_all")))) %>% 
          formatSignif(columns = c("p-value"), digits = 3)
    }else{
      
      req(input$trt)
      
      alltest <- alltest()
      #alltest[, c(3,5)] <- format(alltest[, c(3,5)], digits = 3,nsmall = 2, scientific = TRUE)
      
      alllink <- sapply(strsplit(alltest$`Gene Symbol`,split = "///"), `[`,1)
      alllinks <- paste0("https://www.genecards.org/cgi-bin/carddisp.pl?gene=", alllink)
      alltest <- alltest %>%  mutate(`Gene Symbol` = paste0("<a href='", alllinks,"' target='_blank'>", alllink,"</a>"))
      
      DT::datatable(alltest[,c(6,4,5,2,3)], escape = FALSE, rownames = FALSE,
                    options = list(columnDefs = list(list(className = 'dt-center', targets = "_all")))) %>%
        formatSignif(columns = c('p_TRT', 'p_FC'), digits = 3)
      
    }

  })
  
  # testing plot
  output$alldataplot = renderPlot({
    
    req(input$trt)
    
    s = input$alldata_rows_selected
    
    
    if(length(s) == 1){
      plottestdata(mydata()[[1]][[1]], mydata()[[1]][[3]], alltest(), s)
    }
    
  })
  
  selectgene <- reactive({
    
    req(input$trt)
    s = input$alldata_rows_selected
    
    if(length(s) == 1){
      alltest()[s, 6]
    }
  })
  
  boxplot <-  reactive({
    
    req(input$trt)
    s = input$alldata_rows_selected
    alltest()[s, 6]
    
    if(length(s) == 1){
      plottestdata(mydata()[[1]][[1]], mydata()[[1]][[3]], alltest(), s)
    }
    
  })
  
  
  output$downloadDatatest <- downloadHandler(
    
    
    filename =  function() {
      paste('Differential Abundance Testing Results', "_", trt(),"_",Sys.Date(), ".xlsx", sep='')
    },
    
    content = function(con) {
      
      alltest <- alltest()

      alllink  <- sapply(strsplit(alltest$`Gene Symbol`,split = "///"), `[`,1)
      alllinks <- paste0("https://www.genecards.org/cgi-bin/carddisp.pl?gene=", alllink)
      alltest$Links <- alllinks
      write.xlsx(alltest, con, rowNames =TRUE)          
      
    }
    
  )
  
  
  output$downloadboxplot <- downloadHandler(
    
    filename = function() { paste(selectgene(), '.png', sep='') },
    content = function(file) {
      ggsave(file, plot = boxplot(), device = "png")
    }
  )

  
  # ---------------------------------------- Panel 3. Functional Connection (must do DA analysis first) ---------------------------------------- #
  

  
  cut_fc <- reactive({
    c(input$effetcut_fc, input$pcut_fc)
  })
  
  cut_trt <- reactive({
    c(input$effetcut_trt, input$pcut_trt)
  })
  
  
  fun_gene_fc <- eventReactive(c(input$run, input$effetcut_fc, input$pcut_fc), {

    req(input$trt)
    getfunctiongene(alltest()[,c("probe", "effect_FC", "p_FC")], cut_fc())

  })
  
  fun_gene_trt <- eventReactive(c(input$run,input$effetcut_trt, input$pcut_trt), {
    
    if(is.null(inputdata())){
      req(input$trt)
      getfunctiongene(alltest()[,c("probe", "effect_TRT", "p_TRT")], cut_trt())
      
    }else{
      getfunctiongene(inputdata(), cut_trt())
    }
   
    
  })
  
  
  # volcano plot
  
  output$fc_volcano = renderPlot({
    
    req(input$trt)

    plotvolcano(alltest()[1:3], cut_fc()[1], cut_fc()[2])
    
  })
  
  output$trt_volcano = renderPlot({
    
    
    if(is.null(inputdata())){
      req(input$trt)
      
      plotvolcano(alltest()[c(1,4,5)], cut_trt()[1], cut_trt()[2])
    }else{
      plotvolcano(inputdata(), cut_trt()[1], cut_trt()[2])
      
    }
    
  })
  
  
  
  
  
  fun_res <- eventReactive(input$runfun,{

    shiny::withProgress(
      message = paste0("Give me one second..."),
      value = 0,
      {
        shiny::incProgress(1/20)
        Sys.sleep(2)
        shiny::incProgress(5/20)
        Sys.sleep(2)
        shiny::incProgress(10/20)
        Sys.sleep(2)
        shiny::incProgress(15/20)

        
        if(is.null(inputdata())){
          trt_res <- getfunction(fun_gene_trt(), title = "Based on the treatment")
          fc_res <- getfunction(fun_gene_fc(), title = "Based on the fold-change")
          list(trt_res, fc_res)
          
        }else{
          getfunction(fun_gene_trt(),title = "Based on the treatment")
        }
        

      }
    )
  })
  
  # output gene table for functional analysis
  output$fc_table_gene <- renderDataTable({
    req(input$trt)
    DT::datatable(data.frame(`Selected Gene` = fun_gene_fc()), rownames = FALSE,
                  options = list(columnDefs = list(list(className = 'dt-center', targets = "_all"))))
  })
  output$trt_table_gene <- renderDataTable({
    
    if(is.null(inputdata())){
      req(input$trt)
      DT::datatable(data.frame(`Selected Gene` = fun_gene_trt()), rownames = FALSE,
                    options = list(columnDefs = list(list(className = 'dt-center', targets = "_all"))))
    }else{
      
      DT::datatable(data.frame(`Selected Gene` = fun_gene_trt()), rownames = FALSE,
                    options = list(columnDefs = list(list(className = 'dt-center', targets = "_all"))))
    }
   
  })

  
  # function plot
  output$trt_plot = renderPlot({
    
    if(is.null(inputdata())){
    req(input$trt)
    fun_res()[[1]]$plot
    }else{
      fun_res()$plot
    }
  })
  
  output$fc_plot = renderPlot({
    
    req(input$trt)
    fun_res()[[2]]$plot
  })
  
  
  output$fc_table <- renderDataTable({
    req(input$trt)
    mytable2 <- fun_res()[[2]]$table[,c(1,2,3,5,7)]
    mytable2[,c(4,5)] <- format(mytable2[, c(4,5)], digits = 3,nsmall = 2, scientific = TRUE)
    DT::datatable(mytable2, rownames = FALSE,
                  options = list(columnDefs = list(list(className = 'dt-center', targets = "_all"))))
  })
  output$trt_table <- renderDataTable({
    
    if(is.null(inputdata())){
    req(input$trt)
    mytable3 <- fun_res()[[2]]$table[,c(1,2,3,5,7)]
    mytable3[,c(4,5)] <- format(mytable3[, c(4,5)], digits = 3,nsmall = 2, scientific = TRUE)
    DT::datatable(mytable3, rownames = FALSE,
                  options = list(columnDefs = list(list(className = 'dt-center', targets = "_all"))))
    }else{
      mytable3 <- fun_res()$table[,c(1,2,3,5,7)]
      mytable3[,c(4,5)] <- format(mytable3[, c(4,5)], digits = 3,nsmall = 2, scientific = TRUE)
      DT::datatable(mytable3, rownames = FALSE,
                    options = list(columnDefs = list(list(className = 'dt-center', targets = "_all")))) 
    }
  })
  
  
  
  
  
  # ----------------------------------------- Panel 4. Material Connection  ----------------------------------------- #
  


  cmap_res <- eventReactive(c(input$runcmap, input$cmapgene), {


    shiny::withProgress(
      message = paste0("Give me one second..."),
      value = 0,
      {
        shiny::incProgress(1/20)
        Sys.sleep(2)
        shiny::incProgress(5/20)
        Sys.sleep(2)
        shiny::incProgress(10/20)
        Sys.sleep(2)
        shiny::incProgress(15/20)

        if(is.null(inputdata())){
          req(input$trt)
          pheno <- mydata()[[1]][[3]]$group
          fc <- rowMeans(mydata()[[1]][[1]][, pheno == "Reponse"])
          getcmap(fc, input$cmapgene)
        }else{
          myfc <- inputdata()[,2] %>% pull
          names(myfc) <- inputdata()[,1] %>% pull
          getcmap(myfc, input$cmapgene)
        }
        

      }
    )

  })
  
  
  
  
  output$camptable <- renderDataTable({
    
    if(is.null(inputdata())){
      req(input$trt)
    }
    
    
    myt1     <- cmap_res()[[1]][,c(1,2,3,6)]
    myt1$Rank_GSEA <- base::rank(myt1[,4]*-1, ties.method = "min")
    myt1[,4] <- round(myt1[,4], 3)
    colnames(myt1)[4] <- "Score_GSEA"
    
    myt2     <- cmap_res()[[2]][,c(1,4)]
    myt2$Rank_Sim <- base::rank(myt2[,2]*-1, ties.method = "min")
    myt2[,2] <- round(myt2[,2], 3)
    colnames(myt2)[2] <- "Score_Sim"
    
    myt12 <- myt1 %>% left_join(myt2)
    
    DT::datatable(myt12, rownames = FALSE,
                  options = list(columnDefs = list(list(className = 'dt-center', targets = "_all")),
                                 pageLength = 15))
  })
  
  output$downloadDatamc <- downloadHandler(
    
    
    filename =  function() {
      paste('Cmap',"_",trt(),"_",Sys.Date(), ".xlsx", sep='')
    },
    
    content = function(con) {

      write.xlsx(cmap_res()[[1]], con, rowNames =TRUE)          
      
    }
    
  )
  

  
  
  # ---------------------------------------- Panel 5. Underlying Pattern (unsupervised learning)  ---------------------------------------- #
  
  unsuper_res <- eventReactive(input$rununsuper, {
    
    
    shiny::withProgress(
      message = paste0("Give me one second..."),
      value = 0,
      {
        shiny::incProgress(1/20)
        Sys.sleep(2)
        shiny::incProgress(5/20)
        Sys.sleep(2)
        shiny::incProgress(10/20)
        Sys.sleep(2)
        shiny::incProgress(15/20)
        
        if(input$datatype == "Rank"){
          unsuperdata_fc <- apply(mydata()[[1]][[1]], 2, function(x)rank(-x))
          unsuperdata_bl <- apply(mydata()[[1]][[2]], 2, function(x)rank(-x))
          
        }else{
          unsuperdata_fc <- mydata()[[1]][[1]]
          unsuperdata_bl <- mydata()[[1]][[2]]
        }
        
        metaunsuper <- mydata()[[1]][[3]][,c(1:3)]
        
        if(!is.null(inputdata())){



          myfc <- inputdata()[,2] %>% pull
          names(myfc) <- inputdata()[,1] %>% pull

          name1 <- rownames(unsuperdata_fc)
          name2 <- names(myfc)

          cname <- intersect(name1, name2)

          
          unsuperdata_fc <- cbind(unsuperdata_fc[cname,], newdata = myfc[cname])
          unsuperdata_bl <- cbind(unsuperdata_bl[cname,], newdata = myfc[cname])

          
          metaunsuper <- rbind(metaunsuper, data.frame(name = "newdata", thickdiff = 0, group = "unknow"))
          

        }
        
        upsuper_fc <- getunsupervised(unsuperdata_fc, metaunsuper, "Based on Fold Change")
        upsuper_bl <- getunsupervised(unsuperdata_bl, metaunsuper, "Based on Vehicle Expression")
        
        list(upsuper_fc, upsuper_bl)
      }
    )
    
  })
  
  
  
  # unsupervised learning plot
  output$pcaplot_fc = renderPlotly({
    req(input$trt)
    unsuper_res()[[1]]$pca
  })
  output$tsneplot_fc = renderPlotly({
    req(input$trt)
    unsuper_res()[[1]]$tsne
  })
  output$umapplot_fc = renderPlotly({
    req(input$trt)
    unsuper_res()[[1]]$umap
  })
  output$hcplot_fc = renderPlot({
    
    req(input$trt)

    
    if(input$datatype == "Rank"){
      unsuperdata_fc <- apply(mydata()[[1]][[1]], 2, function(x)rank(-x))

    }else{
      unsuperdata_fc <- mydata()[[1]][[1]]
    }
    
    metaunsuper <- mydata()[[1]][[3]]
    gethclust(unsuperdata_fc, metaunsuper, "Based on Fold Change")

    
  })
  
  output$pcaplot_bl = renderPlotly({
    req(input$trt)
    unsuper_res()[[2]]$pca
  })
  output$tsneplot_bl = renderPlotly({
    req(input$trt)
    unsuper_res()[[2]]$tsne
  })
  output$umapplot_bl = renderPlotly({
    req(input$trt)
    unsuper_res()[[2]]$umap
  })
  output$hcplot_bl = renderPlot({
    req(input$trt)
    #unsuper_res()[[2]]$hc
    
    if(input$datatype == "Rank"){
      unsuperdata_bl <- apply(mydata()[[1]][[2]], 2, function(x)rank(-x))
      
    }else{
      unsuperdata_bl <- mydata()[[1]][[2]]
    }
    
    metaunsuper <- mydata()[[1]][[3]]
    getunsupervised(unsuperdata_bl, metaunsuper, "Based on Vehicle Expression")
  })
  
  
  # ---------------------------------------- Panel 6. Response Prediction (supervised learning)  ---------------------------------------- #
  
 
  pred_res <- eventReactive(input$runpred, {
    
    
    shiny::withProgress(
      message = paste0("Give me one second..."),
      value = 0,
      {
        shiny::incProgress(1/20)
        Sys.sleep(2)
        shiny::incProgress(5/20)
        Sys.sleep(2)
        shiny::incProgress(10/20)
        Sys.sleep(2)
        shiny::incProgress(15/20)
        
        if(input$datatype == "Rank"){
          unsuperdata_fc <- apply(mydata()[[1]][[1]], 2, function(x)rank(-x))
          #unsuperdata_bl <- apply(mydata()[[1]][[2]], 2, function(x)rank(-x))
          
        }else{
          unsuperdata_fc <- mydata()[[1]][[1]]
          #unsuperdata_bl <- mydata()[[1]][[2]]
        }
        
        if(input$numbergene == 49293){
          testtable <- NULL
          topnumber <- 49293
        }else{
          testtable <- alltest()
          topnumber <- input$numbergene
        }
        
        pheno <- mydata()[[1]][[3]]$group
        
        getprediction(t(unsuperdata_fc), pheno, testtable, topnumber)

      }
    )
    
  })
  
  
  output$roc = renderPlot({
    req(input$trt)
    
    getroc(pred_res()[[1]],  mydata()[[1]][[3]]$group)
  })
  
  output$predtable = renderDataTable({
    req(input$trt)
    
    predtable <- getpredtable(pred_res()[[1]], mydata()[[1]][[3]])
    predtable[,2] <- round(predtable[,2],2)
    DT::datatable(predtable[,c(1:3,11)], rownames = FALSE,
                  options = list(columnDefs = list(list(className = 'dt-center', targets = "_all")),
                                 pageLength = 15))
  })
  
    
}

# Run the application 
shinyApp(ui = ui, server = server)

