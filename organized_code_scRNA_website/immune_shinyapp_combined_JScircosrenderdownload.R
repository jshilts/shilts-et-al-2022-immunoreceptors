library(shiny)
source("./shinyCircos.R")
library(tidyverse)

ui <- fluidPage(
  
  titlePanel("Interaction plots"),
  
  # Sidebar layout ----
  fluidRow( 

    conditionalPanel(
      condition = 'input.tabs !== "Circos view"',
        column(4, offset = 1,

           sliderInput(inputId = "thresh",
                       label = "Threshold for expression:",
                       min = 1,
                       max = 50,
                       value = 10))
      ),

    conditionalPanel(
      condition = 'input.tabs == "Circos view"',
      column(width = 2, offset = 1,
             
             sliderInput(inputId = "thresh",
                         label = "Threshold for expression:",
                         min = 1,
                         max = 50,
                         value = 10)
      ),
      column(width = 2, offset = 1,
             sliderInput(inputId ='pctmaxshare',
                  label = div(style = "font-size:10px", "Filter out if interaction shared by more than this percent of cell pairs"),
                  min = 1,
                  max = 100,
                  value = 20))
    ),
    
    column(3, offset = 1, 
           selectInput("dataset", "Human tissue :",
                       c("Spleen" = "spleen_cmin10_intxn_leukmanual", 
                         "Thymus" = "thymus_cmin10_intxn_leukmanual",
                         "Colon" = "colon_Healthy_cmin10_intxn_leukmanual",
                         "Lung" = "lung_Healthy_cmin10_intxn_leukmanual",
                         "Kidney" = "kidney_Healthy_cmin10_intxn_leukmanual",
                         "Bone Marrow" = "bone_cmin10_intxn_leukmanual",
                         "Lymph Node" = "lymph_cmin10_intxn_leukmanual",
                         "Blood" = "pbmc_cmin10_intxn_leukmanual"
                       ))
    ),
    column(1, offset = 0, style = "padding:0px;",
           HTML("<div style='height: 100px;'>"),
           imageOutput(outputId = "iconImg")),
    HTML("</div>")
  ),

  # Main panel for displaying outputs ----
   fluidRow( 
    tabsetPanel(type = "tabs", id = "tabs",
                tabPanel("Sum view", plotOutput("sumplot"), downloadLink("downloadSumplot", "Save plot image")),
                tabPanel("Cell view", plotOutput("gridplot"), downloadLink("downloadGridplot", "Save plot image")),
                tabPanel("Protein view", plotOutput("quiltplot"), downloadLink("downloadQuiltplot", "Save plot image")),
                tabPanel("Circos view", plotOutput("circosplot"), downloadLink("downloadCircosplot", "Save plot image"))
    )
    
  )
  
)




server <- function(input, output) {
  
  output$iconImg <- renderImage({
    tissue_name = input$dataset
    tissue_name = gsub("_.*$", "", tissue_name)
    
    list(src = normalizePath(paste0(tissue_name, ".png")),
         width = 50,
         height = 50)
    
  }, deleteFile = FALSE)
  
  
  data_key_gen <- reactive({  # speeds up by only making reload data when data's input name changes (instead of every time slider moves)
    fname = input$dataset
    
    if_quant_datakey = TRUE
    if (if_quant_datakey) {
      read_csv(paste0("./circos_", fname, "_QKEY.csv"))
      
    }
  })
  
  
  makeSumPlot = function() {
    data_key = data_key_gen()
    
    if_single_condition = TRUE
    
    data_key$Expected = 1
    data_key = mutate(data_key, pair_cells = ifelse(chr1 < chr2, paste0(chr1, " + ", chr2), paste0(chr2, " + ", chr1))) 
    data_key = mutate(data_key, bin_condition = Expected)
    
    user_threshold = input$thresh / 100
    combined_data_key = filter(data_key, pct_expr_1 > user_threshold & pct_expr_2 > user_threshold)
    
    data_key_summary = combined_data_key %>% 
      group_by(chr1, chr2) %>%
      summarize(n_expect = sum(Expected==1), n_draft = sum(Expected==0)) %>%
      ungroup() %>% 
      complete(chr1, chr2, fill = list("n_intxns" = 0))  # puts zero where see no intxns
    
    ggplot(data_key_summary, aes(x = chr1, y = chr2)) + 
      geom_tile(aes(fill = n_expect)) +
      theme_bw() +
      coord_fixed(ratio = 1) +
      theme(axis.text.x =  element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
            axis.text.y =  element_text(vjust = 0.5, size = 8)) +
      xlab("") + ylab("") +
      scale_fill_viridis_c(option = "E", name = "# inxns", na.value = "black")
  }
  output$sumplot <- renderPlot({
    makeSumPlot()
  })
  output$downloadSumplot = downloadHandler(
    filename = function() { paste0(gsub("_.*$", "", input$dataset), "_SumPlot_p", input$thresh, '.png') },
    content = function(file) {
      ggsave(file, plot = makeSumPlot(), width = 7, height = 7)
    }
  )
  
  makeGridPlot = function() {
    
    if_quant_datakey = TRUE
    if_single_condition = TRUE
    
    data_key = data_key_gen()
    
    data_key$Expected = 1
    data_key = mutate(data_key, pair_cells = ifelse(chr1 < chr2, paste0(chr1, " + ", chr2), paste0(chr2, " + ", chr1))) 
    data_key = mutate(data_key, bin_condition = Expected)
    
    user_threshold = input$thresh / 100
    combined_data_key = filter(data_key, pct_expr_1 > user_threshold & pct_expr_2 > user_threshold)
    
    hclust_intxns_tbl = combined_data_key %>%
      complete(pair_cells, pair_label) %>%
      mutate(bin_condition = ifelse(is.na(bin_condition), 0, bin_condition)) %>%
      select(pair_cells, pair_label, bin_condition) %>%
      distinct() %>%
      spread(key = pair_cells, value = bin_condition)
    
    hclust_intxns_m = as.matrix(select(hclust_intxns_tbl, -pair_label))
    rownames(hclust_intxns_m) = hclust_intxns_tbl$pair_label
    hclust_intxns_clust = hclust(dist(hclust_intxns_m))
    combined_data_key$pair_label = factor(combined_data_key$pair_label, levels = hclust_intxns_clust$labels[hclust_intxns_clust$order])
    
    ggplot(combined_data_key, aes(x=chr1, y = pair_label, chr2)) + 
      {if(if_quant_datakey)
        geom_tile(aes(fill=log2(pct_min)))} + 
      {if(!if_quant_datakey)
        geom_tile(aes(fill=factor(bin_condition)))} + 
      facet_wrap(~chr2, nrow=1) +
      theme_bw() +
      theme(axis.text.x =  element_text(size = 3, angle = 90, vjust = 0.5),
            axis.text.y =  element_text(size = 4.5),
            strip.text.x = element_text(size = 5)) + 
      xlab("") + ylab("Interaction") +
      if(if_single_condition){
        if(if_quant_datakey){
          scale_fill_viridis_c(breaks = c(log2(.99), log2(0.75), log2(0.5), log2(0.25), log2(0)),
                               labels = c("100%", "75%", "50%", "25%", "0%"), name ="Express %", na.value = '#440154FF')  # min from viridis::viridis(n = 2)
        } else{
          scale_fill_manual(values = "#c0bfe0", name = "Interaction")
        }
      } else {
        scale_fill_manual(breaks = c(1, 2, 3), labels = c("Asthma", "Control", "Both"), name = "Lung condition",
                          values = c("1" = "#e3c0b1", "2" = "#b1d6c6", "3" = "#c0bfe0"))
      }
  }
  
  output$gridplot <- renderPlot({
    makeGridPlot()
  })
  output$downloadGridplot = downloadHandler(
    filename = function() {paste0(gsub("_.*$", "", input$dataset), "_CellPlot_p", input$thresh, '.png')},
    content = function(file) {
      ggsave(file, plot = makeGridPlot(), width = 18, height = 5)
    }
  )
  
  
  makeQuiltPlot = function() {
    
    data_key = data_key_gen()
    
    if_single_condition = TRUE
    if_quant_datakey = TRUE
    
    data_key$Expected = 1
    data_key = mutate(data_key, pair_cells = ifelse(chr1 < chr2, paste0(chr1, " + ", chr2), paste0(chr2, " + ", chr1)))
    data_key = mutate(data_key, bin_condition = Expected)
    
    user_threshold = input$thresh / 100
    combined_data_key = filter(data_key, pct_expr_1 > user_threshold & pct_expr_2 > user_threshold)
    
    hclust_intxns_tbl = combined_data_key %>%
      complete(pair_cells, pair_label) %>%
      mutate(bin_condition = ifelse(is.na(bin_condition), 0, bin_condition)) %>%
      select(pair_cells, pair_label, bin_condition) %>%
      distinct() %>%
      spread(key = pair_cells, value = bin_condition)
    
    hclust_intxns_m = as.matrix(select(hclust_intxns_tbl, -pair_label))
    rownames(hclust_intxns_m) = hclust_intxns_tbl$pair_label
    hclust_intxns_clust = hclust(dist(hclust_intxns_m)) 
    combined_data_key$pair_label = factor(combined_data_key$pair_label, levels = hclust_intxns_clust$labels[hclust_intxns_clust$order])
    
    if_expbin_coloring = FALSE
    ggplot({if(if_expbin_coloring){combined_data_key_gr1}else{combined_data_key}}, 
           aes(x=chr1, y = chr2)) + 
      geom_tile(aes(fill= log2(pct_min))) +
      facet_wrap(~pair_label, drop = FALSE, nrow = max(3, round(length(unique(combined_data_key$pair_label))/15))) +
      theme_bw() +
      coord_fixed(ratio = 1) +
      theme(axis.text.x =  element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5),
            axis.text.y =  element_text(vjust = 0.5, size = 5),
            strip.text.x = element_text(size = 5)) +
      xlab("") + ylab("") +
      if(if_single_condition){
        if(if_quant_datakey){
          if(if_expbin_coloring){
            scale_fill_manual(values = c("#FFFFFF", RColorBrewer::brewer.pal(9, 'OrRd'), "#54050e"), name = "Percent expressed")
          } else{
            scale_fill_viridis_c(breaks = c(log2(.99), log2(0.75), log2(0.5), log2(0.25), log2(0)),
                                 labels = c("100%", "75%", "50%", "25%", "0%"), name ="Express %", na.value = '#440154FF')  # min from viridis::viridis(n = 2)
          }
        } else{
          scale_fill_manual(values = "#c0bfe0", name = "Interaction")
        }
      } else {
        scale_fill_manual(breaks = c(1, 2, 3), labels = c("Asthma", "Control", "Both"), name = "Lung condition",
                          values = c("1" = "#e3c0b1", "2" = "#b1d6c6", "3" = "#c0bfe0"))
      }
  }
  
  output$quiltplot <- renderPlot({
    makeQuiltPlot()
  })
  output$downloadQuiltplot = downloadHandler(
    filename = function() { paste0(gsub("_.*$", "", input$dataset), "_ProtPlot_p", input$thresh, '.png') },
    content = function(file) {
      ggsave(file, plot = makeQuiltPlot(), width = 15, height = 10)
    }
  )
  
  makeCircosPlot = function() {
    circos.clear()
    data_key = data_key_gen()
    
    if_single_condition = TRUE
    
    data_key$Expected = 1
    data_key = mutate(data_key, pair_cells = ifelse(chr1 < chr2, paste0(chr1, " + ", chr2), paste0(chr2, " + ", chr1)))
    data_key = mutate(data_key, bin_condition = Expected)
    
    user_threshold = input$thresh / 100
    user_pctmaxshare = input$pctmaxshare / 100
    combined_data_key = filter(data_key, pct_expr_1 > user_threshold & pct_expr_2 > user_threshold)
    combined_data_key %<>%
      group_by(pair_label) %>%
      mutate(sum_pair = n(),
             pct_max = sum_pair / length(unique(c(combined_data_key$chr1, combined_data_key$chr2)))**2) %>%
      ungroup() %>%
      filter(pct_max < user_pctmaxshare)
    
    
    data.C <- data.frame(chr = unique(c(combined_data_key$chr1, combined_data_key$chr2)), start = 1, end = max(c(combined_data_key$end1, combined_data_key$end2)))
    data.C[,2] <- as.numeric(data.C[,2])
    data.C[,3] <- as.numeric(data.C[,3])
    data.T.file <- c("")
    data.CN <- NULL
    data.N.file <- c("","","","","","","","","","")
    uploadtrack <- c(1,1,1,1,1,1,1,1,1,1)
    data.N <- lapply(1:10,function(x){
      if(uploadtrack[x] == 2 && nchar(data.N.file[x])>0){	  
        data.frame(fread(data.N.file[x]),stringsAsFactors=F)
      }
    })
    data.T <- NULL
    trackindx <- c()
    data.N <- data.N[trackindx]
    data.N <- NULL
    
    data.L <- data.frame(select(combined_data_key, chr1, start1, end1, chr2, start2, end2, color))
    data.L1 <- data.L[,1:3]
    data.L2 <- data.L[,4:6]
    data.L1[,2] <- as.numeric(data.L1[,2])
    data.L1[,3] <- as.numeric(data.L1[,3])
    data.L2[,2] <- as.numeric(data.L2[,2])
    data.L2[,3] <- as.numeric(data.L2[,3])	  
    data.L1$num <- 1:nrow(data.L1)
    data.L2$num <- 1:nrow(data.L2)
    rownames(data.L1) <- data.L1$num
    rownames(data.L2) <- data.L2$num
    for(i in 1:length(data.T.file)){
      assign(paste("hltdata",i,sep=""),"")
    }
    hltregion.List <- list()
    if(!is.null(data.T)){
      for(k in 1:length(data.T)){
        data.TT <- data.T[[k]]
        hltregion.List[[k]] <- ""
        if(nchar(get(paste("hltdata",k,sep="")))>0){
          tmp <- matrix(strsplit(get(paste("hltdata",k,sep="")), "\n")[[1]])
          myColnames <- c("chr","start","end","color")
          data <- matrix(0, length(tmp), length(myColnames))
          colnames(data) <- myColnames
          for(p in 1:length(tmp)){
            myRow <- strsplit(tmp[p], ",")[[1]]
            if(length(myRow)==4){                                        
              data[p,] <- myRow
            }
          }
          data <- data.frame(data,stringsAsFactors=F)
          data$start <- as.numeric(data$start)
          data$end <- as.numeric(data$end)
          query <- GRanges(seqnames = data$chr,ranges=IRanges(start=data$start,end=data$end),seqinfo=NULL)
          subj <- GRanges(seqnames = data.TT[,1],ranges=IRanges(start=data.TT[,2],end=data.TT[,3]),seqinfo=NULL) 
          indx <- findOverlaps(query,subj)
          indx <- data.frame(indx,stringsAsFactors=F)
          indx$queryHits <- as.numeric(indx$queryHits)
          indx$subjectHits <- as.numeric(indx$subjectHits)
          hltregion <- data.TT[indx$subjectHits,]
          hltregion$color <- data$color[indx[,1]]
          hltregion$id <- paste(hltregion[,1],hltregion[,2],hltregion[,3],sep="")
          hltregion.List[[k]] <- hltregion
        }
      }
    }
    fontSize <- 1
    par(mar=c(0.6, 0.6, 0.6, 0.6), cex=fontSize-0.05)
    trackChr <- "track"
    labelChr <- "labels"
    unitChr <- ""
    rotation <- 0.5
    gap.width <- rep(1, times=nrow(data.C))
    labeltextchr <- 2
    poslabelschr <- "inner"
    heightlabelschr <- 0.06
    marginlabelschr <- 0.01
    colorChr <- rep("grey", times=nrow(data.C))
    
    
    plotcircos(data.C, color=colorChr, plotTypes=unique(c(labelChr,"axis")), units=unitChr, rotation=rotation, gap.width=gap.width, labeltextchr=labeltextchr, poslabelschr=poslabelschr, heightlabelschr=heightlabelschr, marginlabelschr=marginlabelschr, data.CN=data.CN)
    marginLinks <- 0.01  # margin for the gap between lines and outer circle 
    circos.par(track.margin = c(0,marginLinks))
    transparencyLinks <- 0.9  # transparency of connecting link lines
    rou <- get_most_inside_radius()
    rou <- rou[1]
    linkscolor.export = unique(data.L$color)
    groupnum <- length(unique(data.L[,7]))
    randcolorLinks <- data.frame(group=unique(data.L[,7]), cols=linkscolor.export, stringsAsFactors=F)
    data.LC <- merge(data.L,randcolorLinks,by.x="color",by.y="group",all.x=T)	
    colLinks <- adjustcolor(data.LC$cols, alpha.f = transparencyLinks)
    data.L1 <- data.L1[,c(1:3)]
    data.L2 <- data.L2[,c(1:3)]
    circos.genomicLink(data.L1, data.L2, rou = rou, col = colLinks, border = colLinks,
                       lwd = 1)
  }
  
  output$circosplot <- renderPlot({
    makeCircosPlot()
  }, res = 100)  # image resolution 
  
  output$downloadCircosplot = downloadHandler(
    filename = function() { paste0(gsub("_.*$", "", input$dataset), "_CircosPlot_p", input$thresh, "_f", input$pctmaxshare, '.pdf') },
    content = function(file) {
      ggsave(file, plot = makeCircosPlot(), width = 10, height = 10)
    }
  )
}
  

shinyApp(ui = ui, server = server)

