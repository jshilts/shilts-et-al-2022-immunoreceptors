library(org.Hs.eg.db)
library(tidyverse)
library(Seurat)
library(nichenetr)  # devtools::install_github("saeyslab/nichenetr")

setwd("./organized_code_scRNA_analysis")

#### reading data files ####

lr_network = readRDS("./nichenet_internal/lr_network.rds") #readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
weighted_networks = readRDS("./nichenet_internal/weighted_networks.rds")  #readRDS(url("https://zenodo.org/record/3260758/files/weighted_networks.rds"))

sdata = readRDS("sdata_kidney.RDS")


#### building custom receptor network ####

lr_valid_names = unique(c(lr_network$from, lr_network$to))

leuk_intxn_set = read_csv("../organized_code_affinity/Leuk_interaction_affinities.csv")  
leuk_lr_full = leuk_intxn_set %>%
  mutate(from = ifelse(Gene_L %in% lr_valid_names, Gene_L, "DOENTREZCONVERT"),
         to = ifelse(Gene_R %in% lr_valid_names, Gene_R, "DOENTREZCONVERT")) %>%
  mutate(from_entrez = unname(AnnotationDbi::mapIds(org.Hs.eg.db, keys = .$Uniprot_L, column = "ENTREZID", keytype = "UNIPROT", multiVals = "first")),
         to_entrez = unname(AnnotationDbi::mapIds(org.Hs.eg.db, keys = .$Uniprot_R, column = "ENTREZID", keytype = "UNIPROT", multiVals = "first"))) %>%
  mutate(from = ifelse(from == "DOENTREZCONVERT", plyr::mapvalues(from_entrez, geneinfo_human$entrez, geneinfo_human$symbol, warn_missing = F), from),
         to = ifelse(to == "DOENTREZCONVERT", plyr::mapvalues(to_entrez, geneinfo_human$entrez, geneinfo_human$symbol, warn_missing = F), to))
#write_csv(leuk_lr_full, "leuk_lr_full.csv")   # can manually fill in a few gene IDs that this autoconvert tool missed

lr_network_leuk = read_csv("./nichenet_internal/Leuk_custom_LR_list.csv")

ligand_target_matrix_leuk = construct_ligand_target_matrix(weighted_networks = weighted_networks, 
                                                           ligands = as.list(unique(c(lr_network_leuk$from, lr_network_leuk$to))),   # needs to be in format where each item in separate sublist of length 1 for normal ligands, or length 2 for heterodimer ligands
                                                           algorithm = "PPR", damping_factor = hyperparameter_list$damping_factor, ltf_cutoff = hyperparameter_list$ltf_cutoff)
#saveRDS(ligand_target_matrix_leuk, "/nichenet_internal/ligand_target_matrix_leuk.rds")
ligand_target_matrix_leuk = readRDS("./nichenet_internal/ligand_target_matrix_leuk.rds")

pct_expr_threshold = 0.10
logfc_threshold = 0.25 
if_any_sender = FALSE
## receiver
receiver = "Th cell"
expressed_genes_receiver = get_expressed_genes(receiver, sdata, pct = pct_expr_threshold)
background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix_leuk)]
## sender
sender_celltypes = c("Plasmacytoid DC")
list_expressed_genes_sender = sender_celltypes %>% unique() %>% lapply(get_expressed_genes, sdata, pct_expr_threshold)
expressed_genes_sender = list_expressed_genes_sender %>% unlist() %>% unique()


nichenet_tumor_analysis = function(sdata, receiver, sender_celltypes,
                                   expressed_genes_receiver, background_expressed_genes,
                                   list_expressed_genes_sender, expressed_genes_sender,
                                   pct_expr_threshold, logfc_threshold, 
                                   if_any_sender){
  ## define gene set
  seurat_obj_receiver = subset(sdata, idents = receiver)
  seurat_obj_receiver = SetIdent(seurat_obj_receiver, value = seurat_obj_receiver[["Condition"]])
  condition_oi = "Tumor"
  condition_reference = "Healthy" 
  DE_table_receiver = FindMarkers(object = seurat_obj_receiver, ident.1 = condition_oi, ident.2 = condition_reference, min.pct = pct_expr_threshold, logfc.threshold = logfc_threshold, test.use = "wilcox") %>% rownames_to_column("gene")
  geneset_oi = DE_table_receiver %>% filter(p_val_adj <= 0.05 & abs(avg_logFC) >= logfc_threshold) %>% pull(gene)
  geneset_oi = geneset_oi %>% .[. %in% rownames(ligand_target_matrix_leuk)]
  ## define ligand set
  ligands = lr_network_leuk %>% pull(from) %>% unique()
  receptors = lr_network_leuk %>% pull(to) %>% unique()
  expressed_ligands = intersect(ligands,expressed_genes_sender)
  expressed_receptors = intersect(receptors,expressed_genes_receiver)
  if (if_any_sender) {
    potential_ligands = ligands
  } else {
    potential_ligands = lr_network_leuk %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% pull(from) %>% unique()  
  }
  ## ligand activity scores
  ligand_activities = predict_ligand_activities(geneset = geneset_oi, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix_leuk, potential_ligands = potential_ligands)
  ligand_activities = ligand_activities %>% arrange(-pearson) %>% mutate(rank = rank(desc(pearson)))
  best_upstream_ligands = ligand_activities %>% top_n(20, pearson) %>% arrange(-pearson) %>% pull(test_ligand) %>% unique()
  
  p_all_senders = DotPlot(sdata, features = best_upstream_ligands %>% rev(), cols = "YlGn") + 
    scale_color_gradient(low = "#edeada", high = "#346900") + RotatedAxis() + xlab("ligands") + ggtitle("Which cell types top-ranked ligands most found in")
  
  ## targets of receptor signaling
  active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links, geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix_leuk, n = 200) %>% bind_rows() %>% drop_na()
  active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix_leuk, cutoff = 0.33)
  order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev() %>% make.names()
  order_targets = active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links)) %>% make.names()
  rownames(active_ligand_target_links) = rownames(active_ligand_target_links) %>% make.names()
  colnames(active_ligand_target_links) = colnames(active_ligand_target_links) %>% make.names()
  vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()
  
  lr_network_top = lr_network_leuk %>% filter(from %in% best_upstream_ligands & to %in% expressed_receptors) %>% distinct(from,to)
  if (!if_any_sender) {
    rownames(vis_ligand_target) = unname(sapply(rownames(vis_ligand_target), (function(x) paste0(x, " + ", filter(lr_network_top, from == x)$to))))
    vis_ligand_target = vis_ligand_target[apply(vis_ligand_target, 1, function(x){sum(x>0)}) < ncol(vis_ligand_target)/2, ]   # filter out non-specific/ubiquitous predictions for visualization
    vis_ligand_target = vis_ligand_target[ , apply(vis_ligand_target, 2, function(x){sum(x>0)}) < nrow(vis_ligand_target)-1]      
    vis_ligand_target = vis_ligand_target[ , apply(vis_ligand_target, 2, function(x){sum(x)!=0})]  # remove empty columns
    vis_ligand_target = vis_ligand_target[apply(vis_ligand_target, 1, function(x){sum(x)!=0}), ]  # remove empty rows
  }
  p_ligand_target_network = vis_ligand_target %>% make_heatmap_ggplot("Prioritized interactions","Predicted downstream target genes", color = "purple",legend_position = "bottom", x_axis_position = "top",legend_title = "Regulatory potential")  + 
    theme(axis.text.x = element_text(face = "italic")) + 
    {if(if_any_sender) scale_fill_gradient2(low = "whitesmoke",  high = "purple", limits = c(0,max(vis_ligand_target)*0.35), oob = scales::squish)} +  # optimize color range based on spread of data 
    {if(!if_any_sender) scale_fill_gradient2(low = "whitesmoke",  high = "purple", limits = c(0,max(vis_ligand_target)*1.5), oob = scales::squish)}  
  included_prot = gsub(" ", "", unlist(strsplit(rownames(vis_ligand_target), "\\+"))) 
  included_lig = gsub(" \\+.*$", "", rownames(vis_ligand_target))
  ## order the corresponding receptors for the top ligands
  best_upstream_receptors = lr_network_top %>% pull(to) %>% unique()
  
  weighted_networks_lr = weighted_networks$lr_sig %>% inner_join(lr_network_leuk %>% distinct(from,to), by = c("from","to"))
  lr_network_top_df_large = weighted_networks_lr %>% filter(from %in% best_upstream_ligands & to %in% best_upstream_receptors)
  lr_network_top_df = lr_network_top_df_large %>% spread("from","weight",fill = 0)
  lr_network_top_matrix = lr_network_top_df %>% select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df$to)
  
  dist_receptors = dist(lr_network_top_matrix, method = "binary")
  hclust_receptors = hclust(dist_receptors, method = "ward.D2")
  order_receptors = hclust_receptors$labels[hclust_receptors$order]
  dist_ligands = dist(lr_network_top_matrix %>% t(), method = "binary")
  hclust_ligands = hclust(dist_ligands, method = "ward.D2")
  order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]
  order_receptors = order_receptors %>% intersect(rownames(lr_network_top_matrix))
  order_ligands_receptor = order_ligands_receptor %>% intersect(colnames(lr_network_top_matrix))
  
  vis_ligand_receptor_network = lr_network_top_matrix[order_receptors, order_ligands_receptor]
  rownames(vis_ligand_receptor_network) = order_receptors %>% make.names()
  colnames(vis_ligand_receptor_network) = order_ligands_receptor %>% make.names()

  ## fold-change of ligands in sending cells
  if (!if_any_sender) {
    DE_table_all = Idents(sdata) %>% levels() %>% intersect(sender_celltypes) %>% lapply(get_lfc_celltype, seurat_obj = sdata, condition_colname = "Condition", condition_oi = condition_oi, condition_reference = condition_reference, celltype_col = "Celltypes",expression_pct = pct_expr_threshold) %>% reduce(full_join)
    DE_table_all[is.na(DE_table_all)] = 0
    # Combine ligand activities with DE information
    ligand_activities_de = ligand_activities %>% select(test_ligand, pearson) %>% rename(ligand = test_ligand) %>% left_join(DE_table_all %>% rename(ligand = gene), by = "ligand")
    ligand_activities_de[is.na(ligand_activities_de)] = 0
    # make LFC heatmap
    lfc_matrix = ligand_activities_de  %>% select(-ligand, -pearson) %>% as.matrix() %>% magrittr::set_rownames(ligand_activities_de$ligand)
    rownames(lfc_matrix) = rownames(lfc_matrix) %>% make.names()
    order_ligands = order_ligands[order_ligands %in% rownames(lfc_matrix)]
    vis_ligand_lfc = lfc_matrix[order_ligands, , drop = FALSE] # need drop argument in case only have 1 cell type, or else returns vector instead of matrix
    colnames(vis_ligand_lfc) = vis_ligand_lfc %>% colnames() %>% make.names()
    
    vis_ligand_lfc = vis_ligand_lfc[rownames(vis_ligand_lfc) %in% included_lig, , drop = FALSE]  # match entries shown in prior plots
    
    p_ligand_lfc = vis_ligand_lfc %>% make_threecolor_heatmap_ggplot("Prioritized ligands","Sender cells", low_color = "midnightblue",mid_color = "white", mid = median(vis_ligand_lfc), high_color = "red",legend_position = "bottom", x_axis_position = "top", legend_title = "log fold-change") + 
      theme(axis.text.y = element_text(face = "italic"),
            axis.text.x.top = element_text(angle = 0, hjust = 0.5))
    # rbind(as.vector(subset(sdata, idents = "Plasmacytoid DC")@assays$RNA['JAML', ]), subset(sdata, idents = "Plasmacytoid DC")$Condition) %>% t() %>% as_tibble() %>% group_by(V2) %>% summarize(mean(as.numeric(V1)))  # double-check that negative numbers from the DE calculation means down in tumor
  } else {
    vis_ligand_lfc = c()
    p_ligand_lfc = NA
  }
  
  ## combined plot output
  ligand_pearson_matrix = ligand_activities %>% select(pearson) %>% as.matrix() %>% magrittr::set_rownames(ligand_activities$test_ligand)
  rownames(ligand_pearson_matrix) = rownames(ligand_pearson_matrix) %>% make.names()
  colnames(ligand_pearson_matrix) = colnames(ligand_pearson_matrix) %>% make.names()
  
  vis_ligand_pearson = ligand_pearson_matrix[order_ligands, ] %>% as.matrix(ncol = 1) %>% magrittr::set_colnames("Pearson")
  if (!if_any_sender) {
    vis_ligand_pearson = vis_ligand_pearson[rownames(vis_ligand_pearson) %in% included_lig, , drop = FALSE]
  }
  p_ligand_pearson = vis_ligand_pearson %>% make_heatmap_ggplot("Prioritized ligands","Ligand activity", color = "darkorange",legend_position = "bottom", x_axis_position = "top", legend_title = "Pearson correlation coefficient\ntarget gene prediction ability)") + 
    theme(legend.text = element_text(size = 9),
          axis.text.x.top = element_blank(), axis.ticks.x.top = element_blank())
  
  order_ligands_adapted = order_ligands
  order_ligands_adapted[order_ligands_adapted == "H2.M3"] = "H2-M3" 
  order_ligands_adapted[order_ligands_adapted == "H2.T23"] = "H2-T23" 
  if (!if_any_sender) {
    order_ligands_adapted = order_ligands_adapted[order_ligands_adapted %in% included_lig]
  }
  p_lig_express = DotPlot(sdata %>% subset(Celltypes %in% sender_celltypes), features = order_ligands_adapted, cols = "RdYlBu") + coord_flip() + 
    scale_x_discrete(limits = rev) +
    scale_size_area(limits = c(0, 40), oob = scales::squish) + # otherwise makes look like ligand not expressed when really have 10% (because cutoff was set to 0.1)
    scale_color_gradient(low = "#edeada", high = "#346900") +
    theme(legend.text = element_text(size = 10), 
          legend.title = element_text(size = 12),
          legend.box = "horizontal",
          axis.text.x.top = element_text(angle = 0, hjust = 0.5))
  
  return(list(p_all_senders, p_ligand_target_network, p_ligand_lfc, p_ligand_pearson, p_lig_express,
              vis_ligand_lfc, vis_ligand_pearson, vis_ligand_target))
}


nichenet_output = nichenet_tumor_analysis(sdata, receiver, sender_celltypes,
                                          expressed_genes_receiver, background_expressed_genes,
                                          list_expressed_genes_sender, expressed_genes_sender,
                                          pct_expr_threshold, logfc_threshold, 
                                          if_any_sender)
# figures_combined
p_all_senders = nichenet_output[[1]] 
p_ligand_target_network = nichenet_output[[2]]
p_ligand_lfc = nichenet_output[[3]]
p_ligand_pearson = nichenet_output[[4]]
p_lig_express = nichenet_output[[5]]
vis_ligand_lfc = nichenet_output[[6]]
vis_ligand_pearson = nichenet_output[[7]]
vis_ligand_target = nichenet_output[[8]]

# combine plots into figure
figures_without_legend = cowplot::plot_grid(
  p_lig_express + theme(legend.position = "none", axis.ticks = element_blank(), axis.title.x = element_text(size = 12), axis.text.y = element_text(face = "italic", size = 9), axis.text.x = element_text(size = 9,  angle = 90,hjust = 0)) + ylab("Expression in Sender") + xlab("") + scale_y_discrete(position = "right"),
  p_ligand_lfc + theme(legend.position = "none", axis.ticks = element_blank()) + theme(axis.title.x = element_text()) + ylab(""),
  p_ligand_pearson + theme(legend.position = "none", axis.ticks = element_blank()) + theme(axis.title.x = element_text()),
  p_ligand_target_network + theme(legend.position = "none", axis.ticks = element_blank()) + ylab(""),
  align = "hv",
  nrow = 1,
  rel_widths = c(ncol(vis_ligand_lfc) + 7, ncol(vis_ligand_lfc) + 8, ncol(vis_ligand_pearson)+6, ncol(vis_ligand_target)))
legends = cowplot::plot_grid(
  ggpubr::as_ggplot(ggpubr::get_legend(p_lig_express)),
  ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_lfc)),
  ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_pearson)),
  ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_target_network)),
  nrow = 1,
  align = "h", rel_widths = c(1, 1, 1.5, 1))
combined_plot = cowplot::plot_grid(figures_without_legend, legends, rel_heights = c(10,5), nrow = 2, align = "hv")
combined_plot
#ggsave("kidney_scRNA_tumor_nichenet_Th_receiver_pDC_sender_leukonly_combined.png", combined_plot, width = 15, height = 6)


#### analysis for any sender and pDC receiver ####

pct_expr_threshold = 0.10
logfc_threshold = 0.25 
if_any_sender = TRUE
## receiver
receiver = "Plasmacytoid DC"
expressed_genes_receiver = get_expressed_genes(receiver, sdata, pct = pct_expr_threshold)
background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix_leuk)]
## sender
sender_celltypes = unique(sdata@meta.data$Celltypes)
list_expressed_genes_sender = sender_celltypes %>% unique() %>% lapply(get_expressed_genes, sdata, pct_expr_threshold)
expressed_genes_sender = list_expressed_genes_sender %>% unlist() %>% unique()


nichenet_output = nichenet_tumor_analysis(sdata, receiver, sender_celltypes,
                                          expressed_genes_receiver, background_expressed_genes,
                                          list_expressed_genes_sender, expressed_genes_sender,
                                          pct_expr_threshold, logfc_threshold, 
                                          if_any_sender)
# figures_combined
p_all_senders = nichenet_output[[1]] 
p_ligand_target_network = nichenet_output[[2]]
p_ligand_lfc = nichenet_output[[3]]
p_ligand_pearson = nichenet_output[[4]]
p_lig_express = nichenet_output[[5]]
vis_ligand_lfc = nichenet_output[[6]]
vis_ligand_pearson = nichenet_output[[7]]
vis_ligand_target = nichenet_output[[8]]



figures_combined = cowplot::plot_grid(
  p_ligand_pearson,
  p_ligand_target_network, 
  align = "h",
  nrow = 1,
  rel_widths = c(ncol(vis_ligand_pearson)+6, ncol(vis_ligand_target)))
figures_combined
#ggsave("kidney_scRNA_tumor_nichenet_pDC_receiver_any_sender_leukonly_combined.png", figures_combined, width = 15, height = 6)

