library(tidyverse)

# affinity and protein expression analyses from figure 2

setwd("./organized_code_affinity")

affinity_csv = read_csv("Leuk_interaction_affinities.csv")
affinity_list = affinity_csv %>%
  mutate(Uniprot_unique = ifelse(as.character(Uniprot_L) < as.character(Uniprot_R), paste(Uniprot_L, "+", Uniprot_R), paste(Uniprot_R, "+", Uniprot_L))) %>%
  filter(!is.na(Kd_nM))


#### merge with data key files to look at cell-cell pair's interactions ####

fname =  "proteomics_pct05_filt0_intxn_leukmanual"
dkey_fname = paste0("circos_", fname,"_KEY.csv")

which_cellstate = "both"  # "steadystate"   "activated"
if_combine_cellcategories = TRUE   # merges populations for visualizations that focus on simpler subsets of cells focusing on main cell types only
if_keep_fullname = FALSE   # keeps the "_steadystate" etc after end of chr1/chr2 instead of giving cell name only
if_cpdb_proteomics = FALSE  # selection of the interaction dataset to use when picking the datatable formatting (default is false for leuk final draft network)

if (which_cellstate %in% c("steadystate", "activated")) {
  data_key_in = read_csv(dkey_fname) %>%
    mutate(chr1_state = gsub("^.*_", "", chr1),
           chr1 = gsub("_.*$", "", chr1),
           chr2_state = gsub("^.*_", "", chr2),
           chr2 = gsub("_.*$", "", chr2))
  data_key_in = filter(data_key_in,
                       chr1_state == which_cellstate,
                       chr2_state == which_cellstate)
  fname = gsub("_both_", which_cellstate, fname)
  imgheight = 20  # dimensions for plot images
  imgwidth = 20
} else if (which_cellstate == "both") {
  if (if_keep_fullname){
    data_key_in = read_csv(dkey_fname) %>%
      mutate(chr1_state = gsub("^.*_", "", chr1),
             chr2_state = gsub("^.*_", "", chr2))  
  } else {
    data_key_in = read_csv(dkey_fname) %>%
      mutate(chr1_state = gsub("^.*_", "", chr1),
             chr1 = gsub("_.*$", "", chr1),
             chr2_state = gsub("^.*_", "", chr2),
             chr2 = gsub("_.*$", "", chr2)) 
  }
  
  imgheight = 25
  imgwidth = 25 
} else {
  print("invalid cell state entered")
}

if (if_combine_cellcategories) {
  if (which_cellstate != "both") {
    mapping_combine_cellcategories = unique(data_key_in$chr1)
    names(mapping_combine_cellcategories) = mapping_combine_cellcategories # avoid period symbol in names
    mapping_combine_cellcategories = gsub("Th2_|Th17_|Th1_|T4..*_", "T4_", mapping_combine_cellcategories)
    mapping_combine_cellcategories = gsub("B.memory_|B.naive_", "B_", mapping_combine_cellcategories)
    mapping_combine_cellcategories = gsub("NK..*_", "NK_", mapping_combine_cellcategories)
    mapping_combine_cellcategories = gsub("nTregs_|mTregs_", "Treg_", mapping_combine_cellcategories)
    mapping_combine_cellcategories = gsub("T8..*_", "T8_", mapping_combine_cellcategories)
    mapping_combine_cellcategories = gsub("B.plasma_", "B Plasma_", mapping_combine_cellcategories) 
    mapping_combine_cellcategories = gsub("mDC_|pDC_", "DC_", mapping_combine_cellcategories)
    mapping_combine_cellcategories = gsub("MO..*_", "Monocyte_", mapping_combine_cellcategories)
  } else{
    mapping_combine_cellcategories = unique(data_key_in$chr1)
    names(mapping_combine_cellcategories) = mapping_combine_cellcategories
    mapping_combine_cellcategories = gsub("Th2|Th17|Th1|T4..*", "T4", mapping_combine_cellcategories)
    mapping_combine_cellcategories = gsub("B.memory|B.naive", "B", mapping_combine_cellcategories)
    mapping_combine_cellcategories = gsub("NK..*", "NK", mapping_combine_cellcategories)
    mapping_combine_cellcategories = gsub("nTregs|mTregs", "Treg", mapping_combine_cellcategories)
    mapping_combine_cellcategories = gsub("T8..*", "T8", mapping_combine_cellcategories)
    mapping_combine_cellcategories = gsub("B.plasma", "B Plasma", mapping_combine_cellcategories) 
    mapping_combine_cellcategories = gsub("mDC|pDC", "DC", mapping_combine_cellcategories)
    mapping_combine_cellcategories = gsub("MO..*", "Monocyte", mapping_combine_cellcategories)
  }
  
  data_key_in = data_key_in %>%
    mutate(chr1_simple = plyr::revalue(chr1, mapping_combine_cellcategories),
           chr2_simple = plyr::revalue(chr2, mapping_combine_cellcategories))
  fname = paste0(fname, "_simple")
}

if (if_cpdb_proteomics) {
  data_key_kd = data_key_in %>%
    left_join(affinity_list, by = "Uniprot_unique")
} else {
  data_key_kd = data_key_in
}


if (if_keep_fullname) {
  data_key_kd_summary = data_key_kd %>%
    filter(!is.na(Kd_nM)) %>%
    group_by(chr1, chr2) %>%
    summarize(Kd_mean = mean(Kd_nM),
              Kd_med = median(Kd_nM),
              Kd_medlog10 = mean(log10(Kd_nM)),
              Kd_meanlog10 = median(log10(Kd_nM))) %>%
    ungroup()
} else {
  data_key_kd_summary = data_key_kd %>%
    filter(!is.na(Kd_nM)) %>%
    group_by(chr1, chr2, chr1_state, chr2_state) %>%
    summarize(Kd_mean = mean(Kd_nM),
              Kd_med = median(Kd_nM),
              Kd_medlog10 = mean(log10(Kd_nM)),
              Kd_meanlog10 = median(log10(Kd_nM))) %>%
    ungroup()
}

if (if_combine_cellcategories) {
  if (if_keep_fullname) {
    data_key_kd_summary_simple = data_key_kd %>%
      filter(!is.na(Kd_nM)) %>%
      group_by(chr1_simple, chr2_simple) %>%
      summarize(Kd_mean_simple = mean(Kd_nM),
                Kd_med_simple = median(Kd_nM),
                Kd_medlog10_simple = mean(log10(Kd_nM)),
                Kd_meanlog10_simple = median(log10(Kd_nM))) %>%
      ungroup()
  } else {
    data_key_kd_summary_simple = data_key_kd %>%
      filter(!is.na(Kd_nM)) %>%
      group_by(chr1_simple, chr2_simple, chr1_state, chr2_state) %>%
      summarize(Kd_mean_simple = mean(Kd_nM),
                Kd_med_simple = median(Kd_nM),
                Kd_medlog10_simple = mean(log10(Kd_nM)),
                Kd_meanlog10_simple = median(log10(Kd_nM))) %>%
      ungroup()
  }
}

if (if_cpdb_proteomics) {
  data_key_kd %<>%
    mutate(cell_pair = ifelse(chr1 < chr2, paste0(chr1, " + ", chr2), paste0(chr2, " + ", chr1)),
           state_pair = ifelse(chr1_state < chr2_state, paste0(chr1_state, " + ", chr2_state), paste0(chr2_state, " + ", chr1_state)),
           intxn_pair = ifelse(Gene_L < Gene_R, paste0(Gene_L, " + ", Gene_R), paste0(Gene_R, " + ", Gene_L))) %>%
    mutate(cell_upair = paste0(chr1, " + ", chr2),
           state_upair = paste0(chr1_state, " + ", chr2_state)) %>%
    filter(!is.na(Kd_nM)) %>%
    {if (which_cellstate == "both" & if_keep_fullname) left_join(., data_key_kd_summary, by = c("chr1", "chr2")) else left_join(., data_key_kd_summary, by = c("chr1", "chr2", "chr1_state", "chr2_state"))}
} else {
  data_key_kd %<>%
    mutate(cell_pair = ifelse(chr1 < chr2, paste0(chr1, " + ", chr2), paste0(chr2, " + ", chr1)),
           state_pair = ifelse(chr1_state < chr2_state, paste0(chr1_state, " + ", chr2_state), paste0(chr2_state, " + ", chr1_state)),
           intxn_pair = ifelse(gene_name_1 < gene_name_2, paste0(gene_name_1, " + ", gene_name_2), paste0(gene_name_2, " + ", gene_name_1))) %>%
    mutate(cell_upair = paste0(chr1, " + ", chr2),
           state_upair = paste0(chr1_state, " + ", chr2_state)) %>%
    filter(!is.na(Kd_nM)) %>%
    filter(!grepl(",B2M$", gene_name_1) & !grepl(",B2M$", gene_name_2)) %>%  # secreted beta 2 microglobulin not surface
    {if (if_combine_cellcategories & if_keep_fullname) left_join(., data_key_kd_summary_simple, by = c("chr1_simple", "chr2_simple")) else . } %>%
    {if (if_combine_cellcategories & !if_keep_fullname) left_join(., data_key_kd_summary_simple, by = c("chr1_simple", "chr2_simple", "chr1_state", "chr2_state")) else . } %>%
    {if (which_cellstate == "both" & if_keep_fullname) left_join(., data_key_kd_summary, by = c("chr1", "chr2")) else left_join(., data_key_kd_summary, by = c("chr1", "chr2", "chr1_state", "chr2_state"))}
}


#### affinity distribution plots ####

shortening_name_dict = unique(data_key_kd$chr1)
names(shortening_name_dict) = shortening_name_dict  # shorten the names so can fit labels onto the figures
shortening_name_dict = gsub("MO.nonclassical", "MO.nonclass", shortening_name_dict)
shortening_name_dict = gsub("MO.intermediate", "MO.intermed", shortening_name_dict)
shortening_name_dict = gsub("MO.classical",    "MO.classic", shortening_name_dict)
shortening_excl_types = c("Basophil", "Neutrophil", "Eosinophil") # for shortened plot to fit in small figures, only including main cell type groups

if_combstate = FALSE  # plots both active and steady state
if_activeonly = FALSE  # if false and 'combstate' also false, then plots only steady state
if_simplycelltypes = TRUE
which_plot_type = "histogram"  # "density"  "tile"
if_linear_scale = FALSE  # linear or log scale

data_key_kd_selected = data_key_kd %>%
  mutate(chr1 = plyr::revalue(chr1, shortening_name_dict),
         chr2 = plyr::revalue(chr2, shortening_name_dict)) %>%
  mutate(chr1_cond = paste0(chr1, ifelse(chr1_state == "activated", " (active)", " (resting)")),
         chr2_cond = paste0(chr2, ifelse(chr2_state == "activated", " (active)", " (resting)"))) %>%
  mutate(cell_pair_simple = ifelse(chr1_simple < chr2_simple, paste0(chr1_simple," + ",chr2_simple), paste0(chr2_simple," + ",chr1_simple))) %>%
  {if (if_combstate | if_activeonly) . else filter(., state_upair == "steady.state + steady.state")} %>%
  {if (if_activeonly) filter(., state_upair == "activated + activated") else .} %>%
  {if (if_simplycelltypes) ungroup(mutate(group_by(filter(., !chr1_simple %in% shortening_excl_types, !chr2_simple %in% shortening_excl_types), cell_pair_simple), Kd_med_types = median(Kd_med), Kd_mean_types = mean(Kd_meanlog10))) else . }
  
ggplot(data_key_kd_selected) +
    {if(which_plot_type == "histogram" & !if_simplycelltypes & !if_linear_scale) geom_histogram(aes(x = log10(Kd_nM), fill = 10^(Kd_medlog10)), bins = 20)} +
    {if(which_plot_type == "histogram" & if_simplycelltypes& !if_linear_scale) geom_histogram(aes(x = log10(Kd_nM), y = ..density.., fill = (Kd_med_types)), bins = 20)} +
    {if(which_plot_type == "density" & !if_simplycelltypes & !if_linear_scale) geom_density(aes(x = log10(Kd_nM), fill = 10^(Kd_medlog10)), adjust = 1)} +
    {if(which_plot_type == "density" & if_simplycelltypes & !if_linear_scale) geom_density(aes(x = log10(Kd_nM), fill = (Kd_med_types)), adjust = 1)} +
    
    {if(which_plot_type %in% c("histogram", "density", "tile")) theme_minimal() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())} + # for histogram plots 
    {if(which_plot_type %in% c("histogram", "density") & !if_simplycelltypes)facet_grid(chr1 ~ chr2)} +
    {if(which_plot_type %in% c("histogram", "density") & if_simplycelltypes)facet_grid(chr1_simple ~ chr2_simple)} +
    {if(which_plot_type %in% c("histogram", "density")) xlab("log10 Kd (nM)")} +
  
    {if(which_plot_type == "tile" & if_simplycelltypes)geom_tile(aes(x = 1, y = 1, fill = (Kd_med_types)))} +  #fill = (Kd_med_types)  
    {if(which_plot_type == "tile" & !if_simplycelltypes)geom_tile(aes(x = 1, y = 1, fill = 10^(Kd_medlog10)))} +
    #{if(which_plot_type == "tile") theme_dark()} +
    {if(which_plot_type == "tile" & if_simplycelltypes) facet_grid(chr1_simple ~ chr2_simple)} +
    
    {if(if_combstate)facet_grid(chr1_cond ~ chr2_cond)} +
    {if(!if_combstate&!if_simplycelltypes)facet_grid(chr1 ~ chr2)} +
    #scale_fill_viridis_c(name = "Average Kd (nM)", limits = c(1000, 3800), oob = scales::squish) +  # limits = c(1500, 4000)
    scale_fill_gradient(low = "#2a1d59", high = "#d64933", name = "Average Kd (nM)", limits = c(1500, 3500), oob = scales::squish) +  # for when do raw medians 
    {if(!if_combstate&!if_simplycelltypes)theme(panel.spacing = unit(0.10, "lines"))} + #unit(1.005, "lines"))} +  # 0.80 spacing looks good for active only
    {if(!if_combstate&if_simplycelltypes) theme(panel.spacing = unit(0, "lines"))} +
    {if(if_combstate) theme(panel.spacing = unit(0.1, "lines"))} +
    {if(!if_simplycelltypes&!if_combstate)theme(strip.text = element_text(size = 4))} +  # shrinks facet label text
    {if(if_combstate)theme(strip.text = element_text(size = 2.5))}  # shrinks facet label text
if(which_plot_type == "histogram") {
  ggsave(paste0("affinity_", fname, ifelse(if_simplycelltypes, "_simple_","_"), ifelse(if_combstate, "comb", ifelse(if_activeonly, "active","resting")), "_histograms.png"), width = imgwidth, height = imgheight)
} else if (which_plot_type == "density") {
  ggsave(paste0("affinity_", fname, ifelse(if_simplycelltypes, "_simple_","_"), ifelse(if_combstate, "comb", ifelse(if_activeonly, "active","resting")), "_density.png"), width = imgwidth, height = imgheight)
} else if (which_plot_type == "tile") {
  ggsave(paste0("affinity_", fname, ifelse(if_simplycelltypes, "_simple_","_"), ifelse(if_combstate, "comb", ifelse(if_activeonly, "active","resting")), "_linmedtypes.png"), width = imgwidth, height = imgheight)
} else {
  print("Plot type selection not valid")
}

kruskal.test(Kd_nM ~ cell_pair_simple, data = data_key_kd_selected)
pairwise.wilcox.test(data_key_kd_selected$Kd_nM, data_key_kd_selected$cell_pair_simple,
                     p.adjust.method = "BH")

#### same distributions but comparing APC to T cell interaction pairs that already know are important ####

data_key_kd = data_key_kd %>%
  mutate(chr1_andstate = gsub("^.* : ","",chr1_unique),
         chr2_andstate = gsub("^.* : ","",chr2_unique)) %>%
  mutate(cell_pair_state = paste(chr1_andstate, "+", chr2_andstate))
category_T_DC = c(grep(unique(x = data_key_kd$cell_pair), pattern = "^(pDC|mDC)(.*) \\+ (T|mTreg|nTreg)", value = T),  # gathers all the names of dendritic cells (grep just saves time/more versatile then manually writing out all the names)
                  grep(unique(x = data_key_kd$cell_pair), pattern = "^(T|mTreg|nTreg)(.*) \\+ (pDC|mDC)", value = T))
category_T_B = c(grep(unique(x = data_key_kd$cell_pair), pattern = "^B\\.(.*) \\+ (T|mTreg|nTreg)", value = T),
                 grep(unique(x = data_key_kd$cell_pair), pattern = "^(T|mTreg|nTreg)(.*) \\+ B\\.(.*)", value = T))

data_key_category_medKd = data_key_kd %>%
  mutate(pair_category = case_when(
    cell_pair %in% category_T_DC ~ "T_DC",
    cell_pair %in% category_T_B ~ "T_B",
    TRUE ~ "Others"
  )) %>% 
  distinct(cell_pair, state_pair, Kd_mean, Kd_med, Kd_meanlog10, Kd_medlog10, pair_category)  # want one point per cell type pair, not duplicated ones per protein

PrintScientific = function(val) {
  exponent = floor(log10(abs(val)))
  base = val * 10^(-1.0*exponent)
  sprintf("%.1fE%+01d", base, exponent)}

data_key_category_medKd %>% 
  filter(pair_category != "Others") %>%
  ggplot(aes(x = pair_category, y = 10^Kd_meanlog10)) +
  geom_violin(aes(fill = pair_category), draw_quantiles = c(0.25, 0.5, 0.75), adjust = 1.2) +  # adjust to make violin curves more wide
  annotate("text", x = 1, y = 4200, label = sprintf("t = %.2f\np = %s", 
                   t.test(10^(data_key_category_medKd$Kd_meanlog10[data_key_category_medKd$pair_category == "T_DC"]), 10^(data_key_category_medKd$Kd_meanlog10[data_key_category_medKd$pair_category == "T_B"]), var.equal = FALSE)$statistic,
                   PrintScientific(t.test(10^(data_key_category_medKd$Kd_meanlog10[data_key_category_medKd$pair_category == "T_DC"]), 10^(data_key_category_medKd$Kd_meanlog10[data_key_category_medKd$pair_category == "T_B"]), var.equal = FALSE)$p.value))) + 
  #geom_dotplot(binaxis='y', stackdir='center', stackratio = 0.9, dotsize = 0.2, alpha = 0.5)
  geom_jitter(alpha = 0.8, height = 100, width = 0.1) +
  theme_bw() + ylim(c(500, 4500)) + guides(fill = FALSE) +
  xlab("Activated APC and resting T lymphocyte pair") + ylab("Median Kd of cell pair surface interactions (nM)")
#ggsave(paste0("affinity_", fname, "_resting_violin_categories_activeAPCT_kdmeanlog10.png"), width = 4, height = 3.5)


#### merge with quantitative data ####

rieckmann_means = read_csv("Rieckmann_proteomics_summary_mean_leuk.csv") %>%
  pivot_longer(`B.memory_activated`:`Th2_steady.state`, names_to = "Celltype", values_to = "Expression")

affinity_list %>%
  fuzzyjoin::fuzzy_inner_join(rieckmann_means, by = c("Uniprot_unique" = "Uniprot"), match_fun = str_detect) %>%
  filter(Expression > 1) %>%
  ggplot(aes(x = log10(Expression), y = log10(Kd_nM))) +
    geom_point(aes(color = Uniprot), alpha = 0.7) +
    scale_color_grey() +
    facet_wrap( ~ Celltype, nrow = 5) +
    stat_smooth(method = "lm", formula = y ~ x) + 
    ggpubr::stat_cor(method = "pearson", size = 2) + 
    theme_minimal() +
    theme(legend.position = "none") 
#ggsave("affinity_expression_celltypes_correlate_log10.png", width = 20, height = 9)

affinity_list %>%
  inner_join(rieckmann_means, by = c("Uniprot_L" = 'Uniprot')) %>%
  inner_join(rieckmann_means, by = c("Uniprot_R" = 'Uniprot'), suffix = c("_L", "_R")) %>%
  filter(Expression_L > 1, Expression_R > 1) %>%
  mutate(pair_label = paste(Gene_L, "-", Gene_R),
         Expression_sum = Expression_L + Expression_R,
         Expression_min = pmin(Expression_L, Expression_R),
         Expression_max = pmax(Expression_L, Expression_R)) %>%
  ggplot(aes(x = log10(Expression_sum), y = log10(Kd_nM))) +
    geom_hex(bins = 15) +
    scale_color_grey() +
    stat_smooth(method = "lm", formula = y ~ x) + 
    ggpubr::stat_cor(method = "pearson", size = 2) + 
    theme_minimal()
#ggsave("affinity_sumexpression_cellpairs_correlate_hex_log.png", width = 7, height = 6) # remember this has many points per intxn, one point per intxn per cell type pair


#### affinities of differentially expressed proteins ####

library(DESeq2)

rieckmann_raw_csv = read_csv("rieckmann_proteomics_summed.csv")

rieckmann_raw = rieckmann_raw_csv %>%
  rename_with(~gsub("steady-state", "steady.state", .x)) # remove the forbidden dash symbol

rieckmann_raw_matrix = as.matrix(select(rieckmann_raw, -Majority.protein.IDs))
uniprot_row_vec = rieckmann_raw$Majority.protein.IDs
sample_col_vec = colnames(rieckmann_raw_matrix)

rieckmann_metadata = rieckmann_raw %>%
  select(-Majority.protein.IDs) %>%
  pivot_longer(cols = `CopyNumber_B.memory_01_activated`:`CopyNumber_NK.dim_04_steady.state`, names_to = "sample_condition", values_to = "Expression") %>%
  select(-Expression) %>%
  distinct() %>%
  mutate(sample_condition = gsub("^CopyNumber_", "", sample_condition)) %>%
  separate(sample_condition, into = c("Celltype", "Replicate", "Cellstate"), sep = "_", remove = FALSE) %>%
  mutate(Replicate = as.numeric(Replicate)) %>%
  as.data.frame()
rownames(rieckmann_metadata) = rieckmann_metadata$sample_condition

celltype_list = rieckmann_metadata %>%
  group_by(Celltype) %>%
  summarize(if_bothstate = n() > 5, .groups = "drop") %>% # 4 replicates per activation condition, so if have 4 or fewer points that means cell type was only measured in the steady-state condition and lacks paired data of those cells in activated condition
  ungroup() %>%
  filter(if_bothstate) %>%
  select(Celltype) %>%
  as.vector()

deseq_combined = tibble(baseMean = as.numeric(), log2FoldChange = as.numeric(), lfcSE = as.numeric(),
                        stat = as.numeric(), pvalue = as.numeric(), padj = as.numeric(), 
                        Uniprot_list = as.character(), Celltype = as.character())
for (celltype_pick in celltype_list$Celltype) {
  print(celltype_pick)
  picked_matrix = rieckmann_raw_matrix[ , grepl(celltype_pick, colnames(rieckmann_raw_matrix))]
  picked_metadata = rieckmann_metadata[grepl(celltype_pick, rownames(rieckmann_metadata)), ]
  deseq_mat = DESeqDataSetFromMatrix(countData = picked_matrix,
                                     colData = picked_metadata,
                                     design = ~ Cellstate)
  
  deseq_out = DESeq(deseq_mat) # resultsNames(deseq_out)
  deseq_result = results(deseq_out, name = "Cellstate_steady.state_vs_activated")
  
  deseq_tbl = as_tibble(deseq_result) %>%
    mutate(Uniprot_list = uniprot_row_vec) %>%
    mutate(Celltype = celltype_pick)
  
  deseq_combined = bind_rows(deseq_combined, deseq_tbl)
}
#saveRDS(deseq_combined, "DESeq_Rieckmann_proteomics_activation.RDS")

deseq_combined = readRDS("DESeq_Rieckmann_proteomics_activation.RDS")

deseq_receptors = affinity_list %>%
  fuzzyjoin::fuzzy_left_join(deseq_combined, c("Uniprot_L" = "Uniprot_list"), match_fun = str_detect) %>%
  fuzzyjoin::fuzzy_left_join(deseq_combined, c("Uniprot_R" = "Uniprot_list"), match_fun = str_detect) %>%
  rename_with(~gsub(".x$", "_L", .x)) %>% 
  rename_with(~gsub(".y$", "_R", .x))

p_thresh = 0.10
logfc_thresh = 1.1  #fold-change measured on a log2 scale, so to set threshold of more than 2-fold corresponds to number more than 1
if_pthresh = FALSE  # alternatively can select a threshold based on the p value threshold or based on the fold change threshold
deseq_sig_receptors = deseq_receptors %>%
  {if (if_pthresh) filter(., padj_L < p_thresh) else filter(., abs(log2FoldChange_L) > logfc_thresh)} %>%
  select(-matches("_R$")) %>%
  rename_with(~gsub("_L", "", .x))

celltype_has_upanddown = deseq_sig_receptors %>%  # handle situations where all DE is only upreg or only downreg and can't compare
  group_by(Celltype) %>%
  summarize(analyzable = any(stat > 0) & any(stat < 0)) %>%
  ungroup() 

deseq_sig_receptors = deseq_receptors %>%
  {if (if_pthresh) filter(., padj_R < p_thresh) else filter(., abs(log2FoldChange_R) > logfc_thresh)} %>%
  select(-matches("_L$")) %>%
  rename_with(~gsub("_R", "", .x)) %>%
  bind_rows(deseq_sig_receptors) %>%
  filter(!is.na(Celltype)) %>%
  filter(Celltype %in% celltype_has_upanddown[celltype_has_upanddown$analyzable, 'Celltype'][['Celltype']]) %>%
  distinct(Uniprot_unique, Celltype, Uniprot_list, .keep_all = TRUE) %>%
  mutate(Moreactive = log2FoldChange < 0) %>% # if positive, means more expression in resting condition ; if negative, more in active condition
  group_by(Celltype) %>%
  mutate(T_compare_statistic = t.test(Kd_nM[Moreactive], Kd_nM[!Moreactive], var.equal = FALSE)$statistic,
         p_compare_raw = t.test(Kd_nM[Moreactive], Kd_nM[!Moreactive], var.equal = FALSE)$p.value) %>%
  ungroup() %>%
  group_by(Celltype, Moreactive) %>%
  mutate(Kd_nM_mean = mean(Kd_nM),
         Kd_nM_med = median(Kd_nM)) %>%
  ungroup() %>%
  mutate(Moreactive = factor(Moreactive, levels = c(F,T), labels = c("Resting", "Active")))

ggplot(deseq_sig_receptors,
       aes(x = Moreactive, y = log10(Kd_nM))) +
  geom_boxplot(aes(fill = log10(Kd_nM_med)), width = 0.7) + 
  geom_dotplot(binaxis='y', stackdir='center', stackratio = 1.1, dotsize = 0.5, alpha = 0.5, fill = "grey") + 
  geom_text(data = mutate(distinct(deseq_sig_receptors, Celltype, p_compare_raw, T_compare_statistic), p_compare_adj = p.adjust(p_compare_raw)), 
            #mapping = aes(x = -Inf, y = +Inf, label = sprintf("t = %.3f\np = %.3f", T_compare_statistic, p_compare_raw)), hjust = -0.1, vjust = +1.2, size = 2.5) +
            mapping = aes(x = -Inf, y = +Inf, label = sprintf("p = %.3f", p_compare_raw)), hjust = -0.5, vjust = +1.9, size = 2.5) +
  facet_wrap( ~ Celltype, nrow = 1) +
  scale_fill_gradient(low = "#84ad68", high = "#1a65c7") + 
  theme_bw() + xlab("Differential expression condition")
ggsave(paste0("affinity_proteomics_differential_",ifelse(!if_pthresh, paste0("logFC", logfc_thresh), paste0("pval", gsub("\\.","",as.character(sprintf(p_thresh, fmt = "%0.2f"))))),"_celltypes.png"),
       width = 20, height = 6)

deseq_sig_receptors %>%
  mutate(T_compare_statistic = t.test(Kd_nM[Moreactive == "Resting"], Kd_nM[Moreactive == "Active"], var.equal = FALSE)$statistic,
         p_compare_raw = t.test(Kd_nM[Moreactive == "Resting"], Kd_nM[Moreactive == "Active"], var.equal = FALSE)$p.value) %>%
  ggplot(aes(x = Moreactive, y = log10(Kd_nM))) +
  geom_boxplot(notch = TRUE, notchwidth = 0.7) +
  geom_dotplot(binaxis='y', stackdir='center', stackratio = 1.1, dotsize = 0.5, alpha = 0.5, fill = "grey") + 
  geom_text(data = distinct(mutate(deseq_sig_receptors,
                                   T_compare_statistic = t.test(Kd_nM[Moreactive == "Resting"], Kd_nM[Moreactive == "Active"], var.equal = FALSE)$statistic,
                                   p_compare_raw = t.test(Kd_nM[Moreactive == "Resting"], Kd_nM[Moreactive == "Active"], var.equal = FALSE)$p.value), Moreactive, p_compare_raw, T_compare_statistic), 
            #mapping = aes(x = -Inf, y = +Inf, label = sprintf("t = %.3f\np = %.3f", T_compare_statistic, p_compare_raw)), hjust = -0.1, vjust = +1.2, size = 2.5) +
            mapping = aes(x = -Inf, y = +Inf, label = sprintf("t = %.3f\np < 0.001",T_compare_statistic)), hjust = -0.5, vjust = +1.9, size = 2.5) +
  theme_bw() + xlab("Differential expression condition")
ggsave(paste0("affinity_proteomics_differential_",ifelse(!if_pthresh, paste0("logFC", logfc_thresh), paste0("pval", gsub("\\.","",as.character(sprintf(p_thresh, fmt = "%0.2f"))))),"_aggregatecelltypes.png"),
       width = 4, height = 4)

#### same differential expression upon activation analysis as above, except combining cell types ####

cellnames_simple_map = tibble(namefull = c(data_key_kd$chr1, data_key_kd$chr2),
                              namesimple = c(data_key_kd$chr1_simple, data_key_kd$chr2_simple)) %>%
  distinct(namefull, namesimple)

rieckmann_metadata = rieckmann_metadata %>%  # add the shortened cell name labels
  mutate(Celltype_simple = plyr::mapvalues(Celltype, from = cellnames_simple_map$namefull, to = cellnames_simple_map$namesimple)) %>%
  group_by(Celltype_simple, Cellstate) %>%
  mutate(Replicate_simple = row_number()) %>%
  as.data.frame()
rownames(rieckmann_metadata) = rieckmann_metadata$sample_condition

celltype_list = rieckmann_metadata %>%
  filter(!Celltype_simple %in% c("Basophil", "B Plasma", "Eosinophil", "Neutrophil", "Treg")) %>%  # removing categories which don't have paired active/ resting categories
  select(Celltype_simple) %>%
  mutate(Celltype_simple = ifelse(Celltype_simple == "Monocyte", "MO", Celltype_simple)) %>%
  distinct() %>%
  as.vector()

deseq_combined = tibble(baseMean = as.numeric(), log2FoldChange = as.numeric(), lfcSE = as.numeric(),
                        stat = as.numeric(), pvalue = as.numeric(), padj = as.numeric(), 
                        Uniprot_list = as.character(), Celltype = as.character())
for (celltype_pick in celltype_list$Celltype_simple) {
  print(celltype_pick)
  picked_matrix = rieckmann_raw_matrix[ , grepl(paste0("_",celltype_pick,"\\."), colnames(rieckmann_raw_matrix))]
  picked_metadata = rieckmann_metadata[grepl(paste0(celltype_pick,"\\."), rownames(rieckmann_metadata)), ]
  if (celltype_pick == "B") {
    picked_matrix = picked_matrix[ , !grepl("B.plasma", colnames(picked_matrix))]  # regex will return B plasma even though they aren't actually a paired celltype with active and resting
    picked_metadata = picked_metadata[!grepl("B.plasma", rownames(picked_metadata)), ]
  }
  if (celltype_pick %in% c("DC")) {
    picked_matrix = rieckmann_raw_matrix[ , grepl(celltype_pick, colnames(rieckmann_raw_matrix))]
    picked_metadata = rieckmann_metadata[grepl(celltype_pick, rownames(rieckmann_metadata)), ]    
  }
  
  deseq_mat = DESeqDataSetFromMatrix(countData = picked_matrix,
                                     colData = picked_metadata,
                                     design = ~ Cellstate)
  
  deseq_out = DESeq(deseq_mat) # resultsNames(deseq_out)
  deseq_result = results(deseq_out, name = "Cellstate_steady.state_vs_activated")
  
  deseq_tbl = as_tibble(deseq_result) %>%
    mutate(Uniprot_list = uniprot_row_vec) %>%
    mutate(Celltype = celltype_pick)
  
  deseq_combined = bind_rows(deseq_combined, deseq_tbl)
}
#saveRDS(deseq_combined, "DESeq_Rieckmann_proteomics_activation_simplified.RDS")

deseq_combined = readRDS("DESeq_Rieckmann_proteomics_activation_simplified.RDS")

deseq_receptors = affinity_list %>%
  fuzzyjoin::fuzzy_left_join(deseq_combined, c("Uniprot_L" = "Uniprot_list"), match_fun = str_detect) %>%
  fuzzyjoin::fuzzy_left_join(deseq_combined, c("Uniprot_R" = "Uniprot_list"), match_fun = str_detect) %>%
  rename_with(~gsub(".x$", "_L", .x)) %>%  # manually add proper suffix
  rename_with(~gsub(".y$", "_R", .x))

p_thresh = 0.10
logfc_thresh = 1.1
if_pthresh = FALSE
deseq_sig_receptors = deseq_receptors %>%
  {if (if_pthresh) filter(., padj_L < p_thresh) else filter(., abs(log2FoldChange_L) > logfc_thresh)} %>%
  select(-matches("_R$")) %>%
  rename_with(~gsub("_L", "", .x))

celltype_has_upanddown = deseq_sig_receptors %>%  # handle situations where all DE is only upreg or only downreg and can't compare
  group_by(Celltype) %>%
  summarize(analyzable = any(stat > 0) & any(stat < 0)) %>%
  ungroup() 

deseq_sig_receptors = deseq_receptors %>%
  {if (if_pthresh) filter(., padj_R < p_thresh) else filter(., abs(log2FoldChange_R) > logfc_thresh)} %>%
  select(-matches("_L$")) %>%
  rename_with(~gsub("_R", "", .x)) %>%
  bind_rows(deseq_sig_receptors) %>%
  filter(!is.na(Celltype)) %>% # remove intxns with no expn match from when did left join
  filter(Celltype %in% celltype_has_upanddown[celltype_has_upanddown$analyzable, 'Celltype'][['Celltype']]) %>%
  distinct(Uniprot_unique, Celltype, Uniprot_list, .keep_all = TRUE) %>%
  mutate(Moreactive = log2FoldChange < 0) %>% # if positive, means more expression in resting condition ; if negative, more in active condition
  group_by(Celltype) %>%
  mutate(T_compare_statistic = t.test(Kd_nM[Moreactive], Kd_nM[!Moreactive], var.equal = FALSE)$statistic,
         p_compare_raw = t.test(Kd_nM[Moreactive], Kd_nM[!Moreactive], var.equal = FALSE)$p.value) %>%
  ungroup() %>%
  group_by(Celltype, Moreactive) %>%
  mutate(Kd_nM_mean = mean(Kd_nM),
         Kd_nM_med = median(Kd_nM)) %>%
  ungroup() %>%
  mutate(Moreactive = factor(Moreactive, levels = c(F,T), labels = c("Resting", "Active")))

ggplot(deseq_sig_receptors,
       aes(x = Moreactive, y = log10(Kd_nM))) +
  geom_boxplot(aes(fill = log10(Kd_nM_med)), width = 0.7) + 
  geom_dotplot(binaxis='y', stackdir='center', stackratio = 1.1, dotsize = 0.5, alpha = 0.5, fill = "grey") + 
  geom_text(data = distinct(deseq_sig_receptors, Celltype, p_compare_raw, T_compare_statistic), 
            #mapping = aes(x = -Inf, y = +Inf, label = sprintf("t = %.3f\np = %.3f", T_compare_statistic, p_compare_raw)), hjust = -0.1, vjust = +1.2, size = 2.5) +
            mapping = aes(x = -Inf, y = +Inf, label = sprintf("p = %.3f", p_compare_raw)), hjust = -0.5, vjust = +1.9, size = 2.5) +
  facet_wrap( ~ Celltype, nrow = 1) +
  scale_fill_gradient(low = "#84ad68", high = "#1a65c7", name = "log10(Kd) (nM)") + 
  theme_bw() + xlab("Differential expression condition")
ggsave(paste0("affinity_proteomics_differential_",ifelse(!if_pthresh, paste0("logFC", logfc_thresh), paste0("pval", gsub("\\.","",as.character(sprintf(p_thresh, fmt = "%0.2f"))))),"_celltypes_simple.png"),
       width = 11, height = 3.5)

#### loading RNAseq data from pritchard paper instead of Rieckmann proteomics (analysis identical to above) ####

library(biomaRt)
mart <- useMart('ENSEMBL_MART_ENSEMBL')
mart <- useDataset('hsapiens_gene_ensembl', mart) #listDatasets(mart)
annotLookup <- getBM(
  mart = mart,
  attributes = c(
    'ensembl_gene_id',
    'hgnc_symbol',
    'external_gene_name',
    'uniprot_gn_symbol',
    'uniprot_gn_id'),
  uniqueRows=TRUE)

pritchard_DE_list = read_tsv("Supplementary_data_4_RNA_stimulation_DE_genes.txt") %>%  # load so then can slot in to the exact same analysis as above
  dplyr::rename(ENSEMBL = peak_id) %>%
  dplyr::rename(praw = P.Value, padj = adj.P.Val, tval = t, log2FoldChange = logFC) %>%
  mutate(Celltype = gsub("-.*$", "", contrast),
         Celltype = gsub("_S", "", Celltype)) %>%
  mutate(Symbol = ensembldb::mapIds(EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75, keys = `ENSEMBL`, keytype = "GENEID", column = "SYMBOL"))
