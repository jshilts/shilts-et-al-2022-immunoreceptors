library(tidyverse)

setwd("./organized_pharm_analysis")


#### cell activation results analysis ####

combined_cellpct_short = read_csv("combined_cellpct_short.csv")

normtype = "med3"
pvaltype = "all3"
pthresh = 0.10
normname = paste0("Cell_proportion_norm1", normtype)
pvalname = paste0("p_adj_", pvaltype)
combined_cellpct_short %>%
  filter(Timepoint == "24hr") %>%
  filter(!Cell_population %in% c("M14", "M16", "Negs", "T0")) %>%  # activation only measured for B, NK, and T4/T8 cells
  filter(grepl("activat", Cell_population)) %>% # get activation score values
  filter(!Protein %in% c("rCD4", "Buffer", "Carryover")) %>%
  mutate(condition_simpID = paste0(Protein, " (",CellCondition,")")) %>%
  filter(!is.na(!!sym(pvalname))) %>%
  mutate(Cell_population_cond = factor(paste(Cell_population, CellCondition), 
                                       levels = rev(c("T_activation Resting", "T_activation Activated", "NK_activation Resting", "NK_activation Activated", "B_activation Resting", "B_activation Activated")),
                                       labels = rev(c("T cell activation\n(Resting background)", "T cell activation\n(Stimulating background)", "NK cell activation\n(Resting background)", "NK cell activation\n(Stimulating background)", "B cell activation\n(Resting background)", "B cell activation\n(Stimulating background)")))) %>%
  ggplot() +
  geom_point(aes(x = Protein, y = Cell_population_cond, fill = !!sym(normname), 
                 size = -log10(!!sym(pvalname)), color = !!sym(pvalname) < pthresh), pch = 21) +
  scale_size_continuous(limits = c(0,3)) +
  scale_fill_gradient2(midpoint = 0, low = "#363b81", high = "#a80426", mid = "white", limits = c(-0.5,0.5), oob = scales::squish) +
  scale_color_manual(values = c("white", "black"),name = "Significant") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 8, angle = 90, vjust = 0.5),
        axis.text.y = element_text(hjust = 0.5)) +
  xlab("") + ylab("")
#ggsave(paste0("ProtPharm_Trial2_cellpct_norm1", normtype, "_p_adj", pvaltype, "_circles_p", gsub("\\.","",as.character(sprintf(pthresh, fmt = "%0.2f"))),"_comb_shortactivated_longsimp.png"), width = 15, height = 4)


#### cell-to-cell interaction results analysis ####

raw_cellcell_csv = read_csv("ProtPharm_second_CellCellSpatial.csv")

raw_cellcell = raw_cellcell_csv %>%
  mutate(Timepoint = paste0(Timepoint, "hr")) %>% 
  pivot_longer(cols = `B-->B`:`T8-->T8`, names_to = "Interaction_pair", values_to = "Interaction_score") %>%
  mutate(Cell1 = gsub(x = Interaction_pair, pattern = "-->.*$", replacement = ""),
         Cell2 = gsub(x = Interaction_pair, pattern = "^.*-->", replacement = ""))

normavg_cellcell_combdose = raw_cellcell %>%
  mutate(normgroup = 100) %>%
  group_by(normgroup, Protein, CellCondition, Timepoint, Interaction_pair) %>%
  summarize_at(vars(Interaction_score), .funs = median) %>% 
  ungroup() %>%
  group_by(normgroup, CellCondition, Timepoint, Interaction_pair) %>%
  mutate(control_median3 = ifelse(normgroup > 10, median(c(Interaction_score[Protein == "Buffer"], Interaction_score[Protein == "Carryover"], Interaction_score[Protein == "rCD4"])), NA),
         control_mean3 = ifelse(normgroup > 10, mean(c(Interaction_score[Protein == "Buffer"], Interaction_score[Protein == "Carryover"], Interaction_score[Protein == "rCD4"])), NA)) %>%
  mutate(Interaction_score_norm1med3 = ifelse(normgroup == 10, ifelse(grepl("\\+", Protein), (Interaction_score - control_ABprotein) / pmax(Interaction_score, control_ABprotein, na.rm = T), (Interaction_score - control_AB) / pmax(Interaction_score, control_AB, na.rm = T)), (Interaction_score - control_median3) / pmax(Interaction_score, control_median3, na.rm = T)),
         Interaction_score_norm1mean3 = ifelse(normgroup == 10, ifelse(grepl("\\+", Protein), (Interaction_score - control_ABprotein) / pmax(Interaction_score, control_ABprotein, na.rm = T), (Interaction_score - control_AB) / pmax(Interaction_score, control_AB, na.rm = T)), (Interaction_score - control_mean3) / pmax(Interaction_score, control_mean3, na.rm = T))) %>%
  ungroup() %>%
  mutate(conditionID = paste(CellCondition, Timepoint)) %>%
  mutate(conditionID = factor(conditionID, levels = c('Resting 4hr', 'Resting 24hr', 'Activated 4hr', 'Activated 24hr'))) %>%
  mutate(condition_proteinID = paste(conditionID, Protein)) %>%
  rename(Interaction_score_raw = Interaction_score)

#### statistics on cellcell data ####
ctrlavg_raw_cellcell = raw_cellcell %>%
  filter(Protein %in% c("rCD4", "Buffer", "Carryover")) %>%  # the three control conditions (recomb tag, buffer only, and purification column carryover)
  group_by(Protein, CellCondition, Timepoint, Dose, Interaction_pair) %>%
  mutate(Replicate = row_number()) %>%
  ungroup() %>%
  group_by(CellCondition, Timepoint, Dose, Interaction_pair, Replicate) %>%
  summarize(Interaction_score_ctrl_mean = mean(Interaction_score), Interaction_score_ctrl_med = median(Interaction_score), .groups = "drop_last") %>%
  ungroup() %>%
  pivot_longer(cols = contains("ion_score"), names_to = "Protein", values_to = "Interaction_score")

ctrlavg_raw_cellcell = ungroup(mutate(group_by(filter(raw_cellcell, Protein %in% c("rCD4", "Buffer", "Carryover")), Protein, CellCondition, Timepoint, Dose, Interaction_pair), Replicate = row_number())) %>%
  bind_rows(ctrlavg_raw_cellcell) %>%
  rename(ctrl_condition = Protein, ctrl_value = Interaction_score) %>%
  select(-Cell1, -Cell2, -WellID) %>%
  pivot_wider(names_from = ctrl_condition, values_from = ctrl_value)

pvals_cellcell_combdose = raw_cellcell %>%
  group_by(Protein, CellCondition, Timepoint, Dose, Interaction_pair) %>%
  mutate(Replicate = row_number()) %>%
  ungroup() %>%
  filter(Dose != 0) %>%
  filter(!Protein %in% c("rCD4", "Buffer", "Carryover")) %>%
  left_join(ctrlavg_raw_cellcell, by = c("CellCondition", "Timepoint", "Dose", "Interaction_pair", "Replicate")) %>%
  group_by(Protein, CellCondition, Timepoint, Interaction_pair) %>% 
  summarize(
    T_statistic_all3 = t.test(Interaction_score, c(rCD4, `Buffer`, `Carryover`), var.equal = FALSE)$statistic,
    p_raw_all3 = t.test(Interaction_score, c(rCD4, `Buffer`, `Carryover`), var.equal = FALSE)$p.value,
    .groups = "drop_last") %>%
  ungroup() %>%
  distinct() %>%
  filter(!grepl("M16", Interaction_pair)) %>% # population too rare to analyze ; swings from 0% to 100% scores because cell counts too low
  mutate(p_adj_all3 = p.adjust(p_raw_all3, method = "fdr"))

hclust_cellcell_matrix = normavg_cellcell_combdose %>%  # cluster for visualization
  complete(condition_proteinID, Interaction_pair) %>%
  select(condition_proteinID, Interaction_pair, Interaction_score_norm1med3) %>%
  pivot_wider(names_from = Interaction_pair, values_from = Interaction_score_norm1med3)
hclust_cellcell_matrix = as.matrix(data.table::as.data.table(hclust_cellcell_matrix), rownames = 'condition_proteinID')
hclust_cellcell_result = hclust(dist(hclust_cellcell_matrix))
normavg_cellcell_combdose$condition_proteinID = factor(normavg_cellcell_combdose$condition_proteinID, levels = hclust_cellcell_result$labels[hclust_cellcell_result$order])

combined_cellcell_combdose = full_join(filter(normavg_cellcell_combdose, !grepl("M16", Interaction_pair)), pvals_cellcell_combdose,
                                       by = c("Protein", "CellCondition", "Timepoint", "Interaction_pair"))

#### plot of interaction score changes ####
normtype = "med3"
pvaltype = "all3"
normname = paste0("Interaction_score_norm1", normtype)
pvalname = paste0("p_adj_", pvaltype)
combined_cellcell_combdose %>%
  filter(!Protein %in% c("rCD4", "Buffer", "Carryover", "PBS")) %>%
  ggplot() +
  geom_point(aes(x = Interaction_pair, y = conditionID, fill = !!sym(normname), 
                 size = -log10(!!sym(pvalname)), color = !!sym(pvalname) < pthresh), pch = 21) +
  facet_wrap( ~ Protein, ncol = 6) +
  scale_size_continuous(limits = c(0,4.9)) +  # make consistent across all plots
  scale_fill_gradient2(midpoint = 0, low = "#363b81", high = "#a80426", mid = "white", limits = c(-1,1)) +
  scale_color_manual(values = c("white", "black"),name = "Significant") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 8, angle = 90, vjust = 0.5))
#ggsave(paste0("ProtPharm_Trial2_cellcell_combdose_norm1", normtype, "_p_adj", pvaltype, "_circles_p010_noM16_v2b.png"), width = 22, height = 8) #width = 26, height = 8)


#### comparison of experimental perturbation data to model ####

#### combine output files from difeq model ####

mutfile_list = list.files("./difeq_outfiles")
mutfile_list = mutfile_list[grepl("\\.txt$", mutfile_list)&!grepl("^difeq_out\\.", mutfile_list)]

colnames_clean = scan(text = readLines("./difeq_outfiles/difeq_out.txt", 1), what = "", quiet = TRUE)[-1]
mut_tsv = read_delim("./difeq_outfiles/difeq_out.txt", delim = " ", col_names = colnames_clean, skip = 1) %>%
  select(-starts_with("__s")) %>%
  filter(row_number() == n()) %>%
  mutate(condition = "wildtype")

for (mutf in mutfile_list) {  # combining python model outputs
  cond_name = gsub("^difeq_out_", "", gsub(".txt$", "", mutf))
  colnames_clean = scan(text = readLines(paste0("./difeq_outfiles/", mutf), 1), what = "", quiet = TRUE)[-1]  # workaround to get around the # symbol in the file header
  mut_i = read_delim(paste0("./difeq_outfiles/", mutf), delim = " ", col_names = colnames_clean, skip = 1) %>%
    select(-starts_with("__s")) %>%
    filter(row_number() == n()) %>%
    mutate(condition = cond_name)
  mut_tsv = bind_rows(mut_tsv, mut_i)
}

#### normalize cell interaction scores following same calculation procedure as for pharmacoscopy iScores ####

out_tbl = mut_tsv %>%
  group_by(condition) %>%
  mutate(ptotal_T4 = sum(oT4, oT4_T4, oT4_T8, oT4_NK, oT4_DC, oT4_M14, oT4_M16, oT4_B),  # total proportion of T4 cells (% of all cells that are T4 - should equal exact number used as model input)
         ptotal_T8 = sum(oT8, oT4_T8, oT8_T8, oT8_NK, oT8_DC, oT8_M14, oT8_M16, oT8_B),
         ptotal_NK = sum(oNK, oT4_NK, oT8_NK, oNK_NK, oNK_DC, oNK_M14, oNK_M16, oNK_B),
         ptotal_DC = sum(oDC, oT4_DC, oT8_DC, oNK_DC, oDC_DC, oDC_M14, oDC_M16, oDC_B),
         ptotal_M14 = sum(oM14, oT4_M14, oT8_M14, oNK_M14, oDC_M14, oM14_M14, oM14_M16, oM14_B),
         ptotal_M16 = sum(oM16, oT4_M16, oT8_M16, oNK_M16, oDC_M16, oM14_M16, oM16_M16, oM16_B),
         ptotal_B = sum(oB, oT4_B, oT8_B, oNK_B, oDC_B, oM14_B, oM16_B, oB_B),
         ptotal_anybind = 1 - sum(oT4, oT8, oNK, oDC, oM14, oM16, oB)) %>%  # percent of cells bound to another cell
  mutate(fT4_T4 = oT4_T4 / ptotal_T4,  # T4
         fT4_T8 = oT4_T8 / ptotal_T4,
         fT8_T4 = oT4_T8 / ptotal_T8,
         fT4_NK = oT4_NK / ptotal_T4,  # percent of all T4 bound to NK
         fNK_T4 = oT4_NK / ptotal_NK,  # percent of all NK bound to T4
         fT4_DC = oT4_DC / ptotal_T4,
         fDC_T4 = oT4_DC / ptotal_DC,
         fT4_M14 = oT4_M14 / ptotal_T4,
         fM14_T4 = oT4_M14 / ptotal_M14,
         fT4_M16 = oT4_M16 / ptotal_T4,
         fM16_T4 = oT4_M16 / ptotal_M16,
         fT4_B = oT4_B / ptotal_T4,
         fB_T4 = oT4_B / ptotal_B,
         fT8_T8 = oT8_T8 / ptotal_T8,  # T8
         fT8_NK = oT8_NK / ptotal_T8,  
         fNK_T8 = oT8_NK / ptotal_NK,  
         fT8_DC = oT8_DC / ptotal_T8,
         fDC_T8 = oT8_DC / ptotal_DC,
         fT8_M14 = oT8_M14 / ptotal_T8,
         fM14_T8 = oT8_M14 / ptotal_M14,
         fT8_M16 = oT8_M16 / ptotal_T8,
         fM16_T8 = oT8_M16 / ptotal_M16,
         fT8_B = oT8_B / ptotal_T8,
         fB_T8 = oT8_B / ptotal_B,
         fNK_NK = oNK_NK / ptotal_NK,  # NK
         fNK_DC = oNK_DC / ptotal_NK,
         fDC_NK = oNK_DC / ptotal_DC,
         fNK_M14 = oNK_M14 / ptotal_NK,
         fM14_NK = oNK_M14 / ptotal_M14,
         fNK_M16 = oNK_M16 / ptotal_NK,
         fM16_NK = oNK_M16 / ptotal_M16,
         fNK_B = oNK_B / ptotal_NK,
         fB_NK = oNK_B / ptotal_B,
         fDC_DC = oDC_DC / ptotal_DC,   # DC
         fDC_M14 = oDC_M14 / ptotal_DC,
         fM14_DC = oDC_M14 / ptotal_M14,
         fDC_M16 = oDC_M16 / ptotal_DC,
         fM16_DC = oDC_M16 / ptotal_M16,
         fDC_B = oDC_B / ptotal_DC,
         fB_DC = oDC_B / ptotal_B,
         fM14_M14 = oM14_M14 / ptotal_M14,     # M14
         fM14_M16 = oM14_M16 / ptotal_M14,
         fM16_M14 = oM14_M16 / ptotal_M16,
         fM14_B = oM14_B / ptotal_M14,
         fB_M14 = oM14_B / ptotal_B,
         fM16_M16 = oM16_M16 / ptotal_M16,   # M16
         fB_M16 = oM16_B / ptotal_B,
         fM16_B = oM16_B / ptotal_M16,
         fB_B = oB_B / ptotal_B) %>%   # B
  mutate(iT4_T4 = log2(fT4_T4 / (ptotal_T4 * ptotal_T4 * ptotal_anybind)),
         iT4_T8 = log2(fT4_T8 / (ptotal_T4 * ptotal_T8 * ptotal_anybind)),
         iT8_T4 = log2(fT8_T4 / (ptotal_T8 * ptotal_T4 * ptotal_anybind)),
         iT4_NK = log2(fT4_NK / (ptotal_T4 * ptotal_NK * ptotal_anybind)),
         iNK_T4 = log2(fNK_T4 / (ptotal_NK * ptotal_T4 * ptotal_anybind)),
         iT4_DC = log2(fT4_DC / (ptotal_T4 * ptotal_DC * ptotal_anybind)),
         iDC_T4 = log2(fDC_T4 / (ptotal_DC * ptotal_T4 * ptotal_anybind)),
         iT4_M14 = log2(fT4_M14 / (ptotal_T4 * ptotal_M14 * ptotal_anybind)),
         iM14_T4 = log2(fM14_T4 / (ptotal_M14 * ptotal_T4 * ptotal_anybind)),
         iT4_M16 = log2(fT4_M16 / (ptotal_T4 * ptotal_M16 * ptotal_anybind)),
         iM16_T4 = log2(fM16_T4 / (ptotal_M16 * ptotal_T4 * ptotal_anybind)),
         iT4_B = log2(fT4_B / (ptotal_T4 * ptotal_B * ptotal_anybind)),
         iB_T4 = log2(fB_T4 / (ptotal_B * ptotal_T4 * ptotal_anybind)),
         iT8_T8 = log2(fT8_T8 / (ptotal_T8 * ptotal_T8 * ptotal_anybind)),
         iT4_NK = log2(fT4_NK / (ptotal_T8 * ptotal_NK * ptotal_anybind)),
         iNK_T8 = log2(fNK_T8 / (ptotal_NK * ptotal_T8 * ptotal_anybind)),
         iT8_DC = log2(fT8_DC / (ptotal_T8 * ptotal_DC * ptotal_anybind)),
         iDC_T8 = log2(fDC_T8 / (ptotal_DC * ptotal_T8 * ptotal_anybind)),
         iT8_M14 = log2(fT8_M14 / (ptotal_T8 * ptotal_M14 * ptotal_anybind)),
         iM14_T8 = log2(fM14_T8 / (ptotal_M14 * ptotal_T8 * ptotal_anybind)),
         iT8_M16 = log2(fT8_M16 / (ptotal_T8 * ptotal_M16 * ptotal_anybind)),
         iM16_T8 = log2(fM16_T8 / (ptotal_M16 * ptotal_T8 * ptotal_anybind)),
         iT8_B = log2(fT8_B / (ptotal_T8 * ptotal_B * ptotal_anybind)),
         iB_T8 = log2(fB_T8 / (ptotal_B * ptotal_T8 * ptotal_anybind)),
         iNK_NK = log2(fNK_NK / (ptotal_NK * ptotal_NK * ptotal_anybind)),
         iNK_DC = log2(fNK_DC / (ptotal_NK * ptotal_DC * ptotal_anybind)),
         iDC_NK = log2(fDC_NK / (ptotal_DC * ptotal_NK * ptotal_anybind)),
         iNK_M14 = log2(fNK_M14 / (ptotal_NK * ptotal_M14 * ptotal_anybind)),
         iM14_NK = log2(fM14_NK / (ptotal_M14 * ptotal_NK * ptotal_anybind)),
         iNK_M16 = log2(fNK_M16 / (ptotal_NK * ptotal_M16 * ptotal_anybind)),
         iM16_NK = log2(fM16_NK / (ptotal_M16 * ptotal_NK * ptotal_anybind)),
         iNK_B = log2(fNK_B / (ptotal_NK * ptotal_B * ptotal_anybind)),
         iB_NK = log2(fB_NK / (ptotal_B * ptotal_NK * ptotal_anybind)),
         iDC_DC = log2(fDC_DC / (ptotal_DC * ptotal_DC * ptotal_anybind)),
         iDC_M14 = log2(fDC_M14 / (ptotal_DC * ptotal_M14 * ptotal_anybind)),
         iM14_DC = log2(fM14_DC / (ptotal_M14 * ptotal_DC * ptotal_anybind)),
         iDC_M16 = log2(fDC_M16 / (ptotal_DC * ptotal_M16 * ptotal_anybind)),
         iM16_DC = log2(fM16_DC / (ptotal_M16 * ptotal_DC * ptotal_anybind)),
         iDC_B = log2(fDC_B / (ptotal_DC * ptotal_B * ptotal_anybind)),
         iB_DC = log2(fB_DC / (ptotal_B * ptotal_DC * ptotal_anybind)),
         iM14_M14 = log2(fM14_M14 / (ptotal_M14 * ptotal_M14 * ptotal_anybind)), 
         iM14_M16 = log2(fM14_M16 / (ptotal_M14 * ptotal_M16 * ptotal_anybind)),
         iM16_M14 = log2(fM16_M14 / (ptotal_M16 * ptotal_M14 * ptotal_anybind)),
         iM14_B = log2(fM14_B / (ptotal_M14 * ptotal_B * ptotal_anybind)),
         iB_M14 = log2(fB_M14 / (ptotal_B * ptotal_M14 * ptotal_anybind)),
         iM16_M16 = log2(fM16_M16 / (ptotal_M16 * ptotal_M16 * ptotal_anybind)),
         iM16_B = log2(fM16_B / (ptotal_M16 * ptotal_B * ptotal_anybind)),
         iB_M16 = log2(fB_M16 / (ptotal_B * ptotal_M16 * ptotal_anybind)),
         iB_B = log2(fB_B / (ptotal_B * ptotal_B * ptotal_anybind))) %>%
  ungroup()
#write_csv(pivot_wider(select(out_tbl, -proportion), names_from = "cells", values_from = "proportion_scale"), "difeq_v8_mut_combined_propscale_iscore.csv")
#write_csv(pivot_wider(select(out_tbl, -proportion_scale), names_from = "cells", values_from = "proportion"), "difeq_v8_mut_combined_prop_iscore.csv")

#### combine model outputs with pharmacoscopy measurements  ####

map_cellpair_names = c('B-->B' = 'B_B', 'B-->M14' = 'B_M14',
                       'B-->Negs' = 'B_DC', 'B-->NK' = 'B_NK',  
                       'B-->T4' = 'B_T4',
                       'M14-->M14' = 'M14_M14', 'M14-->B' = 'M14_B',
                       'M14-->Negs' = 'M14_DC',
                       'M14-->NK' = 'M14_NK', 'M14-->T4' = 'M14_T4',
                       'Negs-->Negs' = 'DC_DC', 'Negs-->B' = 'DC_B', 
                       'Negs-->M14' = 'DC_M14', 'Negs-->NK' = 'DC_NK', 
                       'Negs-->T4' = 'DC_T4',
                       'NK-->NK' = 'NK_NK', 'NK-->B' = 'NK_B',
                       'NK-->M14' = 'NK_M14', 
                       'NK-->Negs' = 'NK_DC', 'NK-->T4' = 'NK_T4',
                       'T4-->T4' = 'T4_T4', 'T4-->B' = 'T4_B',
                       'T4-->M14' = 'T4_M14', 'T4-->Negs' = 'T4_DC',
                       'T4-->NK' = 'T4_NK')

difeq_model = out_tbl %>%
  pivot_longer(-condition, names_to = "Int_pair_score", values_to = "Interaction_score") %>%
  mutate(Score_type = substr(Int_pair_score, 1, 1)) %>% # first letter indicates type
  mutate(Interaction_pair = ifelse(Score_type != "p", substr(Int_pair_score, 2, nchar(Int_pair_score)), gsub("^.*_", "", Int_pair_score))) %>%
  mutate(Interaction_pair = ifelse(Score_type != "p", plyr::revalue(Interaction_pair, unlist(Biobase::reverseSplit(map_cellpair_names))), Int_pair_score)) %>%  # reverse names and values of the named list
  mutate(Score_type = plyr::mapvalues(Score_type, from = c('o', 'p', 'f', 'i'), to = c('observed proportion', 'total proportion', 'percentage bound', 'iScore'))) %>%
  mutate(Interaction_pair = gsub("_", "-->", Interaction_pair))

normavg_cellcell = raw_cellcell %>% 
  mutate(normgroup = Dose) %>%
  group_by(normgroup, Protein, CellCondition, Timepoint, Interaction_pair) %>%
  summarize_at(vars(Interaction_score), .funs = median) %>% 
  ungroup() %>%
  group_by(normgroup, CellCondition, Timepoint, Interaction_pair) %>%
  mutate(control_median3 = ifelse(normgroup > 10, median(c(Interaction_score[Protein == "Buffer"], Interaction_score[Protein == "Carryover"], Interaction_score[Protein == "rCD4"])), NA),
         control_mean3 = ifelse(normgroup > 10, mean(c(Interaction_score[Protein == "Buffer"], Interaction_score[Protein == "Carryover"], Interaction_score[Protein == "rCD4"])), NA)) %>%
  mutate(Interaction_score_norm1med3 = ifelse(normgroup == 10, ifelse(grepl("\\+", Protein), (Interaction_score - control_ABprotein) / pmax(Interaction_score, control_ABprotein, na.rm = T), (Interaction_score - control_AB) / pmax(Interaction_score, control_AB, na.rm = T)), (Interaction_score - control_median3) / pmax(Interaction_score, control_median3, na.rm = T)),
         Interaction_score_norm1mean3 = ifelse(normgroup == 10, ifelse(grepl("\\+", Protein), (Interaction_score - control_ABprotein) / pmax(Interaction_score, control_ABprotein, na.rm = T), (Interaction_score - control_AB) / pmax(Interaction_score, control_AB, na.rm = T)), (Interaction_score - control_mean3) / pmax(Interaction_score, control_mean3, na.rm = T))) %>%
  ungroup() %>%
  mutate(conditionID = paste(CellCondition, Timepoint)) %>%
  mutate(conditionID = factor(conditionID, levels = c('Resting 4hr', 'Resting 24hr', 'Activated 4hr', 'Activated 24hr'))) %>%
  mutate(condition_proteinID = paste(conditionID, Protein)) %>%
  rename(Interaction_score_raw = Interaction_score)

pvals_cellcell = raw_cellcell %>%
  group_by(Protein, CellCondition, Timepoint, Dose, Interaction_pair) %>%
  mutate(Replicate = row_number()) %>%
  ungroup() %>%
  filter(Dose != 0) %>%
  filter(!Protein %in% c("rCD4", "Buffer", "Carryover")) %>%
  left_join(ctrlavg_raw_cellcell, by = c("CellCondition", "Timepoint", "Dose", "Interaction_pair", "Replicate")) %>%
  group_by(Protein, CellCondition, Timepoint, Interaction_pair) %>% 
  summarize(
    T_statistic_all3 = t.test(Interaction_score, c(rCD4, `Buffer`, `Carryover`), var.equal = FALSE)$statistic,
    p_raw_all3 = t.test(Interaction_score, c(rCD4, `Buffer`, `Carryover`), var.equal = FALSE)$p.value,
    .groups = "drop_last") %>%
  ungroup() %>%
  distinct() %>%
  filter(!grepl("M16", Interaction_pair)) %>% 
  mutate(p_adj_all3 = p.adjust(p_raw_all3, method = "fdr"))

combined_cellcell = full_join(filter(normavg_cellcell, !grepl("M16", Interaction_pair)), pvals_cellcell,
                               by = c("Protein", "CellCondition", "Timepoint", "Interaction_pair"))

combined_cellcell %<>% 
  group_by(Protein, Timepoint, CellCondition) %>%
  mutate(Perturbation_value = abs(Interaction_score_norm1med3)) %>%
  mutate(Perturbation_rank = rank(-Perturbation_value, ties.method = "first")) %>%
  ungroup()

num_total_ranks = difeq_model %>% filter(Score_type == "iScore") %>% filter(!grepl("M16", Interaction_pair)) %>% distinct(Interaction_pair) %>% nrow()

threshold_ranks = round(num_total_ranks * 0.33)  # compare the top third highest predictions
nochange_list = difeq_model %>%
  filter(Score_type == "iScore") %>%
  mutate(Interaction_score_model = (Interaction_score - Interaction_score[condition == "wildtype"]) / pmax(Interaction_score, Interaction_score[condition == "wildtype"], na.rm = T)) %>%
  filter(Interaction_score_model == 0) %>%
  select(condition) %>%
  distinct()
nochange_list = as.vector(nochange_list[[1]])

difeq_pharm_toprank = difeq_model %>%
  filter(Score_type == "iScore") %>%
  filter(!grepl("M16", Interaction_pair)) %>%
  mutate(Interaction_score_model = (Interaction_score - Interaction_score[condition == "wildtype"]) / pmax(Interaction_score, Interaction_score[condition == "wildtype"], na.rm = T)) %>%
  mutate(Interaction_score_model = abs(Interaction_score_model)) %>%
  group_by(condition) %>%
  mutate(condition_max = max(Interaction_score_model)) %>%
  mutate(Score_rank = rank(-Interaction_score_model, ties.method = "first")) %>%
  ungroup() %>%
  filter(condition_max > 2E-4) %>%
  mutate(perturbed_in_model = Score_rank <= threshold_ranks) %>%
  group_by(condition) %>%
  mutate(condition_num_perturbed = sum(perturbed_in_model)) %>%
  ungroup() %>%
  rename(Interaction_score_model_raw = Interaction_score) %>%
  inner_join(combined_cellcell, by = "Interaction_pair") %>%
  mutate(Protein = ifelse(Protein == "SEM4D", "SEMA4D", Protein)) %>%
  filter(CellCondition == "Resting") %>%
  filter(Timepoint == "4hr") %>%
  filter(!condition %in% nochange_list) %>% 
  filter(condition == Protein) %>%
  mutate(Cell1 = gsub("-->.*$", "", Interaction_pair),
         Cell2 = gsub("^.*-->", "", Interaction_pair)) 

difeq_pharm_toprank_pval = difeq_pharm_toprank %>%
  mutate(perturbed_in_model = ifelse(perturbed_in_model, "Predict_perturbed", "Predict_unperturbed")) %>%
  pivot_wider(names_from = perturbed_in_model, values_from = Perturbation_value) %>%
  group_by(condition) %>%
  summarize(t_value = t.test(Predict_perturbed, Predict_unperturbed, var.equal = TRUE, alternative = "greater")$statistic,
            p_raw = t.test(Predict_perturbed, Predict_unperturbed, var.equal = TRUE, alternative = "greater")$p.value,
            .groups = "keep") %>%
  ungroup() %>%
  mutate(p_adj = p.adjust(p_raw, method = "fdr"))

cell_colors = c("B" = "#7EAF9D", "M14" = "#E77F69", "Negs" = "#D196AE", "NK" = "#91A3CC", "T4" = "#4DBCEC", "T8" = "#F2C850")

inner_join(difeq_pharm_toprank, difeq_pharm_toprank_pval, by = "condition") %>%
  ggplot(aes(x = perturbed_in_model, y = Perturbation_value)) +
  geom_boxplot() +
  geom_dotplot(aes(group = Interaction_pair, fill = Cell1, color = Cell2), binaxis='y', stackdir='center', stackratio = 1.1, dotsize = 0.8, stroke = 2.4) + 
  geom_text(data = difeq_pharm_toprank_pval, 
            mapping = aes(x = -Inf, y = +Inf, label = sprintf("p = %.3f\nq = %.3f", p_raw, p_adj)), hjust = -0.1, vjust = +1.2, size = 3.5) +
  scale_fill_manual(values = cell_colors) + 
  scale_color_manual(values = cell_colors) +
  facet_wrap( ~ condition, scales = "free", nrow = 1) +
  theme_bw() + 
  xlab("If cell-cell interaction predicted to be perturbed") +
  ylab("Measured perturbation upon adding recombinant protein")
#ggsave("modelv8_pharmacoscopy_Trial2_protein_KOmatches_toprank_perturb_boxplot_pqvals.png", width = 21, height = 4)
