library('tidyverse')
library('reshape2')

setwd("./organized_code_screen_processing")


#### loading primary screen then median polishing its data ####

data_raw_csv = as.matrix(read.csv("leuk_interactome_180928.csv", row.names = 1, header = T, check.names = F))
data_raw_p1 = data_raw_csv[1:384, ]  # process separate phases of the screen separately, to cover both plate batch effects and time batch effects
data_raw_p2 = data_raw_csv[385:576, ]
data_raw_m = rbind(data_raw_p1, data_raw_p2)

# two-way median polish
data_medpol_full_p1 = medpolish(data_raw_p1, eps = 0.01, trace.iter = F, na.rm = T)
data_medpol_full_p2 = medpolish(data_raw_p2, eps = 0.01, trace.iter = F, na.rm = T)  
data_medpol_p1 = data_medpol_full_p1$residuals  
data_medpol_p2 = data_medpol_full_p2$residuals 
data_medpol_m = rbind(data_medpol_p1, data_medpol_p2)
#write.csv(pmax(data_medpol_m, 0), "Primary_screen_median_polished.csv")

data_tables = reshape2::melt(data_raw_m) %>% # convert matrices to dplyr format
  rename("Bait" = "Var1", "Prey" = "Var2", "Value_raw" = "value") %>%
  left_join(rename(reshape2::melt(data_medpol_m), "Bait" = "Var1", "Prey" = "Var2", "Value_medpol" = "value"), 
            by = c("Bait", "Prey")) %>%
  mutate("Value_medpol_pos" = ifelse(Value_medpol < 0, 0, Value_medpol))
#save(data_tables, file = "data_tables.Rda")

#### same loading and median polish for secondary screen data ####

secondary_screen_raw = as.matrix(read.csv("Leukocyte_2ndary_screen_data.csv", 
                                          row.names = 1, header = T, check.names = F))  # need to avoid checkname conversion of + symbols into . symbols

secondary_screen_medpol_full = medpolish(secondary_screen_raw, eps = 0.01, trace.iter = F, na.rm = T)
secondary_screen_medpol = secondary_screen_medpol_full$residuals  # pull out matrix of residuals
#write.csv(pmax(secondary_screen_medpol, 0), "Secondary_screen_median_polished.csv")

secondary_screen_tbl = reshape2::melt(secondary_screen_raw) %>%
  rename("Bait" = "Var1", "Prey" = "Var2", "Value_raw" = "value") %>%
  left_join(rename(reshape2::melt(secondary_screen_medpol), "Bait" = "Var1", "Prey" = "Var2", "Value_medpol" = "value"), 
            by = c("Bait", "Prey")) %>%
  mutate("Value_medpol_pos" = ifelse(Value_medpol < 0, 0, Value_medpol))
#save(secondary_screen_tbl, file = "secondary_screen_tbl.Rda")


#### combining orientations for each intxn (excluding control wells) ####

control_list_primary = c("pos1", "OX68", "blank", "pos3", "pos1_p2", "OX68_p2", "pos3_p2", "blank_p2", "rCD4_p2")
control_list_secondary = c("pos1", "pos3", "OX68", "empty", "rCD4")
exclude_list = unique(c(control_list_primary, control_list_secondary))

data_filt = filter(data_tables, !(Bait %in% exclude_list))

data_score_label = data_filt %>% # use interaction-specific 'pair_label' to combine orientations of same interaction
  mutate(pair_label = ifelse( as.character(Bait) < as.character(Prey), paste(Bait, "+", Prey), paste(Prey, "+", Bait))) %>%
  mutate(Bait = as.character(Bait)) %>%
  mutate(Prey = as.character(Prey))

data_score_pairgroup = reshape(transform(data_score_label, time=ave(seq_len(nrow(data_score_label)),pair_label,FUN=seq_along)),dir='w',idvar='pair_label')  # pair up the two measured oreintations for each interaction

data_pairgroup_rank = data_score_pairgroup %>%
  rowwise() %>%
  mutate(Sum_medpol = sum(Value_medpol_pos.1, Value_medpol_pos.2, na.rm = T)) %>%
  ungroup() %>%
  mutate(Sum_medpol = ifelse(Bait.1 == Prey.1,  2*Sum_medpol, Sum_medpol)) # double values for homophilic intxns (so on same scale as others)


# same steps for secondary screen data
secondary_screen_filt = filter(secondary_screen_tbl, !(Bait %in% exclude_list)) 
                               
secondary_screen_label = secondary_screen_filt %>%
  mutate(Bait = as.character(Bait)) %>%  # convert to no longer be factor, so can do < comparison on bait names
  mutate(Prey = as.character(Prey)) %>%
  mutate(pair_label = ifelse(Bait < Prey, paste(Bait, "+", Prey), paste(Prey, "+", Bait)))

secondary_screen_label$timevarID = ave(secondary_screen_label$pair_label, as.factor(secondary_screen_label$pair_label), FUN = seq_along)
secondary_screen_label = as.data.frame(secondary_screen_label)
secondary_screen_pairgroup = reshape(data = secondary_screen_label,
                                     idvar = "pair_label", timevar = "timevarID", direction = "wide")
secondary_screen_pairgroup = as_tibble(secondary_screen_pairgroup)

secondary_screen_pairgroup_rank = secondary_screen_pairgroup %>%
  rowwise() %>%
  mutate(Sum_medpol = sum(Value_medpol_pos.1, Value_medpol_pos.2, na.rm = T)) %>%
  ungroup() %>%
  mutate(Sum_medpol = ifelse(Bait.1 == Prey.1,  2*Sum_medpol, Sum_medpol))


#### adding metadata annotations to both sources ####
  
# start with primary screen
leuk_protein_list = read_csv("./protein_metadata/leukocyte_proteins_final.csv") %>%
  mutate(Symbol_unique = ifelse(is.na(Symbol2), Symbol1, paste(Symbol1,"-",Symbol2, sep = ""))) %>%
  mutate(Uniprot_unique = ifelse(is.na(Uniprot2), Uniprot1, paste(Uniprot1,"-",Uniprot2, sep = "")))

leuk_screen_plates = read_csv("./protein_metadata/leukocyte_screen_plates.csv") %>%
  mutate(Symbol_unique = ifelse(is.na(Symbol2), Symbol1, paste(Symbol1,"-",Symbol2, sep = ""))) %>%
  mutate(plate_position = paste(purification_plate,purification_row,purification_col, sep = "."))

expected_csv = read_csv("./protein_metadata/leuk_expected_manual.csv") %>%
  mutate(pair_label = ifelse( as.character(Gene_L) < as.character(Gene_R), paste(Gene_L, "+", Gene_R), paste(Gene_R, "+", Gene_L)))

bradford_data = read_csv("./protein_metadata/bradford_measurements.csv") %>%
  filter(!is.na(Plate)) %>%
  mutate("plate_position" = paste(Plate,Row,Column, sep = "."))

bradford_screen_plates = left_join(bradford_data, leuk_screen_plates, by = "plate_position")  # combine
bradford_screen_leuk = left_join(bradford_screen_plates, leuk_protein_list, by = "Symbol_unique")

expected_csv = read_csv("./protein_metadata/leuk_expected_manual.csv")
expected_list = expected_csv %>%
  mutate(pair_label = ifelse(as.character(Gene_L) < as.character(Gene_R), paste(Gene_L, "+", Gene_R), paste(Gene_R, "+", Gene_L))) 

negatives_csv = read_csv("./protein_metadata/leuk_negatives_manual.csv")  # fyi includes some potentially real but extraordinarily weak nectin family interactions
negatives_list = negatives_csv %>%
  mutate(pair_label = ifelse(as.character(Gene_A) < as.character(Gene_B), paste(Gene_A, "+", Gene_B), paste(Gene_B, "+", Gene_A)))

data_unpaired = data_filt %>%
  mutate(pair_label = ifelse( as.character(Bait) < as.character(Prey), paste(Bait, "+", Prey), paste(Prey, "+", Bait))) %>%
  left_join(expected_list, by = "pair_label") %>%
  mutate(Negative = ifelse(pair_label %in% negatives_list$pair_label, 1, 0)) %>%
  add_column(Bait_U = leuk_protein_list$Uniprot_unique[match(.$Bait, leuk_protein_list$Symbol_unique)], .after = "Bait") %>%  # inserts uniprot for Bait.1
  add_column(Prey_U = leuk_protein_list$Uniprot_unique[match(.$Prey, leuk_protein_list$Symbol_unique)], .after = "Prey") %>%
  add_column(Expected = ifelse(!is.na(.$Gene_L),1,0), .before = "Uniprot_L") %>%
  add_column(Expected_kd = ifelse(!is.na(.$Kd_nM),1,0), .after = "Expected") %>%
  left_join(bradford_screen_plates, by = c("Bait" = "Symbol_unique")) %>%  
  left_join(bradford_screen_plates, by = c("Prey" = "Symbol_unique"), suffix = c(".bait",".prey"))
#save(data_unpaired, file = "data_unpaired.Rda")

data_paired = data_pairgroup_rank %>%
  left_join(expected_list, by = "pair_label") %>%
  mutate(Negative = ifelse(pair_label %in% negatives_list$pair_label, 1, 0)) %>%
  add_column(Bait_U.1 = leuk_protein_list$Uniprot_unique[match(.$Bait.1, leuk_protein_list$Symbol_unique)], .after = "Bait.1") %>%  # inserts uniprot for Bait.1 (first measurement of the pair)
  add_column(Prey_U.1 = leuk_protein_list$Uniprot_unique[match(.$Prey.1, leuk_protein_list$Symbol_unique)], .after = "Prey.1") %>%
  add_column(Bait_U.2 = leuk_protein_list$Uniprot_unique[match(.$Bait.2, leuk_protein_list$Symbol_unique)], .after = "Bait.2") %>%  # inserts uniprot for Bait.2 (second measurement of the pair)
  add_column(Prey_U.2 = leuk_protein_list$Uniprot_unique[match(.$Prey.2, leuk_protein_list$Symbol_unique)], .after = "Prey.2") %>%
  add_column(Expected = ifelse(!is.na(.$Gene_L),1,0), .before = "Uniprot_L") %>%
  add_column(Expected_kd = ifelse(!is.na(.$Kd_nM),1,0), .after = "Expected") %>%
  left_join(select(bradford_screen_plates, Symbol_unique, band_intensity, conc_ngul), by = c("Bait.1" = "Symbol_unique")) %>%  
  left_join(select(bradford_screen_plates, Symbol_unique, band_intensity, conc_ngul), by = c("Prey.1" = "Symbol_unique"), suffix = c(".bait.1",".prey.1")) %>%
  left_join(select(bradford_screen_plates, Symbol_unique, band_intensity, conc_ngul), by = c("Bait.2" = "Symbol_unique")) %>%  
  left_join(select(bradford_screen_plates, Symbol_unique, band_intensity, conc_ngul), by = c("Prey.2" = "Symbol_unique"), suffix = c(".bait.2",".prey.2"))
#save(data_paired, file = "data_paired.Rda")
  

#### repeat same procedure for secondary screen data  ####

leuk_secondary_screen_proteins = read_csv("./protein_metadata/followup_plasmids_FINAL.csv")

secondary_screen_unpaired = secondary_screen_filt %>%
  as_tibble() %>%
  mutate(pair_label = ifelse( as.character(Bait) < as.character(Prey), paste(Bait, "+", Prey), paste(Prey, "+", Bait))) %>%
  left_join(expected_list, by = "pair_label") %>%
  mutate(Negative = ifelse(pair_label %in% negatives_list$pair_label, 1, 0)) %>%
  add_column(Bait_U = leuk_protein_list$Uniprot_unique[match(.$Bait, leuk_protein_list$Symbol_unique)], .after = "Bait") %>%
  add_column(Prey_U = leuk_protein_list$Uniprot_unique[match(.$Prey, leuk_protein_list$Symbol_unique)], .after = "Prey") %>%
  add_column(Expected = ifelse(!is.na(.$Gene_L),1,0), .before = "Uniprot_L") %>%
  add_column(Expected_kd = ifelse(!is.na(.$Kd_nM),1,0), .after = "Expected") %>%
  left_join(leuk_secondary_screen_proteins, by = c("Bait" = "Symbol")) %>%  
  left_join(leuk_secondary_screen_proteins, by = c("Prey" = "Symbol"), suffix = c(".bait",".prey")) %>%
  select(-contains("ul_")) %>%
  rowwise() %>%  
  mutate(expression_min = min(conc_ngul.bait, conc_ngul.prey),
         band_min = min(band_intensity.bait, band_intensity.prey)) %>%
  ungroup()
#save(secondary_screen_unpaired, file = "secondary_screen_unpaired.Rda")

secondary_screen_paired = secondary_screen_pairgroup_rank %>%
  left_join(expected_list, by = "pair_label") %>%
  mutate(Negative = ifelse(pair_label %in% negatives_list$pair_label, 1, 0)) %>%
  add_column(Bait_U.1 = leuk_protein_list$Uniprot_unique[match(.$Bait.1, leuk_protein_list$Symbol_unique)], .after = "Bait.1") %>%
  add_column(Prey_U.1 = leuk_protein_list$Uniprot_unique[match(.$Prey.1, leuk_protein_list$Symbol_unique)], .after = "Prey.1") %>%
  add_column(Bait_U.2 = leuk_protein_list$Uniprot_unique[match(.$Bait.2, leuk_protein_list$Symbol_unique)], .after = "Bait.2") %>%
  add_column(Prey_U.2 = leuk_protein_list$Uniprot_unique[match(.$Prey.2, leuk_protein_list$Symbol_unique)], .after = "Prey.2") %>%
  add_column(Expected = ifelse(!is.na(.$Gene_L),1,0), .before = "Uniprot_L") %>%
  add_column(Expected_kd = ifelse(!is.na(.$Kd_nM),1,0), .after = "Expected") %>%
  left_join(select(leuk_secondary_screen_proteins, Symbol, band_intensity, conc_ngul), by = c("Bait.1" = "Symbol")) %>%  
  left_join(select(leuk_secondary_screen_proteins, Symbol, band_intensity, conc_ngul), by = c("Prey.1" = "Symbol"), suffix = c(".bait.1",".prey.1")) %>%
  left_join(select(leuk_secondary_screen_proteins, Symbol, band_intensity, conc_ngul), by = c("Bait.2" = "Symbol")) %>%  
  left_join(select(leuk_secondary_screen_proteins, Symbol, band_intensity, conc_ngul), by = c("Prey.2" = "Symbol"), suffix = c(".bait.2",".prey.2")) %>%
  rowwise() %>% 
  mutate(expression_min = min(conc_ngul.bait.1, conc_ngul.prey.1),
         band_min = min(band_intensity.bait.1, band_intensity.prey.1)) %>%
  ungroup()
#save(secondary_screen_paired, file = "secondary_screen_paired.Rda")


#### combined score from both primary and secondary ####

leuk_unpaired = data_unpaired %>%
  as_tibble() %>%
  rowwise() %>% 
  mutate(expression_min = min(conc_ngul.bait, conc_ngul.prey),
         band_min = min(band_intensity.bait, band_intensity.prey)) %>%
  ungroup() %>%
  select(-Gene_L, -Gene_R, -Uniprot_L, -Uniprot_R) %>%
  rename_at(.vars = vars(-one_of(c("pair_label", "Bait", "Prey", "Expected_kd", "Kd_nM", "Expected", "Negative")), -contains("Value_")), .funs = list(~paste0(., ".primary"))) %>%
  full_join(secondary_screen_unpaired, 
            by = c("pair_label", "Bait", "Prey"), suffix = c(".primary", ".secondary"), na_matches = 'never') %>%
  mutate(Combined_medpol = ifelse(is.na(Value_medpol_pos.primary), Value_medpol_pos.secondary,
                                    ifelse(is.na(Value_medpol_pos.secondary), Value_medpol_pos.primary,
                                             0.2 * Value_medpol_pos.primary + 0.8 * Value_medpol_pos.secondary)))  # weighted average
#save(leuk_unpaired, file = "leuk_unpaired.Rda")

leuk_paired = data_paired %>%
  distinct() %>%
  rowwise() %>% 
  mutate(expression_min = min(conc_ngul.bait.1, conc_ngul.prey.1),
         band_min = min(band_intensity.bait.1, band_intensity.prey.1)) %>%
  ungroup() %>%
  mutate(Expected_kd = !is.na(Kd_nM)) %>%
  select(-Gene_L, -Gene_R, -Uniprot_L, -Uniprot_R, -Notes) %>%
  rename_at(.vars = vars(-one_of(c("pair_label", "Expected_kd", "Kd_nM", "Expected", "Negative")), -contains("Value_"), -contains("Sum_")), .funs = list(~paste0(., ".primary"))) %>%
  full_join(secondary_screen_paired, 
            by = c("pair_label"), suffix = c(".primary", ".secondary"), na_matches = 'never') %>%
  mutate(Combined_value_medpol.1 = ifelse(is.na(Value_medpol_pos.1.primary), Value_medpol_pos.1.secondary,
                                        ifelse(is.na(Value_medpol_pos.1.secondary), Value_medpol_pos.1.primary,
                                               0.2 * Value_medpol_pos.1.primary + 0.8 * Value_medpol_pos.1.secondary)),
         Combined_value_medpol.2 = ifelse(is.na(Value_medpol_pos.2.primary), Value_medpol_pos.2.secondary,
                                          ifelse(is.na(Value_medpol_pos.2.secondary), Value_medpol_pos.2.primary,
                                                 0.2 * Value_medpol_pos.2.primary + 0.8 * Value_medpol_pos.2.secondary)),
         Combined_sum_medpol = ifelse(is.na(Sum_medpol.primary), Sum_medpol.secondary,
                                      ifelse(is.na(Sum_medpol.secondary), Sum_medpol.primary,
                                             0.2 * Sum_medpol.primary + 0.8 * Sum_medpol.secondary)))
#save(leuk_paired, file = "leuk_paired.Rda")


#### ROC / PRROC analysis ####

library(PRROC)

plot_roc_line <- function(roc_data, tpr_col, fpr_col, cost_col, auc_value, conf_shade = FALSE, lower_shade = NULL, upper_shade = NULL, solid_color = "false", min_alpha = 0.15) {
  
  if (solid_color != "false") {
    #color_ramp = rep(solid_color, times = 100)  # for all one solid color
    color_ramp = sapply(seq(1, min_alpha, length.out = 100), function(a) alpha(solid_color, alpha = a))
    color_ramp = rev(color_ramp)
  }
  else {
    color_ramp = colorRampPalette(c("red","orange","green"))(100)
    #color_ramp = colorRampPalette(c("#d81e33", "#f9dee1"))(100)
    norm_vec = function(v) (v - min(v))/diff(range(v))
    color_by_cost = color_ramp[ceiling(norm_vec((roc_data[[cost_col]]))*99)+1]
  }
  
  auc_label = paste0("AUC = ", formatC(auc_value, digits = 3, format = "f"))
  
  pRoc = ggplot(roc_data, aes_string(x = fpr_col, y = tpr_col)) + 
    {if(conf_shade) geom_ribbon(aes_string(ymin = lower_shade, ymax = upper_shade), alpha = 0.8, fill = "grey")} +
    geom_path(aes_string(color = cost_col), size = 2.5)  +
    scale_color_gradientn(colors = color_ramp, name="Sum_medpol") + 
    #scale_color_gradient2(low = "green", mid = "orange" , high = "red", midpoint=median(roc2$cost)) +
    coord_fixed() + 
    geom_abline(intercept = 0, slope = 1, color=rgb(0,0,1,alpha=0.8), linetype="dotted") + 
    xlab("False Positive Rate") + ylab("True Positive Rate") +
    expand_limits(y=0) +
    theme(panel.background = element_rect(fill = "#f0f0f0")) +
    annotate("label", x = 0.7, y = 0.3, label = auc_label)
    #annotate("label", x = 0.2, y = 0.93, label = auc_label) + theme_bw()
  
  return(pRoc)
}

load("leuk_paired.Rda")   # combined datasets of primary and secondary screens
load("leuk_unpaired.Rda")

if_filter_detectable = TRUE  # consider proteins that were detectably expressed

lectin_stick_list = c("CD209", "CD22", "CLC4M", "CEACAM5", "CD248_duplicate", "CD207")  # omit sugar-binding lectins with many known non-specific interactions and an accidently duplicated bait

if (if_filter_detectable) {
  leuk_paired_expr = filter(leuk_paired, (expression_min.primary > 20)|(band_min.primary > 0),
                          !Bait.1.primary %in% lectin_stick_list, !Bait.2.primary %in% lectin_stick_list)
} else {
  leuk_paired_expr = filter(leuk_paired, !Bait.1.primary %in% lectin_stick_list, !Bait.2.primary %in% lectin_stick_list)  # no expression filter 
}

nonexpect_indexes =  which(leuk_paired_expr$Expected.primary == 0) 
set.seed(10)  # initialize RNG so reproducible
random_ref_500 = sample(nonexpect_indexes, 500)  
random_ref_data = arrange(leuk_paired_expr[random_ref_500, ], desc(Combined_sum_medpol))  # Combined_value_medpol  # for unpaired
leuk_paired_expr = mutate(leuk_paired_expr, random_neg = pair_label %in% random_ref_data$pair_label)

rocset_pos = filter(leuk_paired_expr, Expected.primary == 1)  # Expected_kd.primary   Expected.primary
rocset_neg = filter(leuk_paired_expr, random_neg == 1)  # random_neg    Negative.primary

label = "pair_expressb0_expectedkd_rand500neg_combinedmedpol"

prcurve = pr.curve(scores.class0 = rocset_pos$Combined_sum_medpol,  # Combined_sum_Bscore   Combined_sum_medpol
                   scores.class1 = rocset_neg$Combined_sum_medpol,  # Combined_value_Bscore   Combined_value_medpol
                   curve = T, max.compute = T, min.compute = T, rand.compute = T)
#(paste0("leuk_combined_prroc_",label,".pdf"))
plot(prcurve, rand.plot = TRUE, scale.color = colorRampPalette(c("#f9dee1", "#d81e33"))(20), fill.area = TRUE, fill.color = rgb(0,0,0,0.03), maxminrand.col = "#8c8a8a")
#dev.off()

roccurve = roc.curve(scores.class0 = rocset_pos$Combined_sum_medpol, 
                     scores.class1 = rocset_neg$Combined_sum_medpol, 
                     curve = TRUE)
colnames(roccurve$curve) = c("FPR", "TPR", "Thresh")
plot_roc_line(as.data.frame(roccurve$curve), "TPR", "FPR", auc_value = roccurve$auc,
              solid_color = "#d81e33", conf_shade = F, cost_col = "Thresh", min_alpha = 0.5)
#ggsave(paste0("leuk_combined_roc_",label,".png"), width = 4.5, height = 4)

#### heatmap plot ####

heatmap_matrix <- function(intxn_data, color_fill, bait_col = 'Bait', prey_col = 'Prey', tetcolor_bool = F){
  p_out = ggplot(intxn_data) + geom_tile(aes_string(x=prey_col, y=bait_col, fill=color_fill)) + 
    theme(axis.text.x= element_text(angle=90, hjust=1, vjust = 0.5)) + 
    theme(axis.text.y = element_text(vjust = 0.5)) +
    labs(x = "Prey gene", y = "Bait gene")
  if (tetcolor_bool == 1){
    p_out + scale_fill_manual(values =c("0" = "white", "-1" = "#b5443a", "0.5" = "blue", "1" = "#00c22d"))
  }
  else {
    #p_out + scale_fill_viridis(option="magma", name = "Value")
    p_out + scale_fill_gradient(low = "#636262", high = "#cc1616")
  }
} 
heatmap_at_thresh <- function(intxn_data, thresh_value, value_col, tricolor_expect = FALSE, tetcolor_bool = FALSE) {
  #intxn_data$ranged_value = norm_range(intxn_data[value_col])  # normalized to between 1 and 0
  intxn_data$ranged_value = intxn_data[value_col]
  intxn_data$threshed_value = as.numeric(intxn_data$ranged_value > thresh_value)
  if (tricolor_expect != FALSE) {
    if (tetcolor_bool == TRUE) {
      intxn_data$tetcolor = ifelse(((intxn_data[tricolor_expect] == 0)&(intxn_data$threshed_value == 0)), 0,   # no signal and not expected
                                   ifelse((intxn_data[['pair_label']] %in% hit_list), 0.5,   # novel hits 
                                          #ifelse(((intxn_data[tricolor_expect] == 0)&(intxn_data$threshed_value > 0)), 0.5,  # with signal but not expected
                                          ifelse(((intxn_data[tricolor_expect] == 1)&(intxn_data$threshed_value == 0)), -1,  # no signal but expected
                                                 ifelse(((intxn_data[tricolor_expect] == 1)&(intxn_data$threshed_value == 1)), 1, NA))))  # both signal and expected
      intxn_data$tetcolor_fac = factor(intxn_data$tetcolor)
      heatmap_matrix(intxn_data, color_fill = "tetcolor_fac", tetcolor_bool = TRUE) 
    }
    else {
      intxn_data$threshed_by_known = ifelse(((intxn_data[tricolor_expect] == 0)&(intxn_data$threshed_value > 0)), 0.5, intxn_data$threshed_value) 
      heatmap_matrix(intxn_data, "threshed_by_known") 
    }
  }
  else {
    heatmap_matrix(intxn_data, "threshed_value")
  }
}

hit_list = c("IL1RAP + PTPRF", "CD248 + IL6R", "EFNB1 + EPHA4", "FCER2 + FGFR2",
             "CNTN1 + MCAM", "JAG1 + VASN", "JAG2 + VASN", "EFNB1 + EPHA2",
             "CNTN1 + NRCAM", "PI16 + TNFRSF21", "OLR1 + SLITRK4", 
             "APLP2 + PIGR", "C10orf54 + HLA-F-B2M", "LAIR1 + LAIR1",
             "CDH2 + FCER2", "APLP2 + PLX4", "APP + PLX4", "CD80 + TNR16", 
             "CD93 + IFNGR1", "PLXNB3 + SEMA4D", "AMICA1 + CD320", 
             "CHL1 + L1CAM", "CHL1 + CHL1", "MPZL1 + MPZL1", "C10orf54 + HLA-E-B2M",
             "KIR2DL5A + PVR", "FGFR4 + TGFBR3", "CD1B-B2M + IL6ST")

p_heatmap = ggplot(filter(leuk_unpaired, pair_label %in% hit_list | Expected.primary == 1), 
                   aes(x = Prey, y = Bait)) + 
  geom_tile(aes(alpha = scales::rescale(ifelse(Combined_medpol < 1, 1, Combined_medpol), to=c(0,1)), fill = factor(Expected.primary))) +  # set range for color intensity  
  scale_fill_manual(name = "", breaks = c("1", "0"), labels = c("Expected", "Novel"), values = c("#da4095", "#25783b")) +
  theme(axis.text.x= element_text(angle=90, hjust=1, vjust = 0.5, size = 3.5), # make x labels vertical, centered horizontally
        axis.text.y = element_text(vjust = 0.5, size = 3.5)) +
  labs(x = "Prey gene", y = "Bait gene", alpha = "Value (medpol)")
#ggsave("leuk_combined_unpair_heatmap_expected_or_hitlist.png", p_heatmap, width = 12, height = 11)


