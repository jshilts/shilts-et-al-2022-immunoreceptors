
library(tidyverse)

setwd(r"(C:\Users\jshilts\OneDrive\Documents\Graduate\Lab_Files\Immune_screen_organized\manuscript\leuk_revised2_submission\shilts_custom_code_revised2\organized_code_scRNA_analysis)")

dir_prefix = "../organized_code_scRNA_website"

dataset_list = c("Spleen" = "spleen_cmin10_intxn_leukmanual",  # all tissues
                 "Thymus" = "thymus_cmin10_intxn_leukmanual",
                 "Colon" = "colon_Healthy_cmin10_intxn_leukmanual",
                 "Lung" = "lung_Healthy_cmin10_intxn_leukmanual",
                 "Kidney" = "kidney_Healthy_cmin10_intxn_leukmanual",
                 "Bone Marrow" = "bone_cmin10_intxn_leukmanual",
                 "Lymph Node" = "lymph_cmin10_intxn_leukmanual")

PrintScientific = function(val) {
  exponent = floor(log10(abs(val)))
  base = val * 10^(-1.0*exponent)
  sprintf("%.2fE%+01d", base, exponent)}

#### loop through all tissues ####

myeloid_list = list("Spleen" = c("Monocyte", "DC_2", "Macrophage", "DC_1", "DC_activated"),#, "DC_plasmacytoid"),
                    "Thymus" = c("cDC1", "cDC2", "Macrophage", "mDC", "Monocyte"),#, "pDC"),
                    "Colon" = c("Macrophages", "DC2", "DC1", "Inflammatory Monocytes"),
                    "Lung" = c("Pulmonary dendritic cells", "Luminal macrophages"),
                    "Kidney" = c("MNP1", "MNP2", "MNP3"),#, "Plasmacytoid DC"),
                    "Bone Marrow" = c("Monocyte", "Dendritic cell", "Dendritic/Monocyte", "Macrophage M2", "Macrophage"),#, "Dendritic pDC"),
                    "Lymph Node" = c("Mono", "DC1", "DC2"))#, "pDC"))


data_key_centralities_all = tibble(Celltype = as.character(), Closeness = as.numeric(), Eigencentrality = as.numeric(),
                                   myeloid = as.character(), Tissue = as.character(), Closeness_T_stat = as.numeric(),
                                   Closeness_pval_raw = as.numeric(), Eigencentrality_T_stat = as.numeric(), Eigencentrality_pval_raw = as.numeric(),
                                   Closeness_pval_perm = as.numeric(), Eigencentrality_pval_perm = as.numeric(),
                                   n_intprot_uniq = as.numeric(), Intexpn_T_stat = as.numeric(), Intexpn_pval_raw = as.numeric())
data_permutations_all = tibble(Tissue = as.character(), 
                               Closeness_perm = as.numeric(), Eigencentrality_perm = as.numeric(),
                               Closeness_obs = as.numeric(), Eigencentrality_obs = as.numeric(),
                               Closeness_pval_perm = as.numeric(), Eigencentrality_pval_perm = as.numeric())

for (dataset_i in names(myeloid_list)) {
  
  fname = dataset_list[dataset_i] 
  data_key = read_csv(paste0(dir_prefix, "/circos_", fname, "_QKEY.csv"))
  
  data_key$Expected = 1
  data_key = mutate(data_key, pair_cells = ifelse(chr1 < chr2, paste0(chr1, " + ", chr2), paste0(chr2, " + ", chr1))) 
  data_key = mutate(data_key, bin_condition = Expected)
  
  user_threshold = 10 / 100
  combined_data_key = filter(data_key, pct_expr_1 > user_threshold & pct_expr_2 > user_threshold)
  
  data_key_summary = combined_data_key %>%
    group_by(chr1, chr2) %>%
    summarize(n_expect = sum(Expected==1), .groups = "drop_last") %>%
    ungroup()
  
  data_key_igraph = igraph::graph_from_edgelist(as.matrix(select(data_key_summary, -n_expect)))
  
  closeness_centralities = igraph::closeness(data_key_igraph, weights = data_key_summary$n_expect, mode = "all")
  eigen_centralities = igraph::eigen_centrality(data_key_igraph, directed = FALSE, 
                                                weights = data_key_summary$n_expect)
  
  myeloid_pops = myeloid_list[[dataset_i]]
  
  data_key_centralities = tibble(Celltype = names(closeness_centralities),
                                 Closeness = unname(closeness_centralities),
                                 Eigencentrality = eigen_centralities$vector) %>% 
    mutate(myeloid = factor(Celltype %in% myeloid_pops, levels = c(TRUE, FALSE), labels = c("Myeloid cells", "Other cells")))
  
  closeness_ttest = t.test(filter(data_key_centralities, myeloid == "Other cells")$Closeness,
                           filter(data_key_centralities, myeloid == "Myeloid cells")$Closeness,
                           var.equal = FALSE, alternative = "two.sided")
  eigen_ttest = t.test(filter(data_key_centralities, myeloid == "Other cells")$Eigencentrality,
                       filter(data_key_centralities, myeloid == "Myeloid cells")$Eigencentrality,
                       var.equal = FALSE, alternative = "two.sided")
  
  closeness_mean_dif = mean(filter(data_key_centralities, myeloid == "Myeloid cells")$Closeness) - mean(filter(data_key_centralities, myeloid == "Other cells")$Closeness)
  eigen_mean_dif = mean(filter(data_key_centralities, myeloid == "Myeloid cells")$Eigencentrality) - mean(filter(data_key_centralities, myeloid == "Other cells")$Eigencentrality)
  
  n_permutations = 1000
  permtest_closeness_rand_list = c()
  permtest_eigen_rand_list = c()
  num_myl_to_pick = length(myeloid_pops)
  check_max_permutes = choose(length(names(closeness_centralities)), num_myl_to_pick)
  #print(paste("For", dataset_i, "maxpermute =", check_max_permutes))
  #if (n_permutations > check_max_permutes) {
  #  print(paste("Too many permutations selected for :", dataset_i))
  #}
  for (permute_i in seq(to = n_permutations)) {
    rand_myl = sample(names(closeness_centralities), num_myl_to_pick)
    if (all(rand_myl %in% myeloid_pops)) {
      next
    } 
    data_key_centralities_rand = mutate(data_key_centralities,
                                        myeloid_rand = factor(Celltype %in% rand_myl, levels = c(TRUE, FALSE), labels = c("Myeloid cells", "Other cells")))
    
    closeness_mean_dif_rand = mean(filter(data_key_centralities_rand, myeloid_rand == "Myeloid cells")$Closeness) - mean(filter(data_key_centralities_rand, myeloid_rand == "Other cells")$Closeness)
    eigen_mean_dif_rand = mean(filter(data_key_centralities_rand, myeloid_rand == "Myeloid cells")$Eigencentrality) - mean(filter(data_key_centralities_rand, myeloid_rand == "Other cells")$Eigencentrality)
    permtest_closeness_rand_list = c(permtest_closeness_rand_list, closeness_mean_dif_rand)
    permtest_eigen_rand_list = c(permtest_eigen_rand_list, eigen_mean_dif_rand)
  }
  permute_closeness_pval = sum(closeness_mean_dif > permtest_closeness_rand_list) / n_permutations
  permute_eigen_pval = sum(eigen_mean_dif < permtest_eigen_rand_list) / n_permutations
  
  data_key_centralities = data_key_centralities %>%
    mutate(Tissue = print(dataset_i)) %>%
    mutate(Closeness_T_stat = closeness_ttest$statistic,
           Closeness_pval_raw = closeness_ttest$p.value,
           Eigencentrality_T_stat = eigen_ttest$statistic,
           Eigencentrality_pval_raw = eigen_ttest$p.value,
           Closeness_pval_perm = permute_closeness_pval,
           Eigencentrality_pval_perm = permute_eigen_pval)
  
  data_key_intexpn = data_key %>% 
    select(chr1, Uniprot_1, gene_name_1) %>%
    bind_rows(rename_with(select(data_key, chr2, Uniprot_2, gene_name_2), ~ gsub("2", "1", .x))) %>%
    distinct() %>%
    group_by(chr1) %>%
    summarize(n_intprot_uniq = n(), .groups = "drop_last") %>% 
    ungroup() %>%
    mutate(myeloid = ifelse(chr1 %in% myeloid_pops, "Myeloid cells", "Other cells"))
  
  intexpn_ttest = t.test(filter(data_key_intexpn, myeloid == "Other cells")$n_intprot_uniq,
                         filter(data_key_intexpn, myeloid == "Myeloid cells")$n_intprot_uniq,
                         var.equal = FALSE, alternative = "two.sided")
  
  data_key_intexpn = data_key_intexpn %>%
    rename(Celltype = chr1) %>%
    mutate(Intexpn_T_stat = intexpn_ttest$statistic,
           Intexpn_pval_raw = intexpn_ttest$p.value)
  
  data_key_centralities = data_key_centralities %>%
    left_join(select(data_key_intexpn, - myeloid), by = "Celltype")
  
  data_permutations = tibble(Tissue = dataset_i, 
                             Closeness_perm = permtest_closeness_rand_list, Eigencentrality_perm = permtest_eigen_rand_list,
                             Closeness_obs = closeness_mean_dif, Eigencentrality_obs = eigen_mean_dif,
                             Closeness_pval_perm = permute_closeness_pval,
                             Eigencentrality_pval_perm = permute_eigen_pval)
  
  data_permutations_all = bind_rows(data_permutations_all,
                                    data_permutations)
  data_key_centralities_all = bind_rows(data_key_centralities_all,
                                        data_key_centralities)
}
#write_csv(data_key_centralities_all, "data_key_centralities_all.csv")

data_key_centralities_adj = data_key_centralities_all %>%
  mutate(Closeness_pval_adj = p.adjust(Closeness_pval_raw, method = "BH"),
         Eigencentrality_pval_adj = p.adjust(Eigencentrality_pval_raw, method = "BH"),
         Intexpn_pval_adj = p.adjust(Intexpn_pval_raw, method = "BH"))
           
ggplot(data_key_centralities_adj,
       aes(x = myeloid, y = Eigencentrality, fill = myeloid)) +
  geom_boxplot(alpha = 0.5) +
  geom_dotplot(binaxis='y', stackdir='center', stackratio = 1.1, dotsize = 0.9) +
  geom_text(data = distinct(data_key_centralities_adj, Tissue, .keep_all = T),  # single row for single text label
            mapping = aes(x = -Inf, y = 0, label = sprintf("p = %.3f\nq = %.3f", Eigencentrality_pval_raw, Eigencentrality_pval_adj)),
            color = "black", hjust = -0.4, vjust = -0.2, size = 2.5) +
  theme_bw() +
  facet_wrap( ~ Tissue, scales = "free_y", nrow = 2) +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("")
ggsave(paste0("Centralities_allcond_myeloid_eigencent_boxplot.png"), width = 10, height = 5)

#### permutation test distributions to compare ####

ggplot(data_permutations_all) +
  geom_histogram(aes(x = Eigencentrality_perm, y = ..count../sum(..count..),
                     fill = Eigencentrality_perm > Eigencentrality_obs), alpha = 0.8) +
  geom_text(data = distinct(data_permutations_all, Tissue, .keep_all = T),
            aes(label = sprintf("p = %.3f", Eigencentrality_pval_perm)),
            x = -0.225, y = 0.0175) + 
  geom_vline(aes(xintercept = Eigencentrality_obs), linetype = "longdash") +
  scale_fill_manual(values = c(`TRUE` = "#8a6a53", `FALSE` = "#666b73"), guide = FALSE) +
  theme_bw() + ylab("Frequency in permuted distribution") + xlab("Mean difference in eigenvector centralities") +
  facet_wrap( ~ Tissue, nrow = 2)
ggsave(paste0("Centralities_allcond_myeloid_eigencent_permutations.png"), width = 10, height = 5)


#### control to validate hub centrality is not just artefact ####

ggplot(data_key_centralities_adj,
       aes(x = myeloid, y = n_intprot_uniq)) +
  geom_boxplot() +
  geom_dotplot(binaxis='y', stackdir='center') +
  geom_text(data = distinct(data_key_centralities_adj, Tissue, .keep_all = T),  # single row for single text label
            mapping = aes(x = -Inf, y = 0, label = sprintf("p = %.3f\nq = %.3f", Intexpn_pval_raw, Intexpn_pval_adj)),
            color = "black", hjust = -0.4, vjust = -0.2, size = 2.5) +
  facet_wrap( ~ Tissue, scales = "free_y", nrow = 2) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme_classic() + ylab("Number of expressed interacting surface proteins") + xlab("")
ggsave(paste0("Centralities_allcond_myeloid_intreceptorcount_boxplot.png"), width = 10, height = 5)

