# -*- coding: utf-8 -*-

import scanpy as sc
import numpy as np
import scipy as sp
import pandas as pd
import matplotlib.pyplot as plt
from cycler import cycler as Cycler
import seaborn as sb
from os import chdir as os_chdir

import SpatialDE

dir_path = "/organized_code_spatialRNA"
os_chdir(dir_path)
sc.settings.figdir = dir_path
sc.set_figure_params(dpi_save = 300)

#### loading data from local files ####

adata = sc.read_visium("visium_folder")
adata.var_names_make_unique()
sc.pp.calculate_qc_metrics(adata, inplace = True, percent_top = None)

#### QC processing (following visium tutorial guidelines) ####

adata.var['mt'] = [gene.startswith('MT-') for gene in adata.var_names]
adata.obs['mt_frac'] = adata[:, adata.var['mt']].X.sum(1).A.squeeze()/adata.obs['total_counts']

print(f'Number of cells before filtering: {adata.n_obs}')
sc.pp.filter_cells(adata, min_counts = 4000)
print(f'Number of cells after min count filter: {adata.n_obs}')
sc.pp.filter_cells(adata, max_counts = 36000)
print(f'Number of cells after max count filter: {adata.n_obs}')
adata = adata[adata.obs['mt_frac'] < 0.2]
print(f'Number of cells after MT filter: {adata.n_obs}')
sc.pp.filter_cells(adata, min_genes = 2000)
print(f'Number of cells after gene filter: {adata.n_obs}')
sc.pp.filter_genes(adata, min_cells=5)
print(f'Number of genes after cell filter: {adata.n_vars}')

#### normalizing counts and umaps ####

sc.pp.normalize_total(adata, inplace = True)
sc.pp.log1p(adata)

sc.pp.highly_variable_genes(adata, flavor='seurat', n_top_genes=2000, inplace=True)
sc.pp.pca(adata, n_comps=50, use_highly_variable=True, svd_solver='arpack')
sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.tl.leiden(adata, key_added='clusters')

#### spatial plots showing the histology image ####

sc.pl.spatial(adata,  img_key = "hires", color=['total_counts', 'n_genes_by_counts'], alpha = 0.9)

gene1 = "JAG1"
gene2 = "VASN"
adata_i = adata.copy()
adata_i.obs['gene_bin'] = pd.Categorical(np.array((1*(adata_i[:, gene1].X > 0) + 2*(adata_i[:, gene2].X > 0)).todense()).flatten()) # 0 = neither gene, 1 = gene 1 only, 2 = gene 2 only, 3 = both
adata_i = adata_i[adata_i.obs['gene_bin'].isin([1, 2, 3]), : ]
adata_i.obs['gene_bin'] = adata_i.obs['gene_bin'].cat.rename_categories({1:gene1, 2:gene2, 3:'Both'})
color_pallete_dict = {gene1:"#2eb345", gene2:"#46a4cf", "Both":"#cc9b47"}
adata_i.uns['gene_bin_colors'] = [color_pallete_dict[x] for x in adata_i.obs['gene_bin'].cat.categories.tolist()]
sc.pl.spatial(adata_i,  img_key = "hires", color=['gene_bin'], alpha = 0.9, title = "%s + %s"%(gene1, gene2), crop_coord = [200,1700,1900,150],
              save = "_lymphnode_detected_%s_%s.pdf"%(gene1, gene2), return_fig = False, show = False)


####  co-localization calculations ####

if_full_expect = False      # reported intxns (True) or new list (False)
if_random_intxns = False    # replaces true intxn network pairs with randomized ones to look at
if if_full_expect:
    interaction_network = pd.read_csv("leuk_expected_manual_nolec.csv")
    intxn_network = interaction_network.copy()
else:
    interaction_network = pd.read_csv("leuk_intxns_nohomophil.csv")
    intxn_network = interaction_network.loc[interaction_network.loc[:, 'source'] == "shilts et al", :]
    intxn_network.reset_index(inplace = True)
#end

radius_threshold = 150   # empirical distance between directly neighboring points is 138.34
distmat_num = sp.spatial.distance.cdist(adata.obsm['X_spatial'], adata.obsm['X_spatial'])
distmat_num = pd.DataFrame(distmat_num, columns = [";".join(map(str, x)) for x in adata.obsm['X_spatial'].tolist()], index = [";".join(map(str, x)) for x in adata.obsm['X_spatial'].tolist()])
distmat_bool = (distmat_num < radius_threshold) & (distmat_num > 0)


neighbor_scores = pd.DataFrame()
last_progress = 0
missing_intxns = []
for row_i in range(intxn_network.shape[0]):
    gene1 = intxn_network.loc[row_i, 'Symbol_1']
    gene2 = intxn_network.loc[row_i, 'Symbol_2']
    
    if not(adata.var_names.str.match("^"+gene1+"$").any()):
        missing_intxns.append(gene1+" + "+gene2)
        print("Missing gene: ", gene1)
        continue
    if not(adata.var_names.str.match("^"+gene2+"$").any()):
        missing_intxns.append(gene1+" + "+gene2)
        print("Missing gene: ", gene1)
        continue
    #end
    if gene1 < gene2:
        intxn_pair = gene1 + " + " + gene2
    else:
        intxn_pair = gene2 + " + " + gene1
    #end
    
    if (row_i+1) / intxn_network.shape[0] > last_progress + 0.1:
        last_progress += 0.1
        print("...%d%%"%(last_progress*100))
    #end
    
    spots_gene1 = np.asarray((adata[:, gene1].X > 0).todense().astype(int)).squeeze()
    spots_gene2 = np.asarray((adata[:, gene2].X > 0).todense().astype(int)).squeeze()
    spots_1not2 = spots_gene1 * (~spots_gene2.astype(bool)).astype(int)
    spots_2not1 = spots_gene2 * (~spots_gene1.astype(bool)).astype(int)
    spots_both = spots_gene1 * spots_gene2
    spots_neither = (~(spots_gene1 + spots_gene2).astype(bool)).astype(int)
    
    dist_pair_12 = spots_gene1 * distmat_bool * spots_gene2.T
    dist_pair_21 = spots_gene2 * distmat_bool * spots_gene1.T
    neighbor_pair_count = np.sum(np.sum(dist_pair_12) +
                                 np.sum(dist_pair_21))
    spot_both_count = np.sum(spots_both)
    spot_neither_count = np.sum(spots_neither)
    neighbor_pairnotboth_count = np.sum(np.sum(distmat_bool * spots_gene1 * spots_2not1.T) +
                                        np.sum(distmat_bool * spots_gene2 * spots_1not2.T))
    neighbor_0_count = np.sum(np.sum(spots_gene1 * distmat_bool * (~spots_gene2.astype(bool)).astype(int)) +
                              np.sum(spots_gene2 * distmat_bool * (~spots_gene1.astype(bool)).astype(int)))
    neighbor_11_count = np.sum(np.sum(spots_gene1 * distmat_bool * spots_1not2.T))
    neighbor_22_count = np.sum(np.sum(spots_gene2 * distmat_bool * spots_2not1.T)) 
    neighbor_total_all = neighbor_pair_count + neighbor_0_count + neighbor_11_count + neighbor_22_count 
    neighbor_total_no0 = neighbor_pair_count + neighbor_11_count + neighbor_22_count
    
    
    spots_gene1_ninteract = np.sum(dist_pair_12)[spots_gene1.astype(bool)]
    spots_gene2_ninteract = np.sum(dist_pair_21)[spots_gene2.astype(bool)]
    spot_count_1unpaired = np.sum(spots_gene1_ninteract == 0)
    spot_count_2unpaired = np.sum(spots_gene2_ninteract == 0) 
    spot_count_paired = np.sum(spots_gene1_ninteract > 0)
    
    # check sums
    np.sum(spots_1not2) + np.sum(spots_2not1) + spot_both_count + spot_neither_count == distmat_bool.shape[0]  # all spots have identity assigned into the 4 types
    spot_count_1unpaired + spot_count_2unpaired + spot_count_paired + spot_neither_count == distmat_bool.shape[0]  # all spots are either paired, unpaired, or have no protein
    
    new_row = pd.DataFrame({'gene1':gene1, 'gene2':gene2, 'intxn_pair':intxn_pair,
                            'nspots_1':np.sum(spots_gene1), 'nspots_2':np.sum(spots_gene2),
                            'nspots_1only':np.sum(spots_1not2), 'nspots_2only':np.sum(spots_2not1),
                            'nspots_both':spot_both_count, 'nspots_neither':np.sum((spots_gene1 + spots_gene2) == 0),
                            'nneighbor_pair':neighbor_pair_count, 
                            'pct_neighbors_paired_excl0': neighbor_pair_count / neighbor_total_no0}, index = [row_i])
    
    neighbor_scores = neighbor_scores.append(new_row)
#end        
#neighbor_scores.to_csv("lymph_leuk_neighbor_scores.csv", index = False)


# permutation : keep same actual list of surface proteins but randomize which are paired together
n_permutations = 100
np.random.seed(0)  # set seed
if_leuk_comb_network = True  # considering both surface proteins in the 'known' list and new interaction list
if if_leuk_comb_network:
    interaction_network = interaction_network = pd.read_csv("leuk_intxns_nohomophil.csv")
    interaction_network = interaction_network.loc[interaction_network.loc[:, 'source'] == "shilts et al", :]
    interaction_network = interaction_network[['Symbol_1', 'Symbol_2', 'Uniprot_1', 'Uniprot_2', 'source']]
    intxn_known = pd.read_csv("leuk_expected_manual_nolec.csv")
    intxn_known = intxn_known.dropna(subset = ['Symbol_1', 'Symbol_2'])
    intxn_known.reset_index(inplace = True)
    intxn_known = intxn_known[['Symbol_1', 'Symbol_2', 'Uniprot_L', 'Uniprot_R']]
    intxn_known['source'] = 'leuk known curated'
    intxn_known.rename(columns = {'Uniprot_L':'Uniprot_1', 'Uniprot_R':'Uniprot_2'}, inplace = True)
    interaction_network = interaction_network.append(intxn_known)    
    interaction_network.reset_index(inplace = True)
#end

neighbor_scores_rand = pd.DataFrame()
last_progress = 0
for row_i in range(n_permutations):
    
    if (row_i+1) / n_permutations > last_progress + 0.1:
        last_progress += 0.1
        print("...%d%%"%(last_progress*100))
    #end
    
    rand1 = np.random.randint(interaction_network.shape[0], size = 1).item()
    rand2 = np.random.randint(interaction_network.shape[0], size = 1).item()  #.item to return integer instead of an array of 1 number
    while (rand1 == rand2)|(interaction_network.loc[rand1, 'Symbol_2'] == interaction_network.loc[rand2, 'Symbol_2'])|(interaction_network.loc[rand2, 'Symbol_1'] == interaction_network.loc[rand1, 'Symbol_1']): # remove homophilics (which will have perfect 'colocalization' scores always)
        rand2 = np.random.randint(interaction_network.shape[0], size = 1).item()
    #end
        
    gene1 = interaction_network.loc[rand1, 'Symbol_1']
    gene2 = interaction_network.loc[rand2, 'Symbol_2']
    intxn_pair = gene1 + " + " + gene2
    
    spots_gene1 = np.zeros(adata[:, 0].X.shape[0], dtype = int)
    spots_gene1[np.asarray((adata[:, rand1].X > 0).todense()).flatten()] = 1
    
    spots_gene2 = np.zeros(adata[:, 0].X.shape[0], dtype = int)
    spots_gene2[np.asarray((adata[:, rand2].X > 0).todense()).flatten()] = 1

    spots_1not2 = spots_gene1 * (~spots_gene2.astype(bool)).astype(int)
    spots_2not1 = spots_gene2 * (~spots_gene1.astype(bool)).astype(int)
    spots_both = spots_gene1 * spots_gene2
    spots_neither = (~(spots_gene1 + spots_gene2).astype(bool)).astype(int)
    
    dist_pair_12 = spots_gene1 * distmat_bool * spots_gene2.T
    dist_pair_21 = spots_gene2 * distmat_bool * spots_gene1.T
    neighbor_pair_count = np.sum(np.sum(dist_pair_12) +
                                 np.sum(dist_pair_21))
    spot_both_count = np.sum(spots_both)  
    spot_neither_count = np.sum(spots_neither)
    neighbor_pairnotboth_count = np.sum(np.sum(distmat_bool * spots_gene1 * spots_2not1.T) +
                                        np.sum(distmat_bool * spots_gene2 * spots_1not2.T))
    neighbor_0_count = np.sum(np.sum(spots_gene1 * distmat_bool * (~spots_gene2.astype(bool)).astype(int)) +
                              np.sum(spots_gene2 * distmat_bool * (~spots_gene1.astype(bool)).astype(int)))
    neighbor_11_count = np.sum(np.sum(spots_gene1 * distmat_bool * spots_1not2.T))
    neighbor_22_count = np.sum(np.sum(spots_gene2 * distmat_bool * spots_2not1.T))   
    neighbor_total_all = neighbor_pair_count + neighbor_0_count + neighbor_11_count + neighbor_22_count 
    neighbor_total_no0 = neighbor_pair_count + neighbor_11_count + neighbor_22_count
    
    
    spots_gene1_ninteract = np.sum(dist_pair_12)[spots_gene1.astype(bool)]
    spots_gene2_ninteract = np.sum(dist_pair_21)[spots_gene2.astype(bool)]
    spot_count_1unpaired = np.sum(spots_gene1_ninteract == 0)
    spot_count_2unpaired = np.sum(spots_gene2_ninteract == 0) 
    spot_count_paired = np.sum(spots_gene1_ninteract > 0)
    
    # check sums
    np.sum(spots_1not2) + np.sum(spots_2not1) + spot_both_count + spot_neither_count == distmat_bool.shape[0]  # all spots have identity assigned into the 4 types
    spot_count_1unpaired + spot_count_2unpaired + spot_count_paired + spot_neither_count == distmat_bool.shape[0]  # all spots are either paired, unpaired, or have no protein
    
    new_row = pd.DataFrame({'gene1':gene1, 'gene2':gene2, 'intxn_pair':intxn_pair,
                            'nspots_1':np.sum(spots_gene1), 'nspots_2':np.sum(spots_gene2),
                            'nspots_1only':np.sum(spots_1not2), 'nspots_2only':np.sum(spots_2not1),
                            'nspots_both':spot_both_count, 'nspots_neither':np.sum((spots_gene1 + spots_gene2) == 0),
                            'nneighbor_pair':neighbor_pair_count, 
                            'pct_neighbors_paired_excl0': neighbor_pair_count / neighbor_total_no0}, index = [row_i])

    neighbor_scores_rand = neighbor_scores_rand.append(new_row)
#end        
#neighbor_scores_rand.to_csv("lymph_randpairleuks100_neighbor_scores.csv", index = False)


intxn_known = pd.read_csv("leuk_expected_manual_nolec.csv")
intxn_known = intxn_known.dropna(subset = ['Symbol_1', 'Symbol_2'])
intxn_known.reset_index(inplace = True)
n_sampled = intxn_known.shape[0]
np.random.seed(0)
sampled_rows = np.random.choice(intxn_known.shape[0], n_sampled, replace = False)

neighbor_scores_knowns = pd.DataFrame()  
last_progress = 0
missing_intxns = []
for row_i in sampled_rows:
    
    gene1 = intxn_known.loc[row_i, 'Symbol_1']
    gene2 = intxn_known.loc[row_i, 'Symbol_2']
    
    if not(adata.var_names.str.match("^"+gene1+"$").any()):
        missing_intxns.append(gene1+" + "+gene2)
        print("Missing gene: ", gene1)
        continue
    if not(adata.var_names.str.match("^"+gene2+"$").any()):
        missing_intxns.append(gene1+" + "+gene2)
        print("Missing gene: ", gene1)
        continue
    #end
    if gene1 < gene2:
        intxn_pair = gene1 + " + " + gene2
    else:
        intxn_pair = gene2 + " + " + gene1
    #end
        
    if gene1 == gene2:
        print("Homophilic skipped: ", intxn_pair)
        continue
    #end
    
    if np.where(sampled_rows == row_i)[0].astype(int) / n_sampled > last_progress + 0.1:
        last_progress += 0.1
        print("...%d%%"%(last_progress*100))
    #end
    
    spots_gene1 = np.asarray((adata[:, gene1].X > 0).todense().astype(int)).squeeze()
    spots_gene2 = np.asarray((adata[:, gene2].X > 0).todense().astype(int)).squeeze()
    spots_1not2 = spots_gene1 * (~spots_gene2.astype(bool)).astype(int)
    spots_2not1 = spots_gene2 * (~spots_gene1.astype(bool)).astype(int)
    spots_both = spots_gene1 * spots_gene2
    spots_neither = (~(spots_gene1 + spots_gene2).astype(bool)).astype(int)
    
    dist_pair_12 = spots_gene1 * distmat_bool * spots_gene2.T
    dist_pair_21 = spots_gene2 * distmat_bool * spots_gene1.T
    neighbor_pair_count = np.sum(np.sum(dist_pair_12) +
                                 np.sum(dist_pair_21))
    spot_both_count = np.sum(spots_both)
    spot_neither_count = np.sum(spots_neither)
    neighbor_pairnotboth_count = np.sum(np.sum(distmat_bool * spots_gene1 * spots_2not1.T) +
                                        np.sum(distmat_bool * spots_gene2 * spots_1not2.T))
    neighbor_0_count = np.sum(np.sum(spots_gene1 * distmat_bool * (~spots_gene2.astype(bool)).astype(int)) +
                              np.sum(spots_gene2 * distmat_bool * (~spots_gene1.astype(bool)).astype(int)))
    neighbor_11_count = np.sum(np.sum(spots_gene1 * distmat_bool * spots_1not2.T))
    neighbor_22_count = np.sum(np.sum(spots_gene2 * distmat_bool * spots_2not1.T)) 
    neighbor_total_all = neighbor_pair_count + neighbor_0_count + neighbor_11_count + neighbor_22_count 
    neighbor_total_no0 = neighbor_pair_count + neighbor_11_count + neighbor_22_count
    
    
    spots_gene1_ninteract = np.sum(dist_pair_12)[spots_gene1.astype(bool)]
    spots_gene2_ninteract = np.sum(dist_pair_21)[spots_gene2.astype(bool)]
    spot_count_1unpaired = np.sum(spots_gene1_ninteract == 0)
    spot_count_2unpaired = np.sum(spots_gene2_ninteract == 0) 
    spot_count_paired = np.sum(spots_gene1_ninteract > 0)
    
    # check sums
    np.sum(spots_1not2) + np.sum(spots_2not1) + spot_both_count + spot_neither_count == distmat_bool.shape[0]  # all spots have identity assigned into the 4 types
    spot_count_1unpaired + spot_count_2unpaired + spot_count_paired + spot_neither_count == distmat_bool.shape[0]  # all spots are either paired, unpaired, or have no protein
    
    new_row = pd.DataFrame({'gene1':gene1, 'gene2':gene2, 'intxn_pair':intxn_pair,
                            'nspots_1':np.sum(spots_gene1), 'nspots_2':np.sum(spots_gene2),
                            'nspots_1only':np.sum(spots_1not2), 'nspots_2only':np.sum(spots_2not1),
                            'nspots_both':spot_both_count, 'nspots_neither':np.sum((spots_gene1 + spots_gene2) == 0),
                            'nneighbor_pair':neighbor_pair_count,
                            'pct_neighbors_paired_excl0': neighbor_pair_count / neighbor_total_no0}, index = [row_i])
    
    neighbor_scores_knowns = neighbor_scores_knowns.append(new_row)
#end        
#neighbor_scores_knowns.to_csv("lymph_leukknown_neighbor_scores.csv", index = False)


sb.boxplot(data = pd.melt(pd.DataFrame({'LeukKnown':neighbor_scores_knowns['pct_neighbors_paired_excl0'], 'LeukNew':neighbor_scores['pct_neighbors_paired_excl0'], 'RandomPairs':neighbor_scores_rand['pct_neighbors_paired_excl0']}), value_vars = ['LeukKnown', 'LeukNew','RandomPairs'], value_name = 'percent neighbors interacting', var_name = 'interaction network').dropna(),
             x = 'interaction network', y = 'percent neighbors interacting')  # also looks good as boxplot
plt.savefig("lymph_leuk_vs_known_vs_random_neighbor_scores_boxplot.svg")

from statsmodels.stats.multicomp import pairwise_tukeyhsd

sp.stats.f_oneway(neighbor_scores_rand['pct_neighbors_paired_excl0'],
                  neighbor_scores_knowns['pct_neighbors_paired_excl0'],
                  neighbor_scores['pct_neighbors_paired_excl0'])

neighbor_scores_combined = pd.DataFrame({'dataset':'RandomPairs', 'scores':neighbor_scores_rand['pct_neighbors_paired_excl0']}).append(pd.DataFrame({'dataset':'LeukKnown', 'scores':neighbor_scores_knowns['pct_neighbors_paired_excl0']}).append(pd.DataFrame({'dataset':'LeukNew', 'scores':neighbor_scores['pct_neighbors_paired_excl0']})))
neighbor_scores_combined.reset_index(inplace = True)
print(pairwise_tukeyhsd(neighbor_scores_combined['scores'], neighbor_scores_combined['dataset']))
