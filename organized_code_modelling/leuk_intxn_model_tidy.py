from os import chdir as os_chdir
from scipy import constants
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches  # for custom legends

os_chdir("./organized_code_modelling")


if_vienna_avgtypes = True
if_pharm_avgtypes = False

##### celltype parameters #####
if if_vienna_avgtypes:   # cell types in the already-published Vienna study differ from our own generated dataset
    df_cell_info = pd.read_csv("cell_freq_table_vienna.csv", index_col = 'Celltype')
elif if_pharm_avgtypes:
    df_cell_info = pd.read_csv("cell_freq_table_zurich.csv", index_col = 'Celltype')
else:
    print("Warning: have not picked cell type set")
#end


v_cell_names = df_cell_info.index.to_series()
v_cell_pcts = df_cell_info['Percent_PBMC']
v_cell_rad = df_cell_info['Diameter_um'] / 2

def calculate_SA(cell_radius):
    ''' input cell radius and outputs surface area in the same units squared'''
    return 4 * constants.pi * (cell_radius)**2
#end

cellTotal_conc = 2.5E6                         # cells / mL total cell concentration for whole PBMCs  (ends up not mattering in this model formulation)
v_cell_concs = cellTotal_conc * v_cell_pcts    # cells / mL (not used elsewhere in this version of the model)
v_cell_SA = calculate_SA(v_cell_rad)           # um^2

m_proportions = v_cell_pcts.to_numpy()[:, None] * v_cell_pcts.to_numpy()
m_proportions = pd.DataFrame(m_proportions, columns = v_cell_names.rename(None), index = v_cell_names.rename(None))


##### molecule copy numbers ##### 
prot_condition = "steady.state"
proteomics_csv = pd.read_csv("Rieckmann_proteomics_summary_mean_leuk.csv")
proteomics_csv.set_index('Uniprot', inplace = True)

def molecule_matrices(proteomics_csv,
                      condition = "steady.state",     # "steady.state"   or    "activated"
                      pick_mostabundant = False,      # match proteomics data cell type labales to microscopy labels by picking most abundant single subtype within the match category
                      pick_vienna_avgtypes = True,    # match proteomics data cell type labales to microscopy labels by averaging the subtypes
                      pick_pharm_avgtypes = False,    # match based on averaging for when using our pharmacoscopy cell types
                      pick_null_celltypes = False,        # null model where randomize cell type labels on expression matrices, but keep everything else constant  
                      inp_null_celltypes_vec = np.array([]),      # if pick_null_celltypes = True, optionally can manually specify the randomized vector here, which is useful given the small number of possible permutations 
                      cd = {'pT4':0.70, 'pBn':0.66, 'pBm':0.34, 'pMc':0.84, 'pMn':0.09, 'pMi':0.07, 
                            'pDm':0.65, 'pDp':0.35, 'pT4n':0.41, 'pT4c':0.39, 'pT4e':0.19, 'pT4r':0.01,
                            'pT8n':0.39, 'pT8c':0.07, 'pT8e':0.22, 'pT8r':0.32, 'pNd':0.90, 'pNb':0.10,
                            'pTh1':0.17, 'pTh2':0.01, 'pTh17':0.01, 'pTreg':0.04,
                            'fTh1em':0.45 , 'fTh1cm':0.55, 'fTh2em':0.45, 'fTh2cm':0.55, 'fTregem':0.37, 'fTregcm':0.63, 'fTregmem': 0.80}): 
    '''
    Converts matrix of celltype expression in absolute protein counts (such as Rieckmann data)
    to 2D concentrations formatted for vienna celltypes. \n
    Involves a translation between proteomics data's labels and vienna labels, for which options given.  \n
    Returns matrix of molecule 2D concentrations and vector of molecule rownames
    '''
    
    if pick_mostabundant and pick_vienna_avgtypes:
        m_molecules = proteomics_csv.filter(regex = (("_"+condition+"|").join(df_cell_info.Rieckmann_mostabundant))+"_"+condition).copy()
        renamingcells = dict(zip(df_cell_info.Rieckmann_mostabundant + "_" + condition, df_cell_info.index))
        m_molecules.rename(columns = renamingcells, inplace = True)
        m_molecules = m_molecules.reindex(columns = renamingcells.values(), copy = True)
    elif pick_vienna_avgtypes:
        m_molecules = pd.DataFrame()    
        m_molecules['Bcell'] = cd['pBn']*proteomics_csv['B.naive'+"_"+condition] + cd['pBm']*proteomics_csv['B.memory'+"_"+condition] 
        m_molecules['Mono'] = cd['pMc']*proteomics_csv['MO.classical'+"_"+condition] + cd['pMn']*proteomics_csv['MO.nonclassical_steady.state'] + cd['pMi']*proteomics_csv['MO.intermediate_steady.state']
        m_molecules['DC'] = cd['pDm']*proteomics_csv['mDC'+"_"+condition] + cd['pDp']*proteomics_csv['pDC'+"_"+condition] 
        m_molecules['Tcell'] = (cd['pT4']*(cd['pT4n']*proteomics_csv['T4.naive'+"_"+condition] + cd['pTh1']*proteomics_csv['Th1_steady.state'] + cd['pTh2']*proteomics_csv['Th2_steady.state'] + cd['pTreg']*cd['fTregmem']*proteomics_csv["mTregs"+"_"+condition] + cd['pTreg']*(1-cd['fTregmem'])*proteomics_csv["nTregs"+"_"+condition] + 
                                 (cd['pT4c'] - (cd['pTh1']*cd['fTh1cm']) - (cd['pTh2']*cd['fTh2cm']) - (cd['pTreg']*cd['fTregcm']))*proteomics_csv['T4.CM'+"_"+condition] + (cd['pT4e'] - (cd['pTh1']*cd['fTh1em']) - (cd['pTh2']*cd['fTh2em']) - (cd['pTreg']*cd['fTregem']))*proteomics_csv['T4.EM'+"_"+condition] + cd['pT4r']*proteomics_csv['T4.EMRA'+"_"+condition]) +
                                (1-cd['pT4'])*(cd['pT8n']*proteomics_csv['T8.naive'+"_"+condition] + cd['pT8c']*proteomics_csv['T8.CM'+"_"+condition] + cd['pT8e']*proteomics_csv['T8.EM'+"_"+condition] + cd['pT8r']*proteomics_csv['T8.EMRA'+"_"+condition]))
        m_molecules['NKcell'] = cd['pNd']*proteomics_csv['NK.dim'+"_"+condition] + cd['pNb']*proteomics_csv['NK.bright'+"_"+condition] 
    elif pick_pharm_avgtypes:
        m_molecules = pd.DataFrame()    
        m_molecules['B'] = cd['pBn']*proteomics_csv['B.naive'+"_"+condition] + cd['pBm']*proteomics_csv['B.memory'+"_"+condition] 
        m_molecules['M14'] = cd['pMc']/(1-cd['pMn'])*proteomics_csv['MO.classical'+"_"+condition] + cd['pMi']/(1-cd['pMn'])*proteomics_csv['MO.intermediate_steady.state']
        m_molecules['M16'] = proteomics_csv['MO.nonclassical_steady.state']
        m_molecules['DC'] = cd['pDm']*proteomics_csv['mDC'+"_"+condition] + cd['pDp']*proteomics_csv['pDC'+"_"+condition] 
        m_molecules['T4'] = (cd['pT4n']*proteomics_csv['T4.naive'+"_"+condition] + cd['pTh1']*proteomics_csv['Th1_steady.state'] + cd['pTh2']*proteomics_csv['Th2_steady.state'] + cd['pTreg']*cd['fTregmem']*proteomics_csv["mTregs"+"_"+condition] + cd['pTreg']*(1-cd['fTregmem'])*proteomics_csv["nTregs"+"_"+condition] + 
                                 (cd['pT4c'] - (cd['pTh1']*cd['fTh1cm']) - (cd['pTh2']*cd['fTh2cm']) - (cd['pTreg']*cd['fTregcm']))*proteomics_csv['T4.CM'+"_"+condition] + (cd['pT4e'] - (cd['pTh1']*cd['fTh1em']) - (cd['pTh2']*cd['fTh2em']) - (cd['pTreg']*cd['fTregem']))*proteomics_csv['T4.EM'+"_"+condition] + cd['pT4r']*proteomics_csv['T4.EMRA'+"_"+condition])
        m_molecules['T8'] = cd['pT8n']*proteomics_csv['T8.naive'+"_"+condition] + cd['pT8c']*proteomics_csv['T8.CM'+"_"+condition] + cd['pT8e']*proteomics_csv['T8.EM'+"_"+condition] + cd['pT8r']*proteomics_csv['T8.EMRA'+"_"+condition]
        m_molecules['NK'] = cd['pNd']*proteomics_csv['NK.dim'+"_"+condition] + cd['pNb']*proteomics_csv['NK.bright'+"_"+condition] 
    #end
    v_molec_names = m_molecules.index.to_series()
    
    m_molec_conc = m_molecules.divide(v_cell_SA, axis = 1)
    
    if pick_null_celltypes:
        
        if inp_null_celltypes_vec.any():
            shuffle_vec_indxs = inp_null_celltypes_vec
            
        else:
            shuffle_vec_indxs = np.random.permutation(np.arange(m_molec_conc.shape[1]))
        #end
        m_molec_conc_original = np.array(m_molec_conc)
        m_molec_conc = pd.DataFrame(m_molec_conc_original[:, shuffle_vec_indxs],
                                    index = m_molec_conc.index, columns = m_molec_conc.columns)
    #end
    
    return m_molec_conc, v_molec_names
#end

cellavg_proportions = {'pT4':0.70, 'pBn':0.66, 'pBm':0.34, 'pMc':0.84, 'pMn':0.09, 'pMi':0.07, 
                       'pDm':0.65, 'pDp':0.35, 'pT4n':0.41, 'pT4c':0.39, 'pT4e':0.19, 'pT4r':0.01,
                       'pT8n':0.39, 'pT8c':0.07, 'pT8e':0.22, 'pT8r':0.32, 'pNd':0.90, 'pNb':0.10,
                       'pTh1':0.17, 'pTh2':0.01, 'pTh17':0.01, 'pTreg':0.04, 
                       'fTh1em':0.45 , 'fTh1cm':0.55, 'fTh2em':0.45, 'fTh2cm':0.55, 'fTregem':0.37, 'fTregcm':0.63, 'fTregmem': 0.80}

m_molec_conc, v_molec_names = molecule_matrices(proteomics_csv,
                                                condition = prot_condition,
                                                pick_mostabundant = False, 
                                                pick_vienna_avgtypes = if_vienna_avgtypes, pick_pharm_avgtypes = if_pharm_avgtypes,
                                                cd = cellavg_proportions)

if prot_condition == "steady.state":  # non-activated cells do not put HLA-F on their surface (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3867582/)
    m_molec_conc.loc[m_molec_conc.index == "P30511", :] = np.zeros(m_molec_conc.iloc[0, :].shape)
#end

##### interaction matrix ####
interactions_csv = pd.read_csv("Leuk_interaction_affinities.csv",
                               usecols =[i for i in range(8)])

interactions_csv.Gene_L.replace(r'-[^-]+$', '', regex = True, inplace = True)
interactions_csv.Gene_R.replace(r'-[^-]+$', '', regex = True, inplace = True)

def pair_labeler(str1, str2, sep = " + "):
    ''' Combines two strings into one in alphabetical order (so unique). Returns str'''
    return np.where(str1 < str2, (str1 + sep + str2), (str2 + sep + str1))
#end
interactions_csv['gene_pair'] = interactions_csv.apply(lambda row: pair_labeler(row.Gene_L, row.Gene_R), axis = 1).astype(str)
interactions_csv['pair_label'] = interactions_csv.apply(lambda row: pair_labeler(row.Uniprot_L, row.Uniprot_R), axis = 1).astype(str)

def interaction_matrices(interactions_csv,
                         v_molec_names,
                         pick_excludeintxns = [],         # any interactions to exclude (can pass [] as blank)
                         pick_nullintxn = False,          # randomly permute interaction network to generate a null distribution
                         rand_seed = None):               # if randomly permuting, can optionally provide seed for RNG, otherwise defaults to python's standard rand seed    
    '''
    Reformats interaction table in list format from csv into matrices of all protein-protein pairs and their affinities  \n
    Includes options to create "null models" where interactions or affinities are scrambled, to see how much the specific interaction network 
    contributes to the results (and how much may be explainable by general factors like a celltype's average expression magnitude). \n
    Returns an updated version of the input csv table and a binary intxn matrix and a Kd matrix
    '''
    if len(pick_excludeintxns) > 0:
        interactions_csv = interactions_csv[~interactions_csv['pair_label'].isin(pick_excludeintxns)]  
    #end
    
    m_interaction_names = pd.DataFrame((v_molec_names.values[:,None] + " + " + v_molec_names.values),
                                  columns = v_molec_names, index = v_molec_names)
    
    if pick_nullintxn:
        uniprots_with_intxn = [*interactions_csv['Uniprot_L'].to_list(), *interactions_csv['Uniprot_R'].to_list()]
        uniprots_with_intxn = pd.Series(list(set(uniprots_with_intxn)))
        uniprots_intxn_names_pairs = pd.DataFrame((uniprots_with_intxn.values[:,None] + " + " + uniprots_with_intxn.values),
                              columns = uniprots_with_intxn, index = uniprots_with_intxn) 
        pairlabels_rand = uniprots_intxn_names_pairs.melt().sample(n = interactions_csv.shape[0], random_state = rand_seed)
        interactions_csv.loc[:, 'pair_label'] = pairlabels_rand['value'].values
    #end
    
    remappingkd = dict(zip(interactions_csv[interactions_csv.Kd_nM.notnull()].pair_label, interactions_csv[interactions_csv.Kd_nM.notnull()].Kd_nM))
    remappingkd_recip = dict(zip([" + ".join(intxn_name_i.split(" + ")[::-1]) for intxn_name_i in list(remappingkd.keys())],  # reverses ordering of protein names 
                                 remappingkd.values()))
    remappingkd.update(remappingkd_recip)
    m_interactions_Kd3D = m_interaction_names.applymap(lambda x: remappingkd.get(x, 0))
    
    m_interactions = m_interactions_Kd3D > 0
    
    return interactions_csv, m_interactions, m_interactions_Kd3D
#end

(interactions_csv, 
 m_interactions, m_interactions_Kd3D) = interaction_matrices(interactions_csv,
                                                            v_molec_names,
                                                            pick_excludeintxns = ['C10orf54 + HLA-E', 'C10orf54 + HLA-F', 'TNF + TNR1A'],  # non-surface
                                                            pick_nullintxn = False)

def calculate_3Dto2D(Kd_3D, height = 0.010):
    ''' Converts 3D Kd in Molar units to 2D Kd in units moelcules/um2 (height is in um units)  \n
        uses constant height factor to get rough estimates since very few are empirically measured'''
    return constants.Avogadro * 1E-15 * height * Kd_3D
#end

m_interactions_Kd2D = calculate_3Dto2D(m_interactions_Kd3D * 1E-9, 
                                       height = 0.010)

m_interactions_Kd2D.loc["P06729", "P19256"] = m_interactions_Kd2D.loc["P19256", "P06729"] = 1.1 # rare case of litterature-published 2D Kd can use instead of inferred value (doi:10.1007/s10439-005-2504-5)

def bound_concentration(conc_A_total, conc_B_total, Kd_2D):
    '''
    Calculates amount of bound A-B given the total amount of A and B and their 2D affinity \n
     dervied from the quadratic solution the system of equations \n
     Kd = Rec * Lig / Bound \n
     Rec = Rec_tot - Bound  \n
     Lig = Lig_tot - Bound  \n
    Where all concentrations and Kd are in the same units (e.g. molecules/um2)
    Beware this equation does not handle case if 0 is input as a concentration or Kd (because smaller numbers = stronger binding so zero is perfect conversion to bound state)
    '''
    conc_AB = ((Kd_2D + conc_A_total + conc_B_total) -
               ((Kd_2D + conc_A_total + conc_B_total)**2 -
               (4 * conc_A_total * conc_B_total))**(1/2)) / 2
    return conc_AB
#end

# extract from binary intxn matrix pairs of molecule indexes
v_interaction_tuples = list(zip(*np.asarray(m_interactions_Kd2D).nonzero()))

def bound_matrices(v_interaction_tuples, v_molec_names, 
                   m_molec_conc, m_interactions_Kd2D, v_cell_names):
    '''
    Iterates through all the binary interactions with Kd values, calculates bound concentrations for each.
     Sums up total bound concentration for each cell-cell pair across all the protein interactions. \n
    Returns matrix of cell pairs with sum of bound concentrations, and list of bound concentrations per protein interaction that fed into the sum    
    '''    
    t_bound = []
    sum_bound = pd.DataFrame(np.zeros((len(v_cell_names),len(v_cell_names))), columns = v_cell_names, index = v_cell_names)
    bound_rownames = pd.Series(["null"] * len(v_interaction_tuples), dtype = "str")

    for n, intxn_n in enumerate(v_interaction_tuples):
        intxn_name = "%s + %s"%(v_molec_names[intxn_n[0]], v_molec_names[intxn_n[1]])
        bound_rownames[n] = intxn_name
        
        conc_A_n = m_molec_conc.iloc[intxn_n[0], :]
        conc_B_n = m_molec_conc.iloc[intxn_n[1], :]
        
        conc_A = conc_A_n.to_frame()
        conc_B = conc_B_n.to_frame().transpose()
        
        m_bound = bound_concentration(conc_A.to_numpy()[:, None],
                                       conc_B.to_numpy(),  
                                       Kd_2D = m_interactions_Kd2D.iloc[intxn_n])

        m_bound = np.minimum(m_bound, np.minimum(conc_A.to_numpy()[:, None], conc_B.to_numpy()))
        m_bound = m_bound.squeeze(axis = 1) 
        
        sum_bound += pd.DataFrame(m_bound, columns = v_cell_names, index = v_cell_names)
        t_bound.append(m_bound)
    #end
    return sum_bound, t_bound, bound_rownames
#end

(sum_bound, t_bound, bound_rownames) = bound_matrices(v_interaction_tuples, v_molec_names, 
                                                      m_molec_conc, m_interactions_Kd2D, v_cell_names)


##### fitting kinetic model's matrix of bound sums to vienna data #####

from sklearn import linear_model, metrics
from scipy import stats as stats
import statsmodels.api as sm

m_contact_obs =  pd.DataFrame({'DC':    pd.Series([2.4,     1.4,     np.nan,    0.2]),  # data : https://pubmed.ncbi.nlm.nih.gov/28437395/
                               'Mono':  pd.Series([1.3,     1.5,     0.4,       np.nan]),
                               'Bcell': pd.Series([np.nan,  0.3,     0.2,       0.0]),
                               'Tcell': pd.Series([0.0,     np.nan,  0.0,       0.0])})
m_contact_obs.set_index(m_contact_obs.columns, inplace = True)
m_contact_obs = m_contact_obs[v_cell_names[0:4]]  # re-order columns to match other matrices
m_contact_obs = m_contact_obs.reindex(v_cell_names[0:4], axis = 'index')


x_obs = sum_bound.drop('NKcell', axis=1).drop('NKcell', axis=0).to_numpy().flatten()
y_obs = m_contact_obs.to_numpy().flatten()

linreg = linear_model.LinearRegression()  # import model (ordinary least squares)
mask = np.isfinite(x_obs) & np.isfinite(y_obs)
linreg.fit(X = x_obs[mask, None],   
           y = y_obs[mask, None])

beta = linreg.coef_
intercept = linreg.intercept_
fit_r2 = metrics.r2_score(y_pred = linreg.predict(x_obs[mask, None]),
                          y_true = y_obs[mask, None])  # calculate a few different ways to check correct
fit_r2 = pd.DataFrame({'sumbound':x_obs[mask], 'obscontact':y_obs[mask]}).corr(method = "pearson") ** 2
fit_r2 = linreg.score(X = x_obs[mask, None], 
                      y = y_obs[mask, None])

fit_stats = sm.OLS(y_obs[mask], sm.add_constant(x_obs[mask], prepend = False))
fit_stats = fit_stats.fit()
fit_pval = fit_stats.pvalues[0]  # fit_stats.summary()

colors_dict = {'Bcell':'#de412c', 'Mono':'#ba91bd', 'DC':'#618f59', 'Tcell':'#619fcf'}
color_m_contact_obs = m_contact_obs.rename(columns = colors_dict).rename(index = colors_dict)
m_labels_flat = pd.Series((color_m_contact_obs.columns.values[:, None] + "," + color_m_contact_obs.index.values).flatten())
m_labels_flat.replace(colors_dict, inplace = True)
m_labels = pd.DataFrame({'Cell1': m_labels_flat.replace(r',.*$', '', regex = True),
                         'Cell2': m_labels_flat.replace(r'^.*,', '', regex = True)})

xval = x_obs[mask]   # manually define 95% confidence interval
yval = y_obs[mask]
yfit = x_obs[mask]*fit_stats.params[0] + fit_stats.params[1]
dof = xval .size - 2
t = stats.t.ppf(0.975, df = dof)
resid = yval - yfit 
chi2 = np.sum((resid / yfit)**2)
chi2_red = chi2 / dof
s_err = np.sqrt(np.sum(resid**2) / dof)
x2 = np.linspace(xval.min(),xval.max(),500)
y2 = np.linspace(xval.min(),xval.max(),500)*fit_stats.params[0] + fit_stats.params[1]
ci = t * s_err * np.sqrt(1/xval.size + (x2 - np.mean(xval))**2 / np.sum((xval - np.mean(xval))**2))

fig, ax = plt.subplots(dpi = 300)  # dpi = 150 for console printing  ; else 300 to save
ax.fill_between(x2, y2 + ci, y2 - ci, facecolor ="#918e86", edgecolor = "none", alpha = 0.12)
plt.plot(np.linspace(x_obs[mask].min(),x_obs[mask].max(),500), np.linspace(x_obs[mask].min(),x_obs[mask].max(),500)*fit_stats.params[0] + fit_stats.params[1],
         linestyle = ":", color = "grey")
ax.scatter(x = x_obs[mask],
           y = y_obs[mask],
           s = 50,           # size of points
           c = m_labels.Cell1[mask])
ax.scatter(x = x_obs[mask],
           y = y_obs[mask],
           s = 10,           # size of points
           c = m_labels.Cell2[mask])
fig.text(0.16, 0.56, r'$r^2 = %.3f$' '\n' r'$p = %.3f$'%(fit_r2, fit_pval), fontsize=10)
ax.set_xlabel('Kinetic model interaction score')
ax.set_ylabel('Vienna microscopy contact score')
ax.legend(handles = [mpatches.Patch(color=colors_dict[key], label=key) for key in colors_dict],
          loc = 'upper left')
#fig.savefig('kineticmodel_v8_viennafig_cellavg_linreg.png', bbox_inches='tight')


#### export values to csv ####

names_cellpair = (sum_bound.columns.values[:, None] + " + " + sum_bound.index.values).flatten()  # make square matrix of cell-cell pairs in score table, then flatten to vector as above for 'y_obs' etc
names_intxnpair = [v_molec_names[tuple_i[0]] + " + " + v_molec_names[tuple_i[1]] if v_molec_names[tuple_i[0]] < v_molec_names[tuple_i[1]] else v_molec_names[tuple_i[1]] + " + " + v_molec_names[tuple_i[0]] for tuple_i in v_interaction_tuples]

table_bound_raw = pd.DataFrame([t_sub.flatten() for t_sub in t_bound])  # convert each adhesion matrix into vector
table_bound_raw.columns = names_cellpair   # insert cell pair (col) and protein intxn (row) names
table_bound_raw.index = names_intxnpair
table_bound_raw.index.name = "Interaction"

table_bound_raw['Pair_label'] = table_bound_raw.index.map(dict(zip(interactions_csv.pair_label, interactions_csv.gene_pair)))
#table_bound_raw.to_csv("kineticmodel_v8_zurich_resting_leukmanual_tablebound.csv")
