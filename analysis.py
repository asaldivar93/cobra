import os
import pandas as pd
import cobra
import warnings
import plotly.express as px
import seaborn as sns
import matplotlib.pyplot as plot

from numpy import sign
from numpy import log10
from cobra.flux_analysis import micom, plot_micom_tradeoff
from tools.cobra_tools import add_boundaries_for_simulation
from tools.p_tools import build_genus_database
from tools.plotting import plot_heatmap

warnings.filterwarnings("ignore")
os.chdir('/home/alexis/UAM/cobra/')

# %% codecell
taxonomy = pd.read_csv('data_files/silva_taxonomy.csv', sep ='\t')
taxonomy.set_index('sequence_identifier', inplace = True)
taxonomy.fillna('unknown', inplace = True)

seq_counts = pd.read_csv('data_files/otu-table-corrected.tsv', sep = '\t', index_col = 'OTU_ID')
seq_counts.loc['0281ce539a44c8b1b059ae85c15e3fab', 'My_20'] = 0  # Outlier Mesoaciditoga

samples = ['My_20', 'CIR_19']
taxa_database = build_genus_database(samples, seq_counts, taxonomy)
taxa_database

# %% codecell
community_model = cobra.Model()
community_model.solver = 'cplex'
sample = taxa_database.loc[:, 'sample'] == 'CIR_19'
sample_db = taxa_database[sample].set_index('id')
for k in taxa_database.loc[sample].index:
    new_model = cobra.io.read_sbml_model(taxa_database.loc[k, 'model_path'])
    community_model.merge(new_model)

to_close = ['Ubiquinols', 'CO_A', 'FMN', '10_FORMYL_THF', 'FAD', 'GLUTATHIONE',
            '2_OCTAPRENYLPHENOL', 'METHYLENE_THF', 'PYRIDOXAL_PHOSPHATE',
            'S_ADENOSYLMETHIONINE', 'THF', 'Menaquinols', 'PROTOHEME']
for ex in to_close:
    for rxn in community_model.reactions:
        if 'EX' in rxn.id and ex in rxn.id:
            rxn.lower_bound = 0
            rxn.upper_bound = 0

for taxa_k in sample_db.index:
    for rxn in community_model.reactions:
        if 'EX' in rxn.id and taxa_k in rxn.id:
            proton = ['PROTON' in met.id for met in rxn.metabolites]
            if any(proton):
                rxn.add_metabolites({'_PROTON_pe_' + taxa_k: -1, '_PROTON_cy_' + taxa_k: 1})

rxn = community_model.reactions._ALCOHOL_DEHYDROG_RXN_Meis
rxn.lower_bound = 0
rxn.upper_bound = 100

# %% codecell
media = ['_CH4_ou', '_OXYGEN_MOLECULE_ou', '_NITRATE_ou', '_FE_2_ou',
         '_Pi_ou', '_SULFATE_ou', '_NA__ou', '_MG_2_ou', '_CO_2_ou',
         '_CL__ou']
sample_db = taxa_database[sample].set_index('id')
with community_model as model:
    add_boundaries_for_simulation(model, media=media, carbon_uptake=0.00625)
    nonzero = plot_micom_tradeoff(model, sample_db)

# %% codecell
media = ['_CH4_ou', '_OXYGEN_MOLECULE_ou', '_NITRATE_ou', '_FE_2_ou',
         '_Pi_ou', '_SULFATE_ou', '_NA__ou', '_MG_2_ou', '_CO_2_ou',
         '_CL__ou']
sample_db = taxa_database[sample].set_index('id')

with community_model as model:
    add_boundaries_for_simulation(model, media=media, carbon_uptake=0.00625)
    micom_model, linear_solution, quadratic_solution, micom_solution, taxa_growth_rates = micom(
        model, 0.2, sample_db
    )
print(micom_model.summary(solution = micom_solution))
print(micom_model.metabolites._METOH_pe_Meis.summary(solution = micom_solution))
micom_solution.status

# %% codecell
ex_mets = []
for rxn in community_model.reactions:
    if 'EX' in rxn.id:
        ex_mets.append(rxn.id[4:-5])
external_metabolites = list(pd.Series(ex_mets).unique()[0:2])
b = list(pd.Series(ex_mets).unique()[13:])
external_metabolites.extend(b)
external_metabolites = set(external_metabolites)
to_close = set(to_close)
external_metabolites = external_metabolites - to_close

# %%codecell
import_fluxes = pd.DataFrame(index = external_metabolites, columns = sample_db.index)
export_fluxes = pd.DataFrame(index = external_metabolites, columns = sample_db.index)
external_fluxes = pd.DataFrame(index = external_metabolites, columns = sample_db.index)
for rxn in micom_model.reactions:
    if '_EX' in rxn.id:
        metabolite = rxn.id[4:-5]
        genus_k = rxn.id[-4:]
        abundance_k = sample_db.loc[genus_k, 'relative_abundance']
        if metabolite in import_fluxes.index:
            if micom_solution[rxn.id] > 0:
                import_fluxes.loc[metabolite, genus_k] = 1000 * micom_solution[rxn.id] * abundance_k * list(rxn.metabolites)[0].cmol_by_mol
            elif micom_solution[rxn.id] < 0:
                export_fluxes.loc[metabolite, genus_k] = -1000 * micom_solution[rxn.id] * abundance_k * list(rxn.metabolites)[0].cmol_by_mol

            external_fluxes.loc[metabolite, genus_k] = 1000 * micom_solution[rxn.id] * abundance_k * list(rxn.metabolites)[0].cmol_by_mol

# %%codecell
min_I = 7.787481e-12
min_II = 7.787481e-12
external_fluxes.fillna(0, inplace = True)
import_fluxes.fillna(0, inplace=True)
import_fluxes = import_fluxes + 2 * min_I
import_fluxes = log10(import_fluxes)

export_fluxes.fillna(0, inplace=True)
export_fluxes = export_fluxes + 2 * min_II
export_fluxes = log10(export_fluxes)
exp = import_fluxes - export_fluxes
import_fluxes.reset_index(inplace=True)
export_fluxes.reset_index(inplace=True)

# %%codecell
plot_heatmap(exp, cmap = sns.diverging_palette(263, 244, s = 79, l = 65, sep = 20, as_cmap = True))

external_fluxes.drop(index='CH4', inplace=True)
# %% codecell
interaction_matrix = pd.DataFrame(index = external_fluxes.columns, columns = external_fluxes.columns, dtype = 'float')
for k in interaction_matrix.index:
    imports = external_fluxes[k] > 0
    total_imports_k = external_fluxes.loc[imports, k].sum()
    for j in interaction_matrix.columns:
        flux_from_j = -external_fluxes.loc[imports, j].sum()
        total_external = external_fluxes.loc[imports].drop(columns = [j, k], axis=1).sum()
        imports_others = -total_external[total_external > 0].sum()
        exports_others = -total_external[total_external < 0].sum()
        if sign(flux_from_j) == sign(total_imports_k):
            IC_j_in_k = abs(total_imports_k - flux_from_j) / (flux_from_j + exports_others)
        else:
            if abs(flux_from_j) >= total_imports_k:
                IC_j_in_k = (total_imports_k + flux_from_j) / exports_others
            else:
                IC_j_in_k = -1 + (total_imports_k + flux_from_j) / total_imports_k
        interaction_matrix.loc[k, j] = IC_j_in_k

sns.heatmap(interaction_matrix, cmap = sns.diverging_palette(26, 145, s=75, l=67, sep=1, center='light', as_cmap=True))

# %% codecell
interactions_long = pd.DataFrame(columns=['reciver', 'giver', 'ic'])
for k in interaction_matrix.index:
    for j in interaction_matrix.columns:
        if not k == j:
            interactions_long = interactions_long.append(
                pd.DataFrame(
                    [[k, j, interaction_matrix.loc[k, j]]],
                    columns=['reciver',
                             'giver',
                             'ic'
                             ]
                ),
            )

# %% codecell
import numpy as np
interactions_long['log_ic'] = interactions_long['ic'] + 2
interactions_long['log_ic'] = np.log10(interactions_long['log_ic'])
px.histogram(interactions_long, x='log_ic', template='none')
sns.ecdfplot(interactions_long,x='log_ic')

px.box(
    interactions_long,
    y='ic',
    x='giver',
    color='giver',
    template='none',
    )

px.box(
    interactions_long,
    y='ic',
    x='reciver',
    color='reciver',
    template='none'
)
