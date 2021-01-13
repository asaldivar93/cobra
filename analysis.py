import pandas as pd
import arviz

from stan.run_stan import run_sampling
from cobra.flux_analysis import community
import plotly.express as px

import random
from tqdm import trange
from numpy import sign
from tools.cobra_tools import add_boundaries_for_simulation
from cobra.flux_analysis import micom

# %% codecell
home_path = '/home/alexis/UAM/cobra/'
taxonomy_path = 'data_files/silva_taxonomy.csv'
counts_path = 'data_files/otu-table-corrected.tsv'
results_path = 'results_all_open/'
com = community.community(
    home_path=home_path, taxonomy_path=taxonomy_path, counts_path=counts_path, results_path=results_path
    )
media = ['_CH4_ou', '_OXYGEN_MOLECULE_ou', '_NITRATE_ou', '_FE_2_ou',
         '_Pi_ou', '_SULFATE_ou', '_NA__ou', '_MG_2_ou', '_CO_2_ou',
         '_CL__ou']
yield_mets = dict(EX__OXYGEN_MOLECULE_ou = 'O2', _DM_CARBON_DIOXIDE_ex = 'CO2')
# %% codecell
samples = dict(CIR_19 = 0.0063, My_20 = 0.00913)
com.set_samples(samples)
com.load_models(samples, close_ex=False)

# %% codecell
fig_tradeoff = com.plot_micom_tradeoff(media)

# %% codecell

fig, y_df = com.plot_yields(
    media=media,
    carbon_uptake=0.00625,
    metabolites=yield_mets,
    trade_off=0.9
)
fig.show()

fig = com.plot_metoh_ex(
    media=media,
    carbon_uptake=0.00625
)
fig.show()

com.summary_ex_fluxes(
    media=media,
    carbon_uptake=0.00625,
    trade_off=0.9
)
fig = com.plot_interaction_matrix((9, 4))

# %%codecell
ic_cir = pd.read_csv(results_path + 'interaction_matrix/ic_CIR_19.csv', index_col=0)
ic_my = pd.read_csv(results_path + '/interaction_matrix/ic_My_20.csv', index_col=0)
fit_cir_pos = run_sampling(ic_cir, chains = 6, iter = 2000, positives = True)
fit_cir_neg = run_sampling(ic_cir, chains = 6, iter = 2000, positives = False)
fit_my_pos = run_sampling(ic_my, chains = 6, iter = 2000, positives = True)
fit_my_neg = run_sampling(ic_my, chains = 6, iter = 2000, positives = False)
print(fit_cir_pos)
print(fit_cir_neg)
print(fit_my_pos)
print(fit_my_neg)
arviz.plot_trace(fit_cir_neg, var_names = ['mu', 'phi'])
arviz.plot_posterior(fit_cir_neg, var_names = ['mu', 'phi'])


# %%codecell
exrxns_meis = []
for rxn in com.models[com.samples[0]].reactions:
    if 'EX' in rxn.id and 'Meis' in rxn.id:
        exrxns_meis.append(rxn.id)

exrxns_meis = exrxns_meis[13:]
ex_mets = []
for rxn in com.models[com.samples[0]].reactions:
    if 'EX' in rxn.id:
        ex_mets.append(rxn.id[4:-5])
external_metabolites = list(pd.Series(ex_mets).unique()[0:2])
b = list(pd.Series(ex_mets).unique()[13:])
external_metabolites.extend(b)
external_metabolites = set(external_metabolites)
to_close = set(com.permanent_close)
external_metabolites = external_metabolites - to_close
n_ex = len(exrxns_meis)
columns = ['Giver', 'Reciver', 'Interaction Coefficient', 'Sample', 'N Blocked']
# %%codecell
# columns = ['Giver', 'Reciver', 'Interaction Coefficient', 'Sample', 'N Blocked']
# with open(com.results_path + 'simulations.csv', 'w+') as file:
#     pd.DataFrame(columns=columns).to_csv(file, index=False)
# with open(com.results_path + 'metoh_ex.csv', 'w+') as file:
#     pd.DataFrame(columns=['Sample', 'Exchanged Metanol', 'N Blocked']).to_csv(file, index=False)
# with open(com.results_path + 'yields.csv', 'w+') as file:
#     pd.DataFrame(columns=['Sample', 'Metabolite', 'Yield', 'N Blocked']).to_csv(file, index=False)
# %%codecell
step = 6
n_sim = 24
carbon_uptake = 0.00625
missing_dels = [a for a in range(1, int(n_ex / 2), step)]
missing_dels = [7]
sample = 'CIR_19'
sample_filter = com.taxa_database.loc[:, 'sample'] == sample
sample_db = com.taxa_database[sample_filter].set_index('id')
for max_del in missing_dels:
    for n in trange(n_sim):
        index = [random.randint(0, n_ex - 1) for n in range(max_del)]
        with com.models[sample] as model:
            add_boundaries_for_simulation(
                model, media=media, carbon_uptake=0.00625
                )
            for met in pd.Series(exrxns_meis).iloc[index]:
                rxn = model.reactions.get_by_id(met)
                rxn.lower_bound = -100
                rxn.upper_bound = 0
            micom_model, linear_solution, quadratic_solution, micom_solution, taxa_growth_rates = micom(
                model, 0.9, sample_db
            )

        external_fluxes = pd.DataFrame(index = external_metabolites, columns = sample_db.index)
        for rxn in micom_model.reactions:
            if '_EX' in rxn.id:
                metabolite = rxn.id[4:-5]
                genus_k = rxn.id[-4:]
                abundance_k = sample_db.loc[genus_k, 'relative_abundance']
                if metabolite in external_fluxes.index:
                    external_fluxes.loc[metabolite, genus_k] = 1000 * micom_solution[rxn.id] * abundance_k * list(rxn.metabolites)[0].cmol_by_mol

        external_fluxes.drop(index='CH4', inplace=True)
        external_fluxes.fillna(0, inplace = True)

        interaction_matrix = pd.DataFrame(index = external_fluxes.columns, columns = external_fluxes.columns, dtype = 'float')
        for k in interaction_matrix.index:
            imports = external_fluxes[k] > 0
            total_imports_k = external_fluxes.loc[imports, k].sum()
            for j in interaction_matrix.columns:
                flux_from_j = -external_fluxes.loc[imports, j].sum()
                total_external = external_fluxes.loc[imports].drop(columns = [j, k], axis=1).sum()
                imports_others = -total_external[total_external > 0].sum()
                exports_others = -total_external[total_external < 0].sum()
                if k == j:
                    IC_j_in_k = 0
                elif sign(flux_from_j) == sign(total_imports_k):
                    IC_j_in_k = (total_imports_k - (total_imports_k - flux_from_j)) / (flux_from_j + exports_others)
                else:
                    # IC_j_in_k = -(exports_others - (exports_others + flux_from_j)) / (exports_others + imports_others)
                    IC_j_in_k = flux_from_j / (exports_others + imports_others)
                interaction_matrix.loc[k, j] = IC_j_in_k
        ic_long = com.ic_matrix_melt(interaction_matrix)
        ic_long.loc[:, 'Sample'] = sample
        ic_long.loc[:, 'N Blocked'] = max_del
        with open(com.results_path + 'simulations.csv', 'a+') as simulations:
            ic_long.to_csv(simulations, index=False, header=False)

        summary = micom_model.metabolites._METOH_pe_Meis.summary(solution=micom_solution)
        try:
            metoh_ex = summary.consuming_flux.loc['_EX_METOH_Meis', 'percent']
        except KeyError as e:
            print(e)
            metoh_ex = 0

        with open(com.results_path + 'metoh_ex.csv', 'a+') as file:
            pd.DataFrame(
                [[sample, metoh_ex, max_del]],
                columns=['Sample', 'Exchanged Metanol', 'N Blocked']
                ).to_csv(file, index=False, header=False)

        biomass_fluxes = micom_solution.fluxes.index.str.contains('biomass')
        biomass = 0
        for tax in micom_solution.fluxes.loc[biomass_fluxes].index:
            biomass += micom_solution[tax] * com.models[sample].metabolites.get_by_id(tax.replace('DM_', '')).cmol_by_mol
        biomass = round(biomass / carbon_uptake, 2)

        with open(com.results_path + 'yields.csv', 'a+') as file:
            pd.DataFrame(
                [[sample, 'Biomass', biomass, max_del]],
                columns=['Sample', 'Metabolite', 'Yield', 'N Blocked']
                ).to_csv(file, index=False, header=False)

        for rxn_id in yield_mets.keys():
            met_yield = round(abs(micom_solution[rxn_id] / carbon_uptake), 2)
            with open(com.results_path + 'yields.csv', 'a+') as file:
                pd.DataFrame(
                    [[sample, yield_mets[rxn_id], met_yield, max_del]],
                    columns=['Sample', 'Metabolite', 'Yield', 'N Blocked']
                    ).to_csv(file, index=False, header=False)

# %%codecell
ic_long = pd.read_csv(results_path + 'simulations.csv')
px.box(ic_long.query("Reciver=='Meis'"), x='Giver', y='Interaction Coefficient', template='none')
px.box(ic_long.query("Reciver=='Hyum'"), x='Giver', y='Interaction Coefficient', template='none')
px.box(ic_long.query("Reciver=='Meus'"), x='Giver', y='Interaction Coefficient', template='none')
