import os
import pandas as pd
os.chdir('/home/alexis/UAM/cobra/')
import cobra
import warnings
warnings.filterwarnings("ignore")
from cobra.flux_analysis import micom, add_steadycom_constraints, plot_micom_tradeoff
from tools.cobra_tools import add_boundaries_for_simulation

# %% codecell
taxonomy = pd.read_csv('data_files/silva_taxonomy.csv', sep ='\t')
taxonomy.set_index('sequence_identifier', inplace = True)
seq_counts = pd.read_csv('data_files/otu-table-corrected.tsv', sep = '\t', index_col = 'OTU_ID')
seq_counts.loc['0281ce539a44c8b1b059ae85c15e3fab', 'My_20'] = 0
taxonomy.fillna('unknown', inplace = True)
samples = ['My_20', 'CIR_19']
otus = pd.DataFrame()
for sample in samples:
    sample_counts = seq_counts[sample].loc[seq_counts[sample] != 0].to_frame().rename(columns = {sample: 'seq_counts'}).reset_index()
    sample_counts.loc[:, 'sample'] = sample
    otus = otus.append(sample_counts)

for seq in otus['OTU_ID'].unique():
    seq_filter = otus.loc[:, 'OTU_ID'] == seq
    otus.loc[seq_filter, 'Family'] = taxonomy.loc[seq, 'Family']
    otus.loc[seq_filter, 'Genus'] = taxonomy.loc[seq, 'Genus']


taxa_database = pd.DataFrame(columns = ['id', 'seq_counts', 'sample', 'Genus', 'model_path'])
sample = 'My_20'
genus = 'uncultured'
for sample in samples:
    sample_filter = otus.loc[:, 'sample'] == sample
    new_otus = otus.loc[sample_filter, :].set_index('OTU_ID')
    for genus in new_otus.loc[:, 'Genus'].unique():
        if genus == 'uncultured' or genus == 'unknown':
            gen_filter = new_otus.loc[new_otus['Genus'] == genus].index
            families = new_otus.loc[gen_filter, :]
            for family in families['Family'].unique():
                gen_filter = families.loc[families['Family'] == family].index
                gen_id = family[:2] + family[-2:]
                model_path = 'models/models_gapfilled/{}.xml'.format(gen_id)
                taxa_database = taxa_database.append(pd.DataFrame(
                    [[gen_id,
                     families.loc[gen_filter, 'seq_counts'].sum(),
                     sample,
                     family,
                     model_path]],
                    columns = ['id', 'seq_counts', 'sample', 'Genus', 'model_path']
                    )
                )

        else:
            gen_filter = new_otus['Genus'] == genus
            gen_id = genus[:2] + genus[-2:]
            model_path = 'models/models_gapfilled/{}.xml'.format(gen_id)
            taxa_database = taxa_database.append(pd.DataFrame(
                [[gen_id,
                 new_otus.loc[gen_filter, 'seq_counts'].sum(),
                 sample,
                 genus,
                 model_path]],
                columns = ['id', 'seq_counts', 'sample', 'Genus', 'model_path']
                )
            )

# taxa_database = taxa_database.loc[['Meis', 'Meus'], :]
# taxa_database = taxa_database.loc[~taxa_database.index.isin(['Meus', 'Teas', 'Tala']), :]

taxa_database.reset_index(inplace = True)
for sample in samples:
    sample_filter = taxa_database.loc[:, 'sample'] == sample
    for k in taxa_database.loc[sample_filter, :].index:
        taxa_database.loc[k, 'relative_abundance'] = (taxa_database.loc[k, 'seq_counts'] / taxa_database.loc[sample_filter, 'seq_counts'].sum())
taxa_database
# %% codecell
community_model = cobra.Model()
community_model.solver = 'cplex'
sample = taxa_database.loc[:, 'sample'] == 'My_20'
for k in taxa_database.loc[sample].index:
    new_model = cobra.io.read_sbml_model(taxa_database.loc[k, 'model_path'])
    community_model.merge(new_model)

# %% codecell
media = ['_CH4_ou', '_OXYGEN_MOLECULE_ou', '_NITRATE_ou', '_FE_2_ou',
         '_Pi_ou', '_SULFATE_ou', '_NA__ou', '_MG_2_ou', '_CO_2_ou',
         '_CL__ou']
with community_model as model:
    add_boundaries_for_simulation(model, media)
    objective = {reaction: taxa_database.loc[reaction.id[-4:], 'relative_abundance'] for reaction in model.reactions if 'BIOMASS' in reaction.id}
    model.objective = objective
    solution = model.optimize()
    print(model.summary(solution = solution))

# %% codecell
media = ['_CH4_ou', '_OXYGEN_MOLECULE_ou', '_NITRATE_ou', '_FE_2_ou',
         '_Pi_ou', '_SULFATE_ou', '_NA__ou', '_MG_2_ou', '_CO_2_ou',
         '_CL__ou']
sample_db = taxa_database[sample].set_index('id')
with community_model as model:
    add_boundaries_for_simulation(model, media)
    nonzero = plot_micom_tradeoff(model, sample_db)
nonzero
# %% codecell
media = ['_CH4_ou', '_OXYGEN_MOLECULE_ou', '_NITRATE_ou', '_FE_2_ou',
         '_Pi_ou', '_SULFATE_ou', '_NA__ou', '_MG_2_ou', '_CO_2_ou',
         '_CL__ou']
sample_db = taxa_database[sample].set_index('id')
with community_model as model:
    add_boundaries_for_simulation(model, media)
    micom_model, solution, quadratic_model, qs, ls, gr = micom(
        model, 0.2, sample_db
    )
    print(micom_model.summary(solution = solution))
    print(micom_model.metabolites._METOH_ex.summary(solution = solution))
ls.status

# %%codecell
media = ['_CH4_ou', '_OXYGEN_MOLECULE_ou', '_NITRATE_ou', '_FE_2_ou',
         '_Pi_ou', '_SULFATE_ou', '_NA__ou', '_MG_2_ou', '_CO_2_ou',
         '_CL__ou']
sample_db = taxa_database[sample].set_index('id')
with community_model as model:
    add_boundaries_for_simulation(model, media)
    add_steadycom_constraints(model, sample_db)
    problem = model.problem
    growth_rate = model.variables.get('growth_rate')
    objective = problem.Objective(growth_rate, direction = 'max')
    model.objective = objective
    solution = model.optimize()
    print(model.summary(solution = solution))
    print(model.metabolites._METOH_ex.summary(solution = solution))
solution.status

# %% codecel
for rxn in community_model.reactions:
    if '_EX_' in rxn.id:
        if abs(solution[rxn.id]) > 0:
            print(rxn.id, solution[rxn.id])

range(0, 1, 0.1)
