import cobra
import pandas as pd
import os
from copy import deepcopy

os.chdir('/home/alexis/UAM/cobra')

# %%codecell
taxonomy = pd.read_csv('data_files/silva_taxonomy.csv', sep ='\t')
taxonomy.set_index('sequence_identifier', inplace = True)
seq_counts = pd.read_csv('data_files/otu-table-corrected.tsv', sep = '\t', index_col = 'OTU_ID')
taxonomy.fillna('unknown', inplace = True)
otus = seq_counts['CIR_19'].loc[seq_counts['CIR_19'] != 0].to_frame()

for seq in otus.index:
    otus.loc[seq, 'Family'] = taxonomy.loc[seq, 'Family']
    otus.loc[seq, 'Genus'] = taxonomy.loc[seq, 'Genus']

taxa_database = pd.DataFrame()
for genus in otus.loc[:, 'Genus'].unique():
    gen_id = genus[:2] + genus[-2:]
    model_path = 'models/CIR_19_complete/{}.xml'.format(gen_id)
    gen_filter = otus['Genus'] == genus
    taxa_database.loc[gen_id, 'seq_counts'] = otus.loc[gen_filter, 'CIR_19'].sum()
    taxa_database.loc[gen_id, 'Genus'] = genus
    taxa_database.loc[gen_id, 'model_path'] = model_path

for k in taxa_database.index:
    taxa_database.loc[k, 'relative_abundance'] = (taxa_database.loc[k, 'seq_counts'] / taxa_database.loc[:, 'seq_counts'].sum())
taxa_database
# %%codecell
universal = cobra.io.read_sbml_model('models/cir.xml')
universal.solver = 'cplex'
ex_model = cobra.Model()
for met in universal.metabolites:
    if '_ou' not in met.id and '_ex' not in met.id:
        has_exchange = ['EX_' in rxn.id for rxn in met.reactions]
        if not any(has_exchange):
            metabolite = met

            met_to_ex = met.id[:-3] + '_ex'

            stoichiometry = {}
            stoichiometry[met_to_ex] = -1
            stoichiometry[met.id] = 1

            met_to_add = cobra.Metabolite(
                met_to_ex,
                name = met_to_ex[:-3],
                formula = metabolite.formula,
                charge = metabolite.charge,
                compartment = 'extracellular'
                )

            ex_model.add_metabolites([metabolite, met_to_add])

            rxn_id = '_EX' + met_to_ex
            rxn_to_add = cobra.Reaction(
                rxn_id,
                upper_bound = 1000,
                lower_bound = -1000
                )

            ex_model.add_reactions(
                [rxn_to_add]
            )

            ex_model.reactions.get_by_id(rxn_id).add_metabolites(
                stoichiometry
            )
# %%codecell
for met in ex_model.metabolites:
    if '_ex' not in met.id:
        met.id = met.id + '_Meis'
        ex_model.repair()

for rxn in ex_model.reactions:
    rxn.id = rxn.id + '_Meis'
    ex_model.repair()
cobra.io.write_sbml_model(ex_model, 'models/universal/ex_Meis.xml')
# %%codecell
ex_meis = cobra.io.read_sbml_model('models/universal/ex_Meis.xml')
for taxa_k in taxa_database.index:
    new_model = deepcopy(universal)
    new_model.merge(ex_model)
    for met in new_model.metabolites:
        if '_ou' in met.id or '_ex' in met.id:
            pass
        else:
            met.id = met.id + '_' + taxa_k
        new_model.repair()
    for rxn in new_model.reactions:
        if 'ou_' in rxn.id or ('DM_' and '_ex') in rxn.id:
            pass
        else:
            rxn.id = rxn.id + '_' + taxa_k
        new_model.repair()
    new_model.merge(ex_meis)
    cobra.io.write_sbml_model(new_model, 'models/universal/{}.xml'.format(taxa_k))
