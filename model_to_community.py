import pandas as pd
import p_tools
import os
from tqdm import tqdm
import cobra
from copy import deepcopy
from data_files.corrected_datasets import (multi_comp_rxns,
                                           curated_rxns,
                                           exchange_rxns,
                                           added_if_met,
                                           added_pthwys,
                                           sinks,
                                           true_DM,
                                           possible_product,
                                           corrected_revesibility)
import picrust_parser as pparser
artificial_DM = []
artificial_DM.extend(true_DM)

# %% codecell
pt = p_tools.picrust_tools()
pwy_strat = pt.sample_strat_pathway()

taxonomy = pd.read_csv('/home/alexis/UAM/cobra/data_files/silva_taxonomy.csv', sep ='\t')
taxonomy.set_index('sequence_identifier', inplace = True)
seq_counts = pd.read_csv('/home/alexis/UAM/cobra/data_files/otu-table-corrected.tsv', sep = '\t', index_col = 'OTU_ID')
taxonomy.fillna('unknown', inplace = True)
otus = seq_counts['CIR_19'].loc[seq_counts['CIR_19'] != 0].to_frame()

for seq in otus.index:
    otus.loc[seq, 'Family'] = taxonomy.loc[seq, 'Family']
    otus.loc[seq, 'Genus'] = taxonomy.loc[seq, 'Genus']

pwys_in_gens = {}
for genus in otus['Genus'].unique():
    gen_filter = otus.loc[otus['Genus'] == genus].index
    gen_pwy = pwy_strat['CIR_19'].loc[:, gen_filter].dropna(how = 'all').index.to_list()
    not_minimal = ['added']
    while not_minimal:
        not_minimal = []
        for pwy in gen_pwy:
            pwy_id = '|' + pwy + '|'
            if pwy_id in pparser.all_pathways:
                sub_pathways = pparser.metacyc_db[pwy_id].sub_pathways
                if sub_pathways:
                    gen_pwy.remove(pwy)
                    pwys_to_add = [pwy.replace('|', '') for pwy in sub_pathways]
                    gen_pwy.extend(pwys_to_add)
                    not_minimal.append(['added'])
    gen_pwy = list(dict.fromkeys(gen_pwy))
    pwys_in_gens[genus] = gen_pwy

for genus in pwys_in_gens.keys():
    print(genus, len(pwys_in_gens[genus]))

# %% codecell
os.chdir('/home/alexis/UAM/cobra/')
model = cobra.Model()
pc = pparser.p_model()
pathways = pwys_in_gens['Methylocystis'].copy()
pathways.extend(added_pthwys)
gen_id = genus[:2] + genus[-2:]
model.name = gen_id
for i in tqdm(range(int(len(pathways)))):
    pthwy = pathways[i]
    pthwy_id = '|' + pthwy + '|'
    if pthwy_id in pparser.all_pathways:
        pc.add_pathway(model, pthwy_id)
    else:
        pc.unmatched_pthwys.extend([pthwy_id])

model_mets = [met.id for met in model.metabolites]
to_add_filter = added_if_met['metabolite'].isin(model_mets)
pwys_to_extend = [path for path in added_if_met.loc[to_add_filter, 'pathway']
                  if path not in pathways]

for i in tqdm(range(int(len(pwys_to_extend)))):
    pthwy = pwys_to_extend[i]
    pthwy_id = '|' + pthwy + '|'
    if pthwy_id in pparser.all_pathways:
        pc.add_pathway(model, pthwy_id)
    else:
        pc.unmatched_pthwys.extend([pthwy_id])
c_model = deepcopy(model)
# %% codecell
model = deepcopy(c_model)
biomass_in_model = pc.search_biomass_components(model)
model = pc.add_biomass_rxn(model, biomass_in_model)

print(
    genus,
    model.metabolites.get_by_id('biomass').formula_weight,
    model.metabolites.get_by_id('biomass').formula
)

for rxn_id in curated_rxns.keys():
    if rxn_id not in model.reactions:
        subsystem = curated_rxns[rxn_id]['pathway']
        stoichiometry = curated_rxns[rxn_id]['stoichiometry']
        reversible = curated_rxns[rxn_id]['reversible']
        direction = curated_rxns[rxn_id]['direction']
        pc.add_rxn_from_stoichiometry(model, subsystem, rxn_id, stoichiometry, reversible, direction)

for rxn_id in multi_comp_rxns.keys():
    if rxn_id not in model.reactions:
        subsystem = multi_comp_rxns[rxn_id]['pathway']
        stoichiometry = multi_comp_rxns[rxn_id]['stoichiometry']
        reversible = multi_comp_rxns[rxn_id]['reversible']
        direction = multi_comp_rxns[rxn_id]['direction']
        pc.add_rxn_from_stoichiometry(model, subsystem, rxn_id, stoichiometry, reversible, direction)
        pc.multi_compartment_rxns.append([subsystem, rxn_id, 'added'])

for met in sinks:
    try:
        sink = model.metabolites.get_by_id(met)
    except:
        pass
    else:
        model.add_boundary(
            sink, type = 'sink'
        )

for met in artificial_DM:
    try:
        dm = model.metabolites.get_by_id(met)
    except:
        pass
    else:
        model.add_boundary(
            dm, type = 'demand'
        )

for ex_rxn in exchange_rxns.keys():
    if ex_rxn not in model.reactions:
        subsystem = exchange_rxns[ex_rxn]['pathway']
        stoichiometry = exchange_rxns[ex_rxn]['stoichiometry']
        reversible = exchange_rxns[ex_rxn]['reversible']
        direction = exchange_rxns[ex_rxn]['direction']
        pc.add_rxn_from_stoichiometry(
            model, subsystem, ex_rxn, stoichiometry, reversible, direction
            )

class_filter = biomass_in_model['class'].isin(
    ['Amino_Acids', 'Cofactors', 'Intracellular_Metabolites', 'Carbohydrates']
    )

to_exchange = biomass_in_model.loc[class_filter, 'met_id'].to_list()
to_exchange.extend(possible_product)
for met in to_exchange:
    try:
        metabolite = model.metabolites.get_by_id(met)
    except:
        pass
    else:
        met_to_ex = met[:-3] + '_ex'

        stoichiometry = {}
        stoichiometry[met_to_ex] = -1
        stoichiometry[met] = 1

        met_to_add = cobra.Metabolite(
            met_to_ex,
            name = met_to_ex[:-3].replace('|', ''),
            formula = metabolite.formula,
            charge = metabolite.charge,
            compartment = 'extracellular'
            )

        rxn_id = 'EX_' + met_to_ex[:-3].replace('|', '')
        rxn_to_add = cobra.Reaction(
            rxn_id,
            upper_bound = 1000,
            lower_bound = -1000
            )

        model.add_metabolites(
            met_to_add
        )

        model.add_reaction(
            rxn_to_add
        )

        rxn_to_add.add_metabolites(
            stoichiometry
        )

        model.add_boundary(
            met_to_add,
            type = 'demand'
        )

    for r in corrected_revesibility.keys():
        lower_bound = corrected_revesibility[r]['lower_bound']
        upper_bound = corrected_revesibility[r]['upper_bound']
        try:
            rxn = model.reactions.get_by_id(r)
        except:
            pass
        else:
            rxn.lower_bound = lower_bound
            rxn.upper_bound = upper_bound

    remove_mets = []

    for met in model.metabolites:
        if not met.reactions:
            remove_mets.extend([met])
    model.remove_metabolites(remove_mets)

for met in model.metabolites:
    met.id = met.id + '_' + gen_id

for rxn in model.reactions:
    rxn.id = rxn.id + '_' + gen_id

model.repair()


# %% codecell
ag_model = build_aggregate_model(model)

ag_model.constraints.Amino_Acids_Paas
