# missing method to claculate dg of transport reactions and multicompartment reactions
# missing sparsefba
# %% codecell
from tqdm import tqdm
import pandas as pd
import cobra
from data_files.corrected_datasets import (multi_comp_rxns,
                                           curated_rxns,
                                           exchange_rxns,
                                           added_if_met,
                                           sinks,
                                           true_DM,
                                           possible_product)
import picrust_parser as pparser
artificial_DM = []
artificial_DM.extend(true_DM)

# %% codecell
model = cobra.Model()
pc = pparser.p_model()
for ex_rxn in exchange_rxns.keys():
    if ex_rxn not in model.reactions:
        subsystem = exchange_rxns[ex_rxn]['pathway']
        stoichiometry = exchange_rxns[ex_rxn]['stoichiometry']
        reversible = exchange_rxns[ex_rxn]['reversible']
        direction = exchange_rxns[ex_rxn]['direction']
        pc.add_rxn_from_stoichiometry(model, subsystem, ex_rxn, stoichiometry, reversible, direction)

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

for i in tqdm(range(int(len(pc.pathways)))):
    pthwy = pc.pathways[i]
    pthwy_id = '|' + pthwy + '|'
    if pthwy_id in pparser.all_pathways:
        pc.add_pathway(model, pthwy_id)
    else:
        pc.unmatched_pthwys.extend([pthwy_id])

model_mets = [met.id for met in model.metabolites]
to_add_filter = added_if_met['metabolite'].isin(model_mets)
pwys_to_extend = [path for path in added_if_met.loc[to_add_filter, 'pathway']
                  if path not in pc.pathways]

for i in tqdm(range(int(len(pwys_to_extend)))):
    pthwy = pwys_to_extend[i]
    pthwy_id = '|' + pthwy + '|'
    if pthwy_id in pparser.all_pathways:
        pc.add_pathway(model, pthwy_id)
    else:
        pc.unmatched_pthwys.extend([pthwy_id])

print(str(len(pc.unmatched_pthwys)) + ' pathways do not match')
print(str(len(pc.mets_w_concentration)) + ' metabolites have concentration info')
print(str(len(pc.rxns_w_thermo_error)) + ' reactions got an error')
print(str(len(pc.multi_compartment_rxns)) + ' reactions are multicompartment reactions')
print(str(len(pc.should_be_balanced_rxns)) + ' reactions are wrongly annotated')
print(str(len(pc.one_proton_unbalanced_rxns)) + ' reactions are one H+ unbalanced')
print(str(len(pc.proton_unbalanced_rxns)) + ' reactions are H+ unbalanced')
print(str(len(pc.one_charge_unbalanced_rxns)) + ' reactions are one charge unbalanced')
print(str(len(pc.charge_unbalanced_rxns)) + ' reactions are charge unbalanced')
print(str(len(pc.proton_and_charge_unbalanced_rxns)) + ' reactions are H+ and charge unbalanced')
print(str(len(pc.eqtr_unbalanced_rxns)) + ' reactions are unbalanced on equilibrator')
print(str(len(pc.unbalanced_rxns)) + ' reactions are unbalanced')
print(str(len(pc.corrected_hc_rxns)) + ' reactions had stoichiometry change by nH+')
cobra.io.save_matlab_model(model, 'model.mat')

c_model = model.copy()
model = c_model.copy()

biomass_in_model = pc.search_biomass_components(model)
model = pc.add_biomass_rxn(model, biomass_in_model)

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

class_filter = biomass_in_model['class'].isin(
    ['Amino_Acids', 'Cofactors', 'Intracellular_Metabolites', 'Carbohydrates']
    )

to_exchange = biomass_in_model.loc[class_filter, 'met_id'].to_list()
to_exchange.extend(possible_product)
for met in to_exchange:
    metabolite = model.metabolites.get_by_id(met)
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

c_model = model.copy()
