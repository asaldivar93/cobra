import pandas as pd
import picrust_parser as pparser
from cobra.flux_analysis import find_leaks, find_leak_mode
pc = pparser.p_model()
model = c_model.copy()

# %% codecell
# change constraints for rxns already corrected
from data_files.corrected_datasets import corrected_revesibility
for r in corrected_revesibility.keys():
    lower_bound = corrected_revesibility[r]['lower_bound']
    upper_bound = corrected_revesibility[r]['upper_bound']
    rxn = model.reactions.get_by_id(r)
    rxn.lower_bound = lower_bound
    rxn.upper_bound = upper_bound

rxn = model.reactions.get_by_id('|RXN-16650|')
rxn.add_metabolites({'|C6|_pe': 1, '|C6|_cy': -1})
rxn.build_reaction_string()

# %% codecell
# find leaks and leak modes of the Model
leaks = find_leaks(model)
leak_modes = find_leak_mode(model, leaks, cutoff_mult = 10)

# %% codecell
# sort reactions in leak_mode for leak_i
leak = '|ACET|_cy'
mode = []
for r, flux in leak_modes[leak]:
    rxn = model.reactions.get_by_id(r)
    mode.append([rxn.subsystem, rxn.id, rxn.build_reaction_string(), flux])
mode.sort()
for subs, id, r, flux in mode:
    print(subs, id, r, flux)

# %% codecell
pd.set_option('display.max_rows', 150)
pc.search_biomass_components(model, path_to_biomass = 'data_files/biomass.csv')
dup_filter = pc.biomass_in_model['component'].duplicated()
pc.biomass_in_model[dup_filter]
pc.biomass_in_model.drop(pc.biomass_in_model[dup_filter].index, inplace = True)
pc.biomass_in_model

# %% codecell
media = ['|CH4|_ex', '|OXYGEN-MOLECULE|_ex', '|NITRATE|_ex']
sinks = ['|WATER|_pe', '|WATER|_cy', '|ACYL-ACP|_cy', '|Corrinoid-Adenosyltransferases|_cy', '|CPD-17931|_pe']
# not_sinks = ['|CPD-20903|_cy', '|FADH2|_cy', '|FAD|_cy', '|CPD-15999|_cy', '|ADENOSYL-HOMO-CYS|_cy', '|FE+2|_cy', '|UNDECAPRENYL-DIPHOSPHATE|_pe', '|D-ALANINE|_pe', '|PROTON|_pe']
not_sinks = ['|CPD-20903|_cy', '|FE+2|_cy', '|HS|_cy', '|P3I|_cy']
sinks.extend(not_sinks)
demands = [['|CPD-15999|_cy', '|UNDECAPRENYL-DIPHOSPHATE|_pe', '|ACP|_cy', '|FORMATE|_cy']]
blocked, available, consistent, inconsistent = pc.test_blocked_components(model, media, demands = demands, sinks = sinks, consistency_check = True)
available


# %% codecell
len(consistent)
len(inconsistent)
i_pth = {}
for pthwy in pc.pathways:
    rxn_in_pathwy = []
    i_pth[pthwy] = []
    for rxn in model.reactions:
        if pthwy in rxn.subsystem:
            rxn_in_pathwy.extend([rxn.id])
    for rxn in rxn_in_pathwy:
        if rxn in inconsistent:
            i_pth[pthwy].extend([rxn])
