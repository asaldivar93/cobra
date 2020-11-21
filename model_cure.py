import pandas as pd
import picrust_parser as pparser

from cobra.flux_analysis import (find_leaks,
                                 find_leak_mode,
                                 )
from cobra.flux_analysis import (find_blocked_reactions,
                                 find_blocked_mets,
                                 find_mets_to_connect,
                                 find_dead_ends
                                 )
from data_files.corrected_datasets import corrected_revesibility


# %% codecell
# change constraints for rxns already corrected
del model
model = c_model.copy()

for r in corrected_revesibility.keys():

    lower_bound = corrected_revesibility[r]['lower_bound']
    upper_bound = corrected_revesibility[r]['upper_bound']
    rxn = model.reactions.get_by_id(r)
    rxn.lower_bound = lower_bound
    rxn.upper_bound = upper_bound

remove_mets = []
for met in model.metabolites:
    if not met.reactions:
        remove_mets.extend([met])

model.remove_metabolites(remove_mets)

# %% codecell
# (A) Correct the reversibility of reactions to remove stoichimetrically balanced cycles.
# (1) To do this we'ĺl first find all leaks and leak modes of the Model. The reactions
# found in a leak_mode will need a deeper look to asses reversibility. Note: not all reactions
# in the leak_mode are incorrect

leaks = find_leaks(model)
leak_modes = find_leak_mode(model, leaks, cutoff_mult = 10)
leaks
# %% codecell
# (1.1) Print reactions in leak_mode for leak_i. It is easier to pinpoint incorrect reactions
# if they are seen as components of a pathway rather than as individual reactions.
leak = '|FUM|_cy'
mode = []
for r, flux in leak_modes[leak]:
    rxn = model.reactions.get_by_id(r)
    mode.append([rxn.subsystem, rxn.id, rxn.build_reaction_string(), flux])
    # mode.extend([rxn.id])

mode.sort()
for subs, id, r, flux in mode:
    print(subs, id, r, flux)

# Correct reversivility of reactions and repeat (1) iteratevly until there are no leaks left

# %% codecell
# (Optional) Finding inconsistent reactions in the subset of reactions in leak_mode
# may help to pinpoint mistakes.
w_model = model.copy()
rxns_to_remove = []
for rxn in w_model.reactions:
    if rxn.id in mode:
        pass
    else:
        rxns_to_remove.extend([rxn])
w_model.remove_reactions(rxns_to_remove)
consistent, inconsistent = pc.fastcc(w_model)

# %% codecell
# (3) Identify which biomass compenents are metabolites in the model
pd.set_option('display.max_rows', 150)
biomass_in_model = pc.search_biomass_components(model)
biomass_in_model.drop(index = pc.biomass_in_model[pc.biomass_in_model.loc[:, 'met_id'] == '|ATP|_cy'].index, inplace = True)
biomass_in_model.drop(index = pc.biomass_in_model[pc.biomass_in_model.loc[:, 'met_id'] == '|PPI|_cy'].index, inplace = True)
biomass_in_model = biomass_in_model['met_id'].to_list()

# %% codecell
# (4) Once there are no leaks left in the model we can identify which metabolites
# can be produced by the model. During this step we'll identify blocked reactions
# and dead-end metabolites in an iterative process

# (4.1) Start by defining componets of the media, add or remove sinks iteratevly as needed
media = ['|CH4|_ex', '|OXYGEN-MOLECULE|_ex', '|NITRATE|_ex', '|FE+2|_ex', '|Pi|_ex', '|SULFATE|_ex']

sinks = ['|WATER|_pe', '|WATER|_cy', '|ACP|_cy',
         '|CPD-1302|_cy', '|CPD-1301|_cy', '|Guanine34-in-tRNAs|_cy',
         '|S-ADENOSYL-4-METHYLTHIO-2-OXOBUTANOATE|_cy', '|Hpr-pi-phospho-L-histidines|_cy',
         '|Unsulfurated-Sulfur-Acceptors|_cy', '|Corrinoid-Adenosyltransferases|_cy',
         '|LysW-C-Terminal-L-Glutamate|_cy', '|CPD-17931|_pe',
         '|CoI-Corrinoid-Fe-S-proteins|_cy', '|Cysteine-Desulfurase-L-cysteine|_cy',
         '|D-alanine-carrier-protein|_cy', '|DsrE3A-L-cysteine|_cy', '|Ox-Thioredoxin|_cy',
         '|Sulfur-Carrier-Proteins-ThiI|_cy', '|Thi-S|_cy', '|4Fe-4S+1|_cy',
         '|CPD-381|_cy'
         ]

artificial_EX = ['|NA+|_cy', '|MG+2|_cy', '|CO+2|_cy', '|CL-|_cy']

true_DM = ['|NA+|_pe', '|DTDP-RHAMNOSE|_cy', '|ACP|_pe', '|FORMAMIDE|_cy',
           '|UNDECAPRENYL-DIPHOSPHATE|_pe', '|CPD-10640|_cy', '|UNKNOWN|_cy',
           '|N-ACETYL-D-GLUCOSAMINE|_cy', '|CPD-15999|_cy',
           '|ETHANOL-AMINE|_cy', '|ALLYSINE|_cy', '|Alcohols|_cy',
           '|1-AMINO-PROPAN-2-OL|_cy', '|P3I|_cy', '|Cysteine-Desulfurase-L-cysteine|_cy',
           '|GLYCOLALDEHYDE|_cy', '|4Fe-4S+2|_cy', '|HYDROGEN-PEROXIDE|_cy',
           '|Hpr-Histidine|_cy', '|PPI|_cy', '|Pi|_cy', '|NAD|_cy', '|NADP|_cy'
           ]

possible_product = ['|NITROGEN-MOLECULE|_pe', '|HYDROGEN-MOLECULE|_cy', '|ETOH|_cy',
                    '|CPD-10755|_cy', '|NITROGEN-MOLECULE|_cy', '|Methylketones|_cy',
                    '|CPD-347|_cy', '|PROPIONATE|_cy', '|PROPANE-1-2-DIOL|_cy',
                    '|CPD-10353|_cy', '|BUTANEDIOL|_cy', '|BUTANOL|_cy', '|ACETONE|_cy',
                    '|FORMATE|_cy', 'biomass', '|CARBON-DIOXIDE|_cy', '|ACET|_cy',
                    '|BUTYRIC_ACID|_cy', '|CARBON-MONOXIDE|_cy', '|PUTRESCINE|_cy',
                    '|Poly-Hydroxybutyrate|_cy', '|PROTON|_cy'
                    ]

artificial_DM = []
sinks.extend(artificial_EX)
artificial_DM.extend(true_DM)
artificial_DM.extend(possible_product)
artificial_DM.extend(pc.biomass_in_model['met_id'].to_list())
demands = []
with model as model:
    for substrate in media:
        met = model.metabolites.get_by_id(substrate)
        model.add_boundary(
            met, type = 'exchange'
            )
    model.medium = {rxn.id: 100 for rxn in model.exchanges}

    for met in sinks:
        sink = model.metabolites.get_by_id(met)
        model.add_boundary(
            sink, type = 'sink'
            )

    for met in artificial_DM:
        dm = model.metabolites.get_by_id(met)
        model.add_boundary(
            dm, type = 'demand'
            )
    # Identify blocked reactions
    blocked_rxns = find_blocked_reactions(
        model
        )
    # Find a set of possible dead end metabolites
    dead_ends = find_dead_ends(
        model, carbon_source = '|CH4|_ex'
        )
    # Find a set of blocked metabolites
    blocked_mets, available_mets, produced_from_nothing = find_blocked_mets(
        model, demands = demands, carbon_source = '|CH4|_ex'
        )

# (4.2) Loop through dead end metabolites to classify mets in:
#   Alternative Carbon Sources
#   Posible Products
#   Only produced mets that require artificial demand reactions
#   True sinks of the model (Note: be carefull not to allow artificial carbon production)
#   Sinks required by the model that need correction
#   True dead ends of the model that need to be removed or connected
# Repeat (4.1) until dead_ends is empty
print(len(dead_ends))
print(len(blocked_mets))
print(len(blocked_rxns))

# %% codecell
# (4.3) Once all dead end metabolites have been curated, find the metabolites required
# for production of all other metabolites
with model as model:
    for substrate in media:
        met = model.metabolites.get_by_id(substrate)
        model.add_boundary(
            met, type = 'exchange'
            )
    model.medium = {rxn.id: 100 for rxn in model.exchanges}

    for met in sinks:
        sink = model.metabolites.get_by_id(met)
        model.add_boundary(
            sink, type = 'sink'
            )

    for met in artificial_DM:
        dm = model.metabolites.get_by_id(met)
        model.add_boundary(
            dm, type = 'demand'
            )

    mets_to_connect = find_mets_to_connect(
        model, blocked_mets, carbon_source = '|CH4|_ex', n_workers = 14
        )
# (4.4) Loop through mets_to_connect and associated reactions to identify additional
# required sinks, correct miss annotated stoichiometry and reversibility constraints,
# and find missing reactions. Repeat (4.3) until all mets can be produced
mets_to_connect

# %% codecell
# (5) Check that all biomass components can be produced by the model.
artificial_DM = []
artificial_DM.extend(true_DM)
artificial_DM.extend(possible_product)
demands = pc.biomass_in_model['met_id'].to_list()
# demands = ['|IMP|_cy']
with model as model:
    for substrate in media:
        met = model.metabolites.get_by_id(substrate)
        model.add_boundary(
            met, type = 'exchange'
            )
    medium = {rxn.id: 1000 for rxn in model.exchanges}
    medium['EX_|CH4|_ex'] = 100
    model.medium = medium

    atp_dm = model.reactions.get_by_id('DM_|ATP|_cy')
    atp_dm.lower_bound = 3.5
    atp_dm.upper_bound = 1000

    o2_ex = model.reactions.get_by_id('O2')
    o2_ex.upper_bound = 1000
    o2_ex.lower_bound = 0

    for met in sinks:
        sink = model.metabolites.get_by_id(met)
        model.add_boundary(
            sink, type = 'sink'
            )

    for met in artificial_DM:
        dm = model.metabolites.get_by_id(met)
        model.add_boundary(
            dm, type = 'demand'
            )

    for r in ['|ATP|_cy_syn', '|1.10.2.2-RXN|', '|CYTOCHROME-C-OXIDASE-RXN|', '|NADH-DEHYDROG-A-RXN|']:
        rxn = model.reactions.get_by_id(r)
        rxn.upper_bound = 1000

    # Find the set of blocked metabolites
    blocked_mets, available_mets, produced_from_nothing = find_blocked_mets(
        model, demands = demands, carbon_source = '|CH4|_ex'
        )
    # Once all biomass components can be produced, it is necesary to check for leaks once more
    leaks = find_leaks(
        model, keep_boundaries = True
        )

# (5.1) Loop through blocked_mets and associated reactions, modifying constraints
# and stoichiometry or adding pathways until all biomass components can be produced from CH4
# blocked_mets
blocked_mets
# %% codecell
# (5.2) It is possible to create new stoichiometrich balanced cycles with current
# configuration of sinks in the model. to find leaks, remove CH4 from media, then
# use find_leak_mode on leaks to find which sinks and reactions are causing trouble
media = ['|OXYGEN-MOLECULE|_ex', '|NITRATE|_ex', '|FE+2|_ex', '|Pi|_ex', '|SULFATE|_ex']
leaks = []

with model as model:
    for substrate in media:
        met = model.metabolites.get_by_id(substrate)
        model.add_boundary(
            met, type = 'exchange'
            )
    medium = {rxn.id: 1000 for rxn in model.exchanges}
    model.medium = medium

    atp_dm = model.reactions.get_by_id('DM_|ATP|_cy')
    atp_dm.lower_bound = 3.5
    atp_dm.upper_bound = 1000

    o2_ex = model.reactions.get_by_id('O2')
    o2_ex.upper_bound = 1000
    o2_ex.lower_bound = 0

    for met in sinks:
        sink = model.metabolites.get_by_id(met)
        model.add_boundary(
            sink, type = 'sink'
            )

    for met in artificial_DM:
        dm = model.metabolites.get_by_id(met)
        model.add_boundary(
            dm, type = 'demand'
            )

    for r in ['|ATP|_cy_syn', '|1.10.2.2-RXN|', '|CYTOCHROME-C-OXIDASE-RXN|', '|NADH-DEHYDROG-A-RXN|']:
        rxn = model.reactions.get_by_id(r)
        rxn.upper_bound = 1000
    leak_modes = find_leak_mode(
        model, leaks = leaks, cutoff_mult = 1, keep_boundaries = True
        )

# %% codecell
# Loop through the reactions in leak_model until no metabolites can be produced
# without methane
leak = '|CARBON-DIOXIDE|_cy'
mode = []
super_inconsistent = []
for r, flux in leak_modes[leak]:
    rxn = model.reactions.get_by_id(r)
    if flux < 0:
        mode.append([rxn.subsystem, rxn.id, rxn.build_reaction_string(), flux])
    # if rxn.id in inconsistent:
    # super_inconsistent.extend([rxn.id])
    # mode.extend([rxn.id])

mode.sort()
for subs, id, r, flux in mode:
    # if id not in corrected_revesibility.keys(): # and id not in checked: #and id in inconsistent:
    print(subs, id, r, flux)

pc.rxns_of_metabolite(model, 'CPD-12575')
pc.rxns_in_pthwy(model, '|PWY-7343|')

pparser.metacyc_db['PWY-5433']
model.compartments
