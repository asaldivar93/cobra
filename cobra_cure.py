import pandas as pd
import picrust_parser as pparser
from cobra.flux_analysis import (build_sparse_model,
                                 find_leaks,
                                 find_leak_mode,
                                 find_blocked_reactions
                                 )
from cobra.flux_analysis import (find_blocked_mets,
                                 find_mets_to_connect,
                                 find_dead_ends
                                 )
from data_files.corrected_datasets import corrected_revesibility
pc = pparser.p_model()
model = c_model.copy()

# %% codecell
# change constraints for rxns already corrected
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
sparse_model, met_var = build_sparse_model(model)
leaks = find_leaks(sparse_model, met_var)
leak_modes = find_leak_mode(sparse_model, leaks, cutoff_mult = 10)

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
pc.search_biomass_components(model, path_to_biomass = 'data_files/biomass.csv')
dup_filter = pc.biomass_in_model['component'].duplicated()
pc.biomass_in_model[dup_filter]
pc.biomass_in_model.drop(pc.biomass_in_model[dup_filter].index, inplace = True)
pc.biomass_in_model.drop(index = pc.biomass_in_model[pc.biomass_in_model.loc[:, 'met_id'] == '|ATP|_cy'].index, inplace = True)
pc.biomass_in_model.drop(index = pc.biomass_in_model[pc.biomass_in_model.loc[:, 'met_id'] == '|PPI|_cy'].index, inplace = True)
pc.biomass_in_model['met_id'].to_list()

# %% codecell
# (4) Once there are no leaks left in the model we can identify which metabolites
# can be produced by the model. During this step we'll identify blocked reactions
# and dead-end metabolites in an iterative process


# (4.1) Start by defining componets of the media, add or remove sinks iteratevly as needed
media = ['|CH4|_ex', '|OXYGEN-MOLECULE|_ex', '|NITRATE|_ex', '|FE+2|_ex', '|Pi|_ex', '|SULFATE|_ex']

sinks = ['|WATER|_pe', '|WATER|_cy', '|Menaquinones|_it', '|ACP|_cy',
         '|CPD-1302|_cy', '|CPD-1301|_cy', '|Guanine34-in-tRNAs|_cy', '|CL-|_cy',
         '|CPD-302|_cy', '|Hpr-Histidine|_cy', '|biotin-L-lysine-in-BCCP-dimers|_cy',
         '|carboxybiotin-L-lysine-in-BCCP-dimers|_cy', '|S-ADENOSYL-4-METHYLTHIO-2-OXOBUTANOATE|_cy',
         '|Unsulfurated-Sulfur-Acceptors|_cy', '|Corrinoid-Adenosyltransferases|_cy',
         '|LysW-C-Terminal-L-Glutamate|_cy', '|CPD-17931|_pe', '|Cytochromes-C-Oxidized|_pe',
         '|CoI-Corrinoid-Fe-S-proteins|_cy', '|Cysteine-Desulfurase-L-cysteine|_cy',
         '|D-alanine-carrier-protein|_cy', '|DsrE3A-L-cysteine|_cy', '|Ox-Thioredoxin|_cy',
         '|Sulfur-Carrier-Proteins-ThiI|_cy', '|Thi-S|_cy'
         ]

artificial_EX = ['|NA+|_cy', '|MG+2|_cy', '|CO+2|_cy', '|Fatty-Acids|_cy',
                 '|CPD-12298|_cy', '|CPD-17989|_pe'
                 ]

possible_source = ['|CREATININE|_cy', '|CPD0-1107|_cy', '|D-GALACTARATE|_cy',
                   '|CPD-15633|_cy', '|D-GLUCARATE|_cy', '|3-PHENYLPROPIONATE|_cy',
                   '|3-HYDROXYPHENYL-PROPIONATE|_cy', '|CPD-721|_cy', '|ACETYLENE|_cy',
                   '|PURINE|_cy', '|CPD-148|_cy', '|MYO-INOSITOL|_cy', '|BENZOYLCOA|_cy',
                   '|BETAINE|_cy', '|GLUTARYL-COA|_cy', '|CPD-10663|_cy', '|CPD-10576|_cy',
                   '|2-AMINOPHENOL|_cy', '|CPD-10797|_cy', '|CPD-674|_cy', '|ANDROST4ENE|_cy',
                   '|15-ANHYDRO-D-FRUCTOSE|_cy', '|CPD-3617|_cy', '|L-arabinopyranose|_cy',
                   '|TYRAMINE|_cy', '|CPD-58|_cy', '|1-4-HYDROXYPHENYL-2-METHYLAMINOETHAN|_cy',
                   '|DOPAMINE|_cy', '|CPD0-1068|_cy', '|RS-3-Sulfolactate|_cy',
                   '|PYRIDOXAMINE|_cy', '|PYRIDOXAL|_cy', '|THYMINE|_cy', '|HMP|_cy',
                   '|DEOXYGUANOSINE|_cy', '|AMMONIA|_cy', '|CPD-110|_cy',
                   '|PHENYLETHYLAMINE|_cy', '|Elemental-Sulfur|_cy', '|S2O3|_cy',
                   '|CPD-205|_cy', '|CPD-633|_cy', '|D-SORBITOL-6-P|_cy',
                   '|GALACTITOL|_pe', '|Beta-D-Glucuronides|_cy', '|N-ACETYLNEURAMINATE|_cy',
                   '|CPD-20903|_cy', '|TOLUENE|_cy'
                   ]

true_DM = ['|NA+|_pe', '|DTDP-RHAMNOSE|_cy', '|ACP|_pe', '|FORMAMIDE|_cy',
           '|UNDECAPRENYL-DIPHOSPHATE|_pe', '|CPD-10640|_cy', '|UNKNOWN|_cy',
           '|N-ACETYL-D-GLUCOSAMINE|_cy', '|3-5-ADP|_cy', '|CPD-15999|_cy',
           '|ETHANOL-AMINE|_cy', '|CPD-1091|_cy', '|ALLYSINE|_cy', '|Alcohols|_cy',
           '|1-AMINO-PROPAN-2-OL|_cy', '|P3I|_cy', '|Cysteine-Desulfurase-L-cysteine|_cy'
           ]

possible_product = ['|NITROGEN-MOLECULE|_pe', '|HYDROGEN-MOLECULE|_cy', '|ETOH|_cy',
                    '|CPD-10755|_cy', '|NITROGEN-MOLECULE|_cy', '|Methylketones|_cy',
                    '|CPD-347|_cy'
                    ]

posible_biomass = ['|PROTOHEME|_cy', '|ECTOINE|_cy', '|ADENOSYLCOBALAMIN|_cy',
                   '|PALMITYL-COA|_cy', '|LAUROYLCOA-CPD|_cy', '|CPD-9955|_cy',
                   '|CPD-9957|_cy', '|CPD-9958|_cy', '|URATE|_cy', '|CPD-9247|_pe',
                   '|STEAROYL-COA|_cy', '|CPD-9245|_cy', '|CPD-11994|_cy', '|C3|_cy',
                   '|CPD-12336|_cy', '|tRNAs-with-queuine|_cy', '|CPD-9956|_cy',
                   '|CPD-13167|_cy', '|CPD-10314|_cy', '|TDP-FUC4NAC|_cy',
                   '|11Z-icos-11-enoyl-ACPs|_cy', '|ADP-L-GLYCERO-D-MANNO-HEPTOSE|_cy',
                   '|CPD-12130|_cy', '|CPD-12831|_cy', '|Me-CoM|_cy', '|CPD-1281|_cy',
                   '|CPD-12124|_cy', '|UDP-D-GALACTURONATE|_cy', '|KDO2-LIPID-IVA|_cy',
                   '|CPD-12128|_cy', '|CPD-12125|_cy', '|CPD0-1065|_cy', '|CPD-15644|_cy',
                   '|UDP-MANNACA|_cy', '|CPD-15648|_cy', '|CPD-15647|_cy',
                   '|CPD-15646|_cy', '|CPD-15639|_cy', '|CPD-14795|_cy', '|CPD-9391|_cy',
                   '|CPD-9387|_cy', '|CPD-353|_cy', '|CPD-354|_cy', '|CPD-12121|_cy',
                   '|CPD-12127|_cy', '|CPD-12126|_cy', '|CPD-9067|_cy', '|SPERMIDINE|_cy',
                   '|CPD-13118|_cy', '|BIOTIN|_cy', '|Gro-P-Teichoic-peptidoglycan|_pe',
                   '|CPD-12129|_cy'
                   ]

artificial_DM = []
sinks.extend(artificial_EX)
sinks.extend(possible_source)
artificial_DM.extend(true_DM)
artificial_DM.extend(possible_product)
artificial_DM.extend(posible_biomass)
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
#   True dead ends of the model that need to be removed
# Repeat (4.1) until dead_ends is empty
print(len(dead_ends))
print(len(blocked_mets))
print(len(blocked_rxns))
# %% codecell
# (optional) Remove biomass components from dead_end metabolites
for met in dead_ends:
    if met in pc.biomass_in_model['met_id'].to_list():
        print(met)
        posible_biomass.extend([met])
        dead_ends.remove(met)

dead_ends
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
# (5) Check that all biomass components can be produced by the model. Do a first run
# with additional carbon sources open, then close all external reactions other than media.

media = ['|CH4|_ex', '|OXYGEN-MOLECULE|_ex', '|NITRATE|_ex', '|FE+2|_ex', '|Pi|_ex', '|SULFATE|_ex']

sinks = ['|WATER|_pe', '|WATER|_cy', '|Menaquinones|_it', '|ACP|_cy',
         '|CPD-1302|_cy', '|CPD-1301|_cy', '|Guanine34-in-tRNAs|_cy', '|CL-|_cy',
         '|CPD-302|_cy', '|Hpr-Histidine|_cy', '|biotin-L-lysine-in-BCCP-dimers|_cy',
         '|carboxybiotin-L-lysine-in-BCCP-dimers|_cy', '|S-ADENOSYL-4-METHYLTHIO-2-OXOBUTANOATE|_cy',
         '|Unsulfurated-Sulfur-Acceptors|_cy', '|Corrinoid-Adenosyltransferases|_cy',
         '|LysW-C-Terminal-L-Glutamate|_cy', '|CPD-17931|_pe', '|Cytochromes-C-Oxidized|_pe',
         '|CoI-Corrinoid-Fe-S-proteins|_cy', '|Cysteine-Desulfurase-L-cysteine|_cy',
         '|D-alanine-carrier-protein|_cy', '|DsrE3A-L-cysteine|_cy', '|Ox-Thioredoxin|_cy',
         '|Sulfur-Carrier-Proteins-ThiI|_cy', '|Thi-S|_cy'
         ]

artificial_EX = ['|NA+|_cy', '|MG+2|_cy', '|CO+2|_cy', '|Fatty-Acids|_cy',
                 '|CPD-12298|_cy', '|CPD-17989|_pe'
                 ]

true_DM = ['|NA+|_pe', '|DTDP-RHAMNOSE|_cy', '|ACP|_pe', '|FORMAMIDE|_cy',
           '|UNDECAPRENYL-DIPHOSPHATE|_pe', '|CPD-10640|_cy', '|UNKNOWN|_cy',
           '|N-ACETYL-D-GLUCOSAMINE|_cy', '|3-5-ADP|_cy', '|CPD-15999|_cy',
           '|ETHANOL-AMINE|_cy', '|CPD-1091|_cy', '|ALLYSINE|_cy', '|Alcohols|_cy',
           '|1-AMINO-PROPAN-2-OL|_cy', '|P3I|_cy', '|Cysteine-Desulfurase-L-cysteine|_cy',
           '|PPI|_cy', '|ADENOSYL-HOMO-CYS|_cy'
           ]

possible_product = ['|NITROGEN-MOLECULE|_pe', '|HYDROGEN-MOLECULE|_cy', '|ETOH|_cy',
                    '|CPD-10755|_cy', '|NITROGEN-MOLECULE|_cy', '|Methylketones|_cy',
                    '|CPD-347|_cy'
                    ]


artificial_DM = []
sinks.extend(artificial_EX)
artificial_DM.extend(true_DM)
artificial_DM.extend(possible_product)
demands = pc.biomass_in_model['met_id'].to_list()
# demands = ['|R-1-AMINOPROPAN-2-YL-PHOSPHATE|_cy']
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

    # Find a set of blocked metabolites
    blocked_mets, available_mets, produced_from_nothing = find_blocked_mets(
        model, demands = demands, carbon_source = '|CH4|_ex'
        )

# (5.1) Loop through blocked_mets and associated reactions, modifying constraints
# and stoichiometry or adding pathways until all biomass components can be produced from CH4
blocked_mets
produced_from_nothing
# %% codecell
# (5.2) It is possible to create new stoichiometrich balanced cycles with current
# configuration of sinks in the model, use find_leak_mode on produced_from_nothing
# metabolites to find which sinks and reactions are causing trouble

media = ['|OXYGEN-MOLECULE|_ex', '|NITRATE|_ex', '|FE+2|_ex', '|Pi|_ex', '|SULFATE|_ex']

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

    sparse_model, met_vars = build_sparse_model(
        model, keep_boundaries = True
        )
    leak_modes = find_leak_mode(
        sparse_model, leaks = produced_from_nothing, cutoff_mult = 10
        )

# %% codecell
leak = '|PYRUVATE|_cy'
mode = []
for r, flux in leak_modes[leak]:
    rxn = sparse_model.reactions.get_by_id(r)
    # mode.append([rxn.subsystem, rxn.id, rxn.build_reaction_string(), flux])
    mode.extend([rxn.id])

mode.sort()
for subs, id, r, flux in mode:
    print(subs, id, r, flux)

# %% codecell
w_model = sparse_model.copy()
rxns_to_remove = []
for rxn in w_model.reactions:
    if rxn.id in mode:
        pass
    else:
        rxns_to_remove.extend([rxn])
w_model.remove_reactions(rxns_to_remove)
consistent, inconsistent = pc.fastcc(w_model)
len(mode)
inconsistent
