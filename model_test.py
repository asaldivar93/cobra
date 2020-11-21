import cobra
from data_files.corrected_datasets import corrected_revesibility

# %%codecell
# for rxn in model.reactions:
#    if rxn.id not in multi_comp_rxns.keys() or rxn.id not in curated_rxns.keys() or rxn.id not in exchange_rxns.keys():
#        rxn.lower_bound = -100
#        rxn.upper_bound = 100

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

media = ['|CH4|_ou', '|OXYGEN-MOLECULE|_ou', '|NITRATE|_ou', '|FE+2|_ou',
         '|Pi|_ou', '|SULFATE|_ou', '|NA+|_ou', '|MG+2|_ou', '|CO+2|_ou',
         '|CL-|_ou']
met
with model as model:
    for substrate in media:
        met = model.metabolites.get_by_id(substrate)
        model.add_boundary(
            met, type = 'exchange'
            )

    medium = {rxn.id: 1000 for rxn in model.exchanges}
    medium['EX_|CH4|_ou'] = 100
    model.medium = medium

    # dm = model.add_boundary(
    #    model.metabolites.get_by_id('Intracellular_Metabolites'), type = 'demand'
    # )

    atp_dm = model.reactions.get_by_id('DM_|ATP|_cy')
    atp_dm.lower_bound = 3.5
    atp_dm.upper_bound = 1000

    o2_ex = model.reactions.get_by_id('O2')
    o2_ex.upper_bound = 1000
    o2_ex.lower_bound = 0
    o2_ex = model.reactions.get_by_id('EX_O2')
    o2_ex.upper_bound = 1000
    o2_ex.lower_bound = 0
    dm = model.metabolites.get_by_id('biomass')
    model.add_boundary(
        dm, type = 'demand'
        )

    for r in ['|ATP|_cy_syn', '|1.10.2.2-RXN|', '|CYTOCHROME-C-OXIDASE-RXN|', '|NADH-DEHYDROG-A-RXN|', 'DM_biomass']:
        rxn = model.reactions.get_by_id(r)
        rxn.upper_bound = 1000

    model.objective = {
        # model.reactions.get_by_id('DM_Intracellular_Metabolites'): 1.0
        model.reactions.get_by_id('BIOMASS'): 1.0
        }

    sol = cobra.flux_analysis.pfba(model)
    print(model.summary())
    print(model.metabolites.get_by_id('|OXYGEN-MOLECULE|_cy').summary())
