import cobra
from tools.cobra_tools import add_boundaries_for_simulation, search_biomass_components

# %%codecell
cir = cobra.io.read_sbml_model('models/models_gapfilled/Meis.xml')
cir = cobra.io.read_sbml_model('models/cir.xml')

# %% codecell
# Reactions with corrected reversibility for Meis
for reaction in cir.reactions:
    if '_RXN_2961' in reaction.id:
        reaction.lower_bound = 0
        reaction.upper_bound = 0
    if '_PYRUFLAVREDUCT_RXN' in reaction.id:
        reaction.lower_bound = -100
        reaction.upper_bound = 0
    if '_CITSYN_RXN' in reaction.id:
        reaction.lower_bound = 0
        reaction.upper_bound = 100
    if '_2OXOGLUTARATEDEH_RXN' in reaction.id:
        reaction.lower_bound = 0
        reaction.upper_bound = 100
    if '_BUTYRYL_COA_DEHYDROGENASE_RXN' in reaction.id:
        reaction.lower_bound = 0
        reaction.upper_bound = 0
    if '_PROPIONATE__COA_LIGASE_RXN' in reaction.id:
        reaction.lower_bound = 0
        reaction.upper_bound = 0
    if '_METHYLENETHFDEHYDROG_NADP_RXN_Meis' in reaction.id:
        reaction.lower_bound = 0
        reaction.upper_bound = 0
    if '_2131_RXN_Meis' in reaction.id:
        reaction.lower_bound = 0
        reaction.upper_bound = 0
    if '_RXN0_268_Meis' in reaction.id:
        reaction.lower_bound = 0
        reaction.upper_bound = 0
    if '_PEPSYNTH_RXN_Mei' in reaction.id:
        reaction.lower_bound = 0
        reaction.upper_bound = 0

# %% codecell
# Reactions with corrected reversibility for CIR
for rxn in cir.reactions:
    if '_RXN_2961' in rxn.id:
        print(rxn.id, rxn.lower_bound, rxn.upper_bound)
        rxn.lower_bound = 0
        rxn.upper_bound = 0
    if 'RXN_2802' in rxn.id:
        print(rxn.id, rxn.lower_bound, rxn.upper_bound)
        rxn.lower_bound = 0
        rxn.upper_bound = 0
    if 'R10_RXN' in rxn.id:
        print(rxn.id, rxn.lower_bound, rxn.upper_bound)
        rxn.lower_bound = 0
        rxn.upper_bound = 0
    if 'ISOCIT_CLEAV_RXN' in rxn.id:
        print(rxn.id, rxn.lower_bound, rxn.upper_bound)
        rxn.lower_bound = 0
        rxn.upper_bound = 0
    if 'RXN_11489' in rxn.id:
        print(rxn.id, rxn.lower_bound, rxn.upper_bound)
        rxn.lower_bound = -100
        rxn.upper_bound = 0
    if '_GLY3KIN_RXN' in rxn.id:
        print(rxn.id, rxn.lower_bound, rxn.upper_bound)
        rxn.lower_bound = 0
        rxn.upper_bound = 0
    if 'CITSYN_RXN' in rxn.id:
        print(rxn.id, rxn.lower_bound, rxn.upper_bound)
        rxn.lower_bound = 0
        rxn.upper_bound = 100
    if '_RXN_20917' in rxn.id:
        print(rxn.id, rxn.lower_bound, rxn.upper_bound)
        rxn.lower_bound = 0
        rxn.upper_bound = 0
    if '_EX_METOH_1' in rxn.id:
        rxn.lower_bound = -100
        rxn.upper_bound = 0
    if '_GLUTAMATE_DEHYDROGENASE_NADP__RXN' in rxn.id:
        rxn.lower_bound = 0
        rxn.upper_bound = 100
    if '_2_1_3_1_RXN' in rxn.id:
        rxn.lower_bound = 0
        rxn.upper_bound = 0
    if '_RXN_12168' in rxn.id:
        rxn.lower_bound = 0
        rxn.upper_bound = 0
    if '_R125_RXN' in rxn.id:
        rxn.lower_bound = 0
        rxn.upper_bound = 0
    if '_GLUTACONYL_COA_DECARBOXYLASE_RXN' in rxn.id:
        rxn.lower_bound = 0
        rxn.upper_bound = 0
    if '_PUTTRANSAM_RXN' in rxn.id:
        rxn.lower_bound = 0
        rxn.upper_bound = 0
    if '_RXN_15125' in rxn.id:
        rxn.lower_bound = 0
        rxn.upper_bound = 0
    if '_RXNI_3' in rxn.id:
        rxn.lower_bound = 0
        rxn.upper_bound = 0
    if '_R601_RXN' in rxn.id:
        rxn.lower_bound = 0
        rxn.upper_bound = 0
    if '_RXN_7774' in rxn.id:
        rxn.lower_bound = 0
        rxn.upper_bound = 0
    if '_RXNI_2' in rxn.id:
        rxn.lower_bound = 0
        rxn.upper_bound = 0
    if '_GCVMULTI_RXN' in rxn.id:
        rxn.lower_bound = -100
        rxn.upper_bound = 100
    if '_RXN_7566' in rxn.id:
        rxn.lower_bound = 0
        rxn.upper_bound = 0
    if '_ASPARTASE_RXN' in rxn.id:
        rxn.lower_bound = -100
        rxn.upper_bound = 0
    if '_GLUTAMATE_DEHYDROGENASE_RXN' in rxn.id:
        rxn.lower_bound = -100
        rxn.upper_bound = 100
    if '_KETOGLUTREDUCT_RXN' in rxn.id:
        rxn.lower_bound = 0
        rxn.upper_bound = 0
    if '_RXN_11662' in rxn.id:
        rxn.lower_bound = 0
        rxn.upper_bound = 0
    if '_2_METHYLCITRATE_SYNTHASE_RXN' in rxn.id:
        rxn.lower_bound = 0
        rxn.upper_bound = 0
    if '_PROPKIN_RXN' in rxn.id:
        rxn.lower_bound = 0
        rxn.upper_bound = 0
    if '_3_HYDROXYBUTYRATE_DEHYDROGENASE_RXN' in rxn.id:
        rxn.lower_bound = -100
        rxn.upper_bound = 100
    if '_ACETOACETYL_COA_TRANSFER_RXN' in rxn.id:
        rxn.lower_bound = -100
        rxn.upper_bound = 100
    if '_RXN_8807' in rxn.id:
        rxn.lower_bound = -100
        rxn.upper_bound = 100
    if '_ACETATEKIN_RXN' in rxn.id:
        rxn.lower_bound = 0
        rxn.upper_bound = 0
    if '_ACETATE__COA_LIGASE_ADP_FORMING_RXN' in rxn.id:
        rxn.lower_bound = 0
        rxn.upper_bound = 0

    if '_ASPAMINOTRANS_RXN' in rxn.id:
        rxn.lower_bound = 0
        rxn.upper_bound = 0
    if '_GLYCINE_AMINOTRANSFERASE_RXN' in rxn.id:
        rxn.lower_bound = 0
        rxn.upper_bound = 0






# %% codecell
media = ['_CH4_ou', '_OXYGEN_MOLECULE_ou', '_NITRATE_ou', '_FE_2_ou',
         '_Pi_ou', '_SULFATE_ou', '_NA__ou', '_MG_2_ou', '_CO_2_ou',
         '_CL__ou']

with cir as model:
    add_boundaries_for_simulation(model, 0.00625, media)
    model.objective = {rxn: 1 for rxn in model.reactions if 'BIOMASS' in rxn.id}
    # for metabolite in model.metabolites:
    #     if '_PUTRESCINE_cy' in metabolite.id:
    #         obj = model.add_boundary(
    #             model.metabolites._PUTRESCINE_cy_Meis,
    #             type = 'demand'
    #         )
    sol = model.optimize()
    print(model.summary(solution = sol))
    print(model.metabolites.get_by_id('_METHYLENE_THF_cy').summary(solution = sol))
4.115E-05 + 0.0023907 -1.779E-05-2.908E-07-5.784E-06
0.001402+0.0002437
for rxn in cir.reactions:
    if cir.metabolites.get_by_id('_ACETYL_COA_cy') in rxn.metabolites.keys():
        if cir.metabolites.get_by_id('_PYRUVATE_cy') in rxn.metabolites.keys():
            print(rxn.id, rxn.lower_bound, rxn.upper_bound, rxn.build_reaction_string())
cir.reactions.get_by_id('_GCVMULTI_RXN').build_reaction_string()
for met in cir.metabolites:
    if '_ACET_cy' in met.id:
        for rxn in met.reactions:
            print(rxn.id,rxn.lower_bound, rxn.upper_bound, rxn.build_reaction_string())
# %%codecell
0.0006309*2
a = pd.DataFrame(
    [['O<sub>2</sub>', round(0.009458 / 0.00625, 3)], ['CO<sub>2</sub>', round(0.00424 / 0.00625,3)], ['Biomasa', round(0.05421 * 0.03772482056203139 / 0.00625, 3)]],
    columns=['Metabolite', 'Rendimiento (C-mol C-mol<sup>-1</sup>)']
)
fig = px.bar(
    a,
    x='Metabolite',
    y='Rendimiento (C-mol C-mol<sup>-1</sup>)',
    template='none',
    text='Rendimiento (C-mol C-mol<sup>-1</sup>)'
    )
fig.show()
fig.write_image(results_path + 'plots/yields_general.svg')

model.metabolites._biomass.cmol_by_mol
# %% codecell

bm = search_biomass_components(cir)
class_filter = bm['class'].isin(
    ['Intracellular_Metabolites', 'Carbohydrates']
)

to_exchange = bm.loc[class_filter, 'met_id'].to_list()

for metabolite in to_exchange:
    met_to_ex = metabolite[:-5].replace('cy', 'ex')
    print(metabolite, met_to_ex)
    stoichiometry = {}
    stoichiometry[met_to_ex] = -1
    stoichiometry[metabolite] = 1

    met_to_add = cobra.Metabolite(
        met_to_ex,
        name = met_to_ex[:-3],
        compartment = 'extracellular'
        )

    rxn_id = '_EX' + met_to_ex[:-3]
    rxn_to_add = cobra.Reaction(
        rxn_id,
        upper_bound = 1000,
        lower_bound = -1000
        )

    cir.add_metabolites(
        met_to_add
    )

    cir.add_reaction(
        rxn_to_add
    )

    rxn_to_add.add_metabolites(
        stoichiometry
    )

    cir.add_boundary(
        met_to_add,
        type = 'demand'
    )

for rxn in cir.reactions:
    if 'EX' in rxn.id:
        print(rxn.id)

cobra.io.write_sbml_model(cir, 'model/Meis.xml')
