import pandas as pd
import p_tools
import os
import cobra
import warnings
warnings.filterwarnings("ignore")
from tqdm import trange
from data_files.corrected_datasets import (multi_comp_rxns,
                                           curated_rxns,
                                           exchange_rxns,
                                           added_if_met,
                                           added_pthwys,
                                           true_DM,
                                           possible_product,
                                           corrected_revesibility,
                                           sinks)
from cobra.flux_analysis import (gapfill,
                                 find_blocked_mets)
from cobra_tools import (search_biomass_components,
                         add_boundaries_for_simulation)
import picrust_parser as pparser

pc = pparser.p_model()
pd.set_option('display.max_rows', 150)
artificial_DM = []
artificial_DM.extend(true_DM)

# %% codecell
pt = p_tools.picrust_tools()
pwy_strat = pt.sample_strat_pathway()
taxonomy = pd.read_csv('/home/alexis/UAM/cobra/data_files/silva_taxonomy.csv', sep ='\t')
taxonomy.set_index('sequence_identifier', inplace = True)
seq_counts = pd.read_csv('/home/alexis/UAM/cobra/data_files/otu-table-corrected.tsv', sep = '\t', index_col = 'OTU_ID')
taxonomy.fillna('unknown', inplace = True)
samples = ['My_20', 'CIR_19']
otus = pd.DataFrame()
for sample in samples:
    otus = otus.join(seq_counts[sample].loc[seq_counts[sample] != 0].to_frame(), how = 'outer')

for seq in otus.index:
    otus.loc[seq, 'Family'] = taxonomy.loc[seq, 'Family']
    otus.loc[seq, 'Genus'] = taxonomy.loc[seq, 'Genus']


# %% codecel
pathways = pd.DataFrame()
for sample in samples:
    pathways = pathways.join(pwy_strat[sample], how = 'outer', rsuffix='r')

pwys_in_gens = {}
for genus in otus['Genus'].unique():
    if genus == 'uncultured' or genus == 'unknown':
        gen_filter = otus.loc[otus['Genus'] == genus].index
        families = otus.loc[gen_filter, :]
        for family in families['Family'].unique():
            fam_filter = otus.loc[otus['Family'] == family].index
            gen_pwy = pc.parse_pathways_in_genus(pathways, fam_filter)
            pwys_in_gens[family] = gen_pwy
    else:
        gen_filter = otus.loc[otus['Genus'] == genus].index
        gen_pwy = pc.parse_pathways_in_genus(pathways, gen_filter)
        pwys_in_gens[genus] = gen_pwy

for genus in pwys_in_gens.keys():
    print(genus, len(pwys_in_gens[genus]))

# %% codecell
os.chdir('/home/alexis/UAM/cobra/')

def build_genus_model(genus):
    model = cobra.Model()
    pc = pparser.p_model()
    pathways = pwys_in_gens[genus].copy()
    pathways.extend(added_pthwys)
    gen_id = genus[:2] + genus[-2:]
    model.name = genus
    for i in trange(int(len(pathways))):
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

    for i in trange(int(len(pwys_to_extend))):
        pthwy = pwys_to_extend[i]
        pthwy_id = '|' + pthwy + '|'
        if pthwy_id in pparser.all_pathways:
            pc.add_pathway(model, pthwy_id)
        else:
            pc.unmatched_pthwys.extend([pthwy_id])

    for rxn_id in curated_rxns.keys():
        if rxn_id in model.reactions:
            model.remove_reactions(
                model.reactions.get_by_id(rxn_id)
            )
        subsystem = curated_rxns[rxn_id]['pathway']
        stoichiometry = curated_rxns[rxn_id]['stoichiometry']
        reversible = curated_rxns[rxn_id]['reversible']
        direction = curated_rxns[rxn_id]['direction']
        pc.add_rxn_from_stoichiometry(model, subsystem, rxn_id, stoichiometry, reversible, direction)

    for rxn_id in multi_comp_rxns.keys():
        if rxn_id in model.reactions:
            model.remove_reactions(
                model.reactions.get_by_id(rxn_id)
            )
        subsystem = multi_comp_rxns[rxn_id]['pathway']
        stoichiometry = multi_comp_rxns[rxn_id]['stoichiometry']
        reversible = multi_comp_rxns[rxn_id]['reversible']
        direction = multi_comp_rxns[rxn_id]['direction']
        pc.add_rxn_from_stoichiometry(model, subsystem, rxn_id, stoichiometry, reversible, direction)
        pc.multi_compartment_rxns.append([subsystem, rxn_id, 'added'])

    for ex_rxn in exchange_rxns.keys():
        if ex_rxn in model.reactions:
            model.remove_reactions(
                model.reactions.get_by_id(ex_rxn)
            )
        subsystem = exchange_rxns[ex_rxn]['pathway']
        stoichiometry = exchange_rxns[ex_rxn]['stoichiometry']
        reversible = exchange_rxns[ex_rxn]['reversible']
        direction = exchange_rxns[ex_rxn]['direction']
        pc.add_rxn_from_stoichiometry(
            model, subsystem, ex_rxn, stoichiometry, reversible, direction
            )

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

    remove_mets = []
    for met in model.metabolites:
        if not met.reactions:
            remove_mets.extend([met])
    model.remove_metabolites(remove_mets)
    biomass_in_model = pc.search_biomass_components(model)
    model = pc.add_biomass_rxn(model, biomass_in_model)

    print(
        genus,
        model.metabolites.get_by_id('biomass').formula
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

    for met in model.metabolites:
        if '_ou' in met.id or '_ex' in met.id:
            met.id = '_' + met.id.replace('|', '').replace('.', '').replace('-', '_').replace('+', '_')
        else:
            met.id = '_' + met.id.replace('|', '').replace('.', '').replace('-', '_').replace('+', '_')
        model.repair()
    for met in model.metabolites:
        if not met.name:
            met.name = met.id
        if not met.compartment:
            met.compartment = 'cytosol'
        met.charge = 0
        met.annotation = {}

    for rxn in model.reactions:
        if 'ou_' in rxn.id:
            rxn.id = '_' + rxn.id.replace('|', '').replace('.', '').replace('-', '_').replace('+', '_')
        else:
            rxn.id = '_' + rxn.id.replace('|', '').replace('.', '').replace('-', '_').replace('+', '_')
        model.repair()
    for rxn in model.reactions:
        rxn.annotation = {}

    cobra.io.write_sbml_model(model, 'models/models/{}.xml'.format(gen_id))


# %% codecell
for genus in pwys_in_gens.keys():
    build_genus_model(genus)

# %% codecell
model_list = os.listdir('models/models/')
added_to_genus = {}
need_manual_check = []
for k in model_list:
    media = ['_OXYGEN_MOLECULE_ou', '_NITRATE_ou', '_FE_2_ou',
             '_Pi_ou', '_SULFATE_ou', '_NA__ou', '_MG_2_ou', '_CO_2_ou',
             '_CL__ou']
    universal = cobra.io.read_sbml_model('models/cir.xml')
    universal.solver = 'cplex'
    genus_model = cobra.io.read_sbml_model('models/models/{}'.format(k))
    genus_model.solver = 'cplex'
    genus = genus_model.name
    if genus == 'Methylocystis':
        print(genus)
        carbon_source = ['_CH4_ou']
        sinks = []
        exchange_reactions = False
        for rxn in genus_model.reactions:
            if '_RXN_2961' in rxn.id:
                rxn.lower_bound = -100
                rxn.upper_bound = 100
            if '_PYRUFLAVREDUCT_RXN' in rxn.id:
                rxn.lower_bound = -100
                rxn.upper_bound = 0
            if '_CITSYN_RXN' in rxn.id:
                rxn.lower_bound = 0
                rxn.upper_bound = 100
            if '_2OXOGLUTARATEDEH_RXN' in rxn.id:
                rxn.lower_bound = 0
                rxn.upper_bound = 100
            if '_BUTYRYL_COA_DEHYDROGENASE_RXN' in rxn.id:
                rxn.lower_bound = 0
                rxn.upper_bound = 0
            if '_PROPIONATE__COA_LIGASE_RXN' in rxn.id:
                rxn.lower_bound = 0
                rxn.upper_bound = 0
            if '_METHYLENETHFDEHYDROG_NADP_RXN_Meis' in rxn.id:
                rxn.lower_bound = 0
                rxn.upper_bound = 0
            if '_2131_RXN_Meis' in rxn.id:
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
                rxn.lower_bound = 0
                rxn.upper_bound = 0
            if '_RXN_7566' in rxn.id:
                rxn.lower_bound = 0
                rxn.upper_bound = 0
            if '_ASPARTASE_RXN' in rxn.id:
                rxn.lower_bound = -100
                rxn.upper_bound = 0
            if '_GLUTAMATE_DEHYDROGENASE_RXN' in rxn.id:
                rxn.lower_bound = 0
                rxn.upper_bound = 0
            if '_KETOGLUTREDUCT_RXN' in rxn.id:
                rxn.lower_bound = 0
                rxn.upper_bound = 0
            if '_RXN_11662' in rxn.id:
                rxn.lower_bound = 0
                rxn.upper_bound = 0
            if '_2_METHYLCITRATE_SYNTHASE_RXN' in rxn.id:
                rxn.lower_bound = 0
                rxn.upper_bound = 0
    elif genus == 'Methylophilus' or genus == 'Hyphomicrobium':
        sinks = []
        print(genus)
        carbon_source = ['_METOH_ou']
        exchange_reactions = False
        for rxn in genus_model.reactions:
            if '_RXN_12219' in rxn.id:
                rxn.lower_bound = 0
                rxn.upper_bound = 0
            if '_METHANE_MONOOXYGENASE_RXN' in rxn.id:
                rxn.lower_bound = 0
                rxn.upper_bound = 0
            if '_EX_METOH' in rxn.id:
                rxn.lower_bound = 0
                rxn.upper_bound = 100
            if '_EX_METOH_1' in rxn.id:
                rxn.lower_bound = 0
                rxn.upper_bound = 100
        for rxn in universal.reactions:
            if '_RXN_12219' in rxn.id:
                rxn.lower_bound = 0
                rxn.upper_bound = 0
            if '_METHANE_MONOOXYGENASE_RXN' in rxn.id:
                rxn.lower_bound = 0
                rxn.upper_bound = 0
            if '_EX_METOH' in rxn.id:
                rxn.lower_bound = 0
                rxn.upper_bound = 100
            if '_EX_METOH_1' in rxn.id:
                rxn.lower_bound = 0
                rxn.upper_bound = 100
    else:
        print(genus)
        sinks = ['_ETOH_cy', '_CPD_10755_cy', '_Methylketones_cy', '_CPD_347_cy',
                 '_PROPIONATE_cy', '_PROPANE_1_2_DIOL_cy', '_CPD_10353_cy',
                 '_BUTANEDIOL_cy', '_BUTANOL_cy', '_ACETONE_cy', '_FORMATE_cy',
                 '_ACET_cy', '_BUTYRIC_ACID_cy', '_PUTRESCINE_cy']
        carbon_source = ['_CH4_ou']
        for rxn in genus_model.reactions:
            if '_RXN_12219' in rxn.id:
                rxn.lower_bound = 0
                rxn.upper_bound = 0
            if '_METHANE_MONOOXYGENASE_RXN' in rxn.id:
                rxn.lower_bound = 0
                rxn.upper_bound = 0
            if 'EX_METOH' in rxn.id:
                rxn.lower_bound = 0
                rxn.upper_bound = 0
            if '_EX_METOH_1' in rxn.id:
                rxn.lower_bound = 0
                rxn.upper_bound = 0
        for rxn in universal.reactions:
            if '_RXN_12219' in rxn.id:
                rxn.lower_bound = 0
                rxn.upper_bound = 0
            if '_METHANE_MONOOXYGENASE_RXN' in rxn.id:
                rxn.lower_bound = 0
                rxn.upper_bound = 0
            if '_EX_METOH' in rxn.id:
                rxn.lower_bound = 0
                rxn.upper_bound = 0
            if '_EX_METOH_1' in rxn.id:
                rxn.lower_bound = 0
                rxn.upper_bound = 0

    media.extend(carbon_source)
    bm_in_model = search_biomass_components(genus_model)
    demands = bm_in_model['met_id'].to_list()

    already_added = set()
    to_add = set()
    with genus_model as g_model:
        add_boundaries_for_simulation(g_model, media, sinks = sinks)
        blocked_mets, available_mets, produced_from_nothing = find_blocked_mets(
            g_model, demands = demands, carbon_source = carbon_source[0]
            )
        for rxn in g_model.reactions:
            if 'EX_METOH' in rxn.id:
                print(rxn.id, rxn.lower_bound, rxn.upper_bound)
        with universal as u_model:
            iter = 0
            while blocked_mets and iter < len(demands):
                print(blocked_mets)
                for blocked in blocked_mets:
                    g_model.objective = g_model.add_boundary(g_model.metabolites.get_by_id(blocked), type = 'demand')
                    added_from_gapfill = gapfill(g_model, u_model, demand_reactions=True)
                    added_from_gapfill = set(rxn.id for rxn in added_from_gapfill[0])
                    to_add = added_from_gapfill.difference(already_added)
                    already_added.update(added_from_gapfill)
                    for rxn in to_add:
                        if 'DM_' in rxn:
                            metabolite = rxn.replace('DM_', '')
                            if metabolite not in g_model.metabolites:
                                g_model.add_metabolites(
                                    u_model.metabolites.get_by_id(metabolite).copy()
                                )
                            g_model.add_boundary(
                                g_model.metabolites.get_by_id(metabolite),
                                type = 'demand'
                                )
                        else:
                            reaction = u_model.reactions.get_by_id(rxn)
                            for met in reaction.metabolites:
                                if met not in g_model.metabolites:
                                    g_model.add_metabolites(met.copy())
                            g_model.add_reaction(reaction)
                blocked_mets, available_mets, produced_from_nothing = find_blocked_mets(
                    g_model, demands = demands, carbon_source = carbon_source[0]
                    )
                iter += 1
        g_model.objective = g_model.reactions._BIOMASS
        biomass_flux = g_model.slim_optimize()

    if biomass_flux > 1:
        for rxn in already_added:
            if 'DM_' in rxn:
                metabolite = rxn.replace('DM_', '')
                if metabolite not in genus_model.metabolites:
                    genus_model.add_metabolites(
                        universal.metabolites.get_by_id(metabolite).copy()
                    )
                genus_model.add_boundary(
                    genus_model.metabolites.get_by_id(rxn.replace('DM_', '')),
                    type = 'demand'
                )
            else:
                reaction = universal.reactions.get_by_id(rxn)
                for met in reaction.metabolites:
                    if met not in genus_model.metabolites:
                        genus_model.add_metabolites(met.copy())
                genus_model.add_reaction(reaction)
        added_to_genus[genus] = already_added.copy()
        for met in genus_model.metabolites:
            if '_ou' in met.id or '_ex' in met.id:
                pass
            else:
                met.id = met.id + '_' + k[:4]
            genus_model.repair()
        for rxn in genus_model.reactions:
            if 'ou_' in rxn.id or ('DM_' and '_ex') in rxn.id:
                pass
            else:
                rxn.id = rxn.id + '_' + k[:4]
            genus_model.repair()
        cobra.io.write_sbml_model(genus_model, 'models/models_gapfilled/{}'.format(k))
    else:
        need_manual_check.append(genus)

# %%codecell
media = ['_CH4_ou', '_OXYGEN_MOLECULE_ou', '_NITRATE_ou', '_FE_2_ou',
         '_Pi_ou', '_SULFATE_ou', '_NA__ou', '_MG_2_ou', '_CO_2_ou',
         '_CL__ou']
Meis = cobra.io.read_sbml_model('models/Meis.xml')
k_model = cobra.io.read_sbml_model('models/CIR_19_complete/Meus.xml')
k_model.merge(Meis)
universal = cobra.io.read_sbml_model('models/universal/Meus.xml')
u_meis = cobra.io.read_sbml_model('models/universal/Meis.xml')
universal.merge(u_meis)
with k_model as model:
    add_boundaries_for_simulation(model, media)
    solution = gapfill(model, universal, community = True, taxa_database = taxa_database, demand_reactions = True, exchange_reactions = True)

solution

universal = cobra.io.read_sbml_model('models/universal/Meus.xml')

for rxn in universal.reactions:
    if 'EX' in rxn.id:
        print(rxn.id)
