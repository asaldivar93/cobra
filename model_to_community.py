import pandas as pd
import p_tools
from tqdm import tqdm
import cobra
from data_files.corrected_datasets import (multi_comp_rxns,
                                           curated_rxns,
                                           exchange_rxns,
                                           added_if_met,
                                           added_pthwys,
                                           sinks,
                                           artificial_EX,
                                           true_DM,
                                           possible_product)
import picrust_parser as pparser
balanced_rxns = pd.read_csv('data_files/balanced_rxns.csv', header=None, dtype='str')
balanced_rxns = balanced_rxns.iloc[:, 0].to_list()
ignore_massbalance = pd.read_csv('data_files/ignore_massbalance.csv', header=None, dtype='str')
ignore_massbalance = ignore_massbalance.iloc[:, 0].to_list()

artificial_DM = []
sinks.extend(artificial_EX)
artificial_DM.extend(true_DM)
# %% codecell
pt = p_tools.picrust_tools()
pwy_strat = pt.sample_strat_pathway()

# %% codecell
taxonomy = pd.read_csv('/home/alexis/UAM/cobra/data_files/silva_taxonomy.csv', sep ='\t')
taxonomy.set_index('sequence_identifier', inplace = True)
seq_counts = pd.read_csv('/home/alexis/UAM/cobra/data_files/otu-table-corrected.tsv', sep = '\t', index_col = 'OTU_ID')
taxonomy.fillna('unknown', inplace = True)
otus = seq_counts['CIR_19'].loc[seq_counts['CIR_19'] != 0].to_frame()

for seq in otus.index:
    otus.loc[seq, 'Family'] = taxonomy.loc[seq, 'Family']
    otus.loc[seq, 'Genus'] = taxonomy.loc[seq, 'Genus']


# %% codecell
pwys_in_gens = {}
for gen in otus['Genus'].unique():
    gen_filter = otus.loc[otus['Genus'] == gen].index
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
    pwys_in_gens[gen] = gen_pwy

# %% codecell
pwys_in_gens[gen]
for gen in pwys_in_gens.keys():
    print(gen, len(pwys_in_gens[gen]))

# %% codecell
model = cobra.Model()
pc = pparser.p_model()
pathways = pwys_in_gens['Methylocystis'].copy()
pathways.append(added_pthwys)

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

for i in tqdm(range(int(len(pathways)))):
    pthwy = pathways[i]
    pthwy_id = '|' + pthwy + '|'
    if pthwy_id in pparser.all_pathways:
        for rxn_id in pparser.metacyc_db[pthwy_id].reaction_list:
            substrates = pparser.metacyc_db.reaction_reactants_and_products(rxn_id, pwy=pthwy_id)
            if not all(substrates):
                pass
            else:
                # Create a list of compartments of reaction
                comps_of_rxn = []
                for compartment in pparser.metacyc_db.compartments_of_reaction(rxn_id):
                    compartment = compartment.replace('|', '').replace('CCO-', '').lower()
                    if compartment == 'peri-bac' or compartment == 'out' or compartment == 'pm-bac-neg':
                        compartment = 'periplasm'
                    elif compartment == 'in':
                        compartment = 'cytosol'
                    comps_of_rxn.extend([compartment])
                    if compartment not in pc.compartments:
                        pc.compartments.extend([compartment])

                # Multicompartment and transport rxns need to be parsed manually
                # If rxn is multicompartment
                if len(comps_of_rxn) > 1:
                    if rxn_id not in multi_comp_rxns.keys():
                        pc.multi_compartment_rxns.append([pthwy_id, rxn_id])
                else:
                    # Create a Dictionary of reactants and products
                    substrates = pparser.metacyc_db.reaction_reactants_and_products(rxn_id, pwy=pthwy_id)
                    if not substrates:
                        substrates = [pparser.metacyc_db[rxn_id].left,
                                      pparser.metacyc_db[rxn_id].right]
                    mets_in_reaction = dict(zip(['reactants', 'products'],
                                                substrates))

                    # Create a List of metabolites in the reaction
                    metabolites = mets_in_reaction['reactants'].copy()
                    metabolites.extend(mets_in_reaction['products'])
                    # Replace metabolites that fall under the redox_pairs class
                    metabolites = pc.replace_redox_pairs(metabolites)

                    comp = '_' + comps_of_rxn[0][0:2]
                    # For every metabolite in reaction
                    for met in metabolites:
                        # Create id with corresponding compartment
                        met_id = met + comp
                        # If the metabolite isn't in the model already
                        if met_id not in model.metabolites:
                            # Parse metabolite info
                            met_to_add = pc.parse_metabolite(met, comp)
                            met_to_add.compartment = comps_of_rxn[0]
                            # add metabolite
                            model.add_metabolites(met_to_add)
                    # Update metabolites list with the corresponding compartment
                    metabolites = [met + comp for met in metabolites]

                    # if reaction already in model
                    if rxn_id in model.reactions:
                        # add_ subsystem to reaction
                        rxn = model.reactions.get_by_id(rxn_id)
                        if pthwy_id not in rxn.subsystem:
                            rxn.subsystem = rxn.subsystem + ' ' + pthwy_id
                    else:
                        # Add reaction
                        rxn_to_add = cobra.Reaction(rxn_id)
                        model.add_reaction(rxn_to_add)

                        # Add Stoichiometry
                        stoichiometry = pc.parse_stoichiometry(model, rxn_id, mets_in_reaction, metabolites)
                        rxn_to_add.add_metabolites(stoichiometry)

                        # Contra-intuitive check_mass_balance returns false when reaction is mass_balanced
                        # If rxn is in list of known balanced rxns
                        if rxn_id in balanced_rxns or rxn_id in ignore_massbalance:
                            # then mass_balanced = True
                            mass_balanced = {}
                        else:
                            # else check mass balace of rxn
                            mass_balanced = rxn_to_add.check_mass_balance()

                        # update Stoichiometry
                        # If reaction is unbalanced by equal H and charge
                        unbalanced_elements = list(mass_balanced.keys())
                        if len(unbalanced_elements) == 2 and ('charge' in unbalanced_elements and 'H' in unbalanced_elements):
                            if mass_balanced['charge'] == mass_balanced['H']:
                                # and if H+ is one of the metabolites
                                if '|PROTON|' + comp in stoichiometry.keys():
                                    # Then Rxn can be balanced by the addition of H+
                                    rxn_to_add.add_metabolites({'|PROTON|' + comp: -mass_balanced['H']})
                                    mass_balanced = rxn_to_add.check_mass_balance()
                                    pc.corrected_hc_rxns.extend([rxn_id])

                        # If rxn is mass_balanced
                        if not mass_balanced:
                            lower_bound, upper_bound, dG_lb, dG_ub, ri_lb, ri_ub = pc.add_thermo_constraints(rxn_to_add)
                        # if rxn is not mass balanced
                        else:
                            pc.classify_unbalanced_rxn(rxn_to_add)
                            # set rxn as reversible
                            dG_lb = []
                            dG_ub = []
                            ri_lb = []
                            ri_ub = []
                            lower_bound = -100
                            upper_bound = 100

                        rxn_to_add.lower_bound = lower_bound
                        rxn_to_add.upper_bound = upper_bound
                        genes_of_reaction = [gene.replace('|', '')
                                             for gene in pparser.metacyc_db.genes_of_reaction(rxn_id)]
                        if genes_of_reaction:
                            rxn_to_add.gene_reaction_rule = '( ' + ' or '.join(genes_of_reaction) + ' )'

                        rxn_to_add.subsystem = pthwy_id
                        atom_mappings = pparser.metacyc_db[rxn_id].atom_mappings
                        ec_number = pparser.metacyc_db[rxn_id].ec_number
                        dblinks = pparser.metacyc_db[rxn_id].dblinks
                        rxn_to_add.annotation = {'atom_mappings': atom_mappings,
                                                 'ec_number': ec_number,
                                                 'dblinks': dblinks,
                                                 'dG': [dG_lb, dG_ub],
                                                 'rIndex': [ri_lb, ri_ub]}
    else:
        pc.unmatched_pthwys.extend([pthwy_id])

    add_pwys = []
    for b_met in added_if_met['metabolite']:
        for met in model.metabolites:
            if b_met in met.id:
                add_pwys.append([True])
            else:
                add_pwys.append([False])

    for path in added_if_met.loc[add_pwys, 'pathway']:
        if path not in pathways:
            pathways.append(path)

biomass_in_model = pc.search_biomass_components(model)
model = pc.add_biomass_rxn(model, biomass_in_model)

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

for r in ['|ATP|_cy_syn', '|1.10.2.2-RXN|', '|CYTOCHROME-C-OXIDASE-RXN|', '|NADH-DEHYDROG-A-RXN|', 'DM_biomass']:
    rxn = model.reactions.get_by_id(r)
    rxn.upper_bound = 1000
