
# missing method to claculate dg of transport reactions and multicompartment reactions
# %% codecell
import pandas as pd
import os
import cobra
os.chdir('/home/alexis/UAM/cobra')
import picrust_parser as pparser
balanced_rxns = pd.read_csv('data_files/balanced_rxns.csv', header=None, dtype='str')
balanced_rxns = balanced_rxns.iloc[:, 0].to_list()
ignore_massbalance = pd.read_csv('data_files/ignore_massbalance.csv', header=None, dtype='str')
ignore_massbalance = ignore_massbalance.iloc[:, 0].to_list()

# %% codecell
model = cobra.Model()
pc = pparser.p_model()
n = 0
step = 0.1
for pthwy in pc.pathways:
    n += 1
    if n / len(pc.pathways) > step:
        step += 0.1
        print(str(n) + ' pathways parsed %3.2f completed' % (n / len(pc.pathways) * 100))

    pthwy_id = '|' + pthwy + '|'
    if pthwy_id in pparser.all_pathways:
        for rxn_id in pparser.metacyc_db[pthwy_id].reaction_list:
            # Create a list of compartments of reaction
            comps_of_rxn = []
            for compartment in pparser.metacyc_db.compartments_of_reaction(rxn_id):
                compartment = compartment.replace('|', '').replace('CCO-', '').lower()
                comps_of_rxn.extend([compartment])
                if compartment not in pc.compartments:
                    pc.compartments.extend([compartment])

            # Multicompartment and transport rxns need to be parsed manually
            # If rxn is multicompartment
            if len(comps_of_rxn) > 1 and (not pparser.metacyc_db.reaction_type(rxn_id) == '|TRANSPORT|'):
                pc.multi_compartment_rxns.extend([rxn_id])
            # Else if reaction is transport
            elif pparser.metacyc_db.reaction_type(rxn_id) == '|TRANSPORT|':
                pc.transport_rxns.extend([rxn_id])
            # Only add non-transport rxns
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
                    rxn.subsystem = rxn.subsystem + ' ' + pthwy_id
                else:
                    # Add reaction
                    rxn_to_add = cobra.Reaction(rxn_id)
                    model.add_reaction(rxn_to_add)

                    # Add Stoichiometry
                    stoichiometry = pc.parse_stoichiometry(model, rxn_id, mets_in_reaction, metabolites, comp)
                    rxn_to_add.add_metabolites(stoichiometry)

                    # Check if metabolites have thermodynamic info
                    mets_w_thermo_info = [met.id not in pc.mets_wo_thermo_info and met.id not in pc.mets_wo_info for met in rxn_to_add.metabolites]

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
                        # And if all metabolites have thermo info
                        if all(mets_w_thermo_info):
                            # set reversibility constrains
                            lower_bound, upper_bound, dG_lb, dG_ub = pc.add_thermo_constraints(rxn_to_add)
                        # if not
                        else:
                            # set reaction as reversible
                            pc.rxns_wo_thermo_info.extend([rxn_to_add.id])
                            dG_lb = []
                            dG_ub = []
                            lower_bound = -100
                            upper_bound = 100
                    # if rxn is not mass balanced
                    else:
                        # parse reaction to one of the categories
                        if all(mets_w_thermo_info):
                            pass
                        else:
                            pc.rxns_wo_thermo_info.extend([rxn_to_add.id])
                        pc.classify_unbalanced_rxn(rxn_to_add)

                        # set rxn as irreversible
                        dG_lb = []
                        dG_ub = []
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
                                             'dG_lb': dG_lb,
                                             'dG_ub': dG_ub}
    else:
        pc.unmatched_pthwys.extend([pthwy_id])

print(str(len(pc.unmatched_pthwys)) + ' pathways do not match')
print(str(len(pc.mets_w_concentration)) + ' metabolites have concentration info')
print(str(len(pc.mets_wo_thermo_info)) + ' metabolites don`t have thermo info')
print(str(len(pc.mets_wo_info)) + ' metabolites have no info')
print(str(len(pc.rxns_wo_thermo_info)) + ' reactions have missing thermo info')
print(str(len(pc.rxns_w_thermo_error)) + ' reactions got an error')
print(str(len(pc.transport_rxns)) + ' reactions are transport reactions')
print(str(len(pc.multi_compartment_rxns)) + ' reactions are multicompartment reactions')
print(str(len(pc.should_be_balanced_rxns)) + ' reactions are wrongly annotated')
print(str(len(pc.one_proton_unbalanced_rxns)) + ' reactions are one H+ unbalanced')
print(str(len(pc.proton_unbalanced_rxns)) + ' reactions are H+ unbalanced')
print(str(len(pc.one_charge_unbalanced_rxns)) + ' reactions are one charge unbalanced')
print(str(len(pc.charge_unbalanced_rxns)) + ' reactions are charge unbalanced')
print(str(len(pc.proton_and_charge_unbalanced_rxns)) + ' reactions are H+ and charge unbalanced')
print(str(len(pc.eqtr_unbalanced_rxns)) + ' reactions are unbalanced')
print(str(len(pc.unbalanced_rxns)) + ' reactions are unbalanced')
print(str(len(pc.corrected_hc_rxns)) + ' reactions had stoichiometry change by nH+')
