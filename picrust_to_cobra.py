# %% codecell
from __future__ import print_function
from corrected_datasets import eqtr_metacyc_map, corrected_stoichiometry, corrected_metabolites
import os
import pandas as pd
import cobra
import pythoncyc as pcyc
import equilibrator_api as eqtr


class idNotFoundError(Exception):
    pass


class dbLinksMissingError(Exception):
    pass


cc = eqtr.ComponentContribution()
cc.p_h = eqtr.Q_(7.4)
cc.p_mg = eqtr.Q_(3.0)
cc.ionic_strength = eqtr.Q_("0.25M")
cc.temperature = eqtr.Q_("298.15K")


metacyc_db = pcyc.select_organism('meta')
all_pathways = metacyc_db.all_pathways()

os.chdir('/home/alexis/UAM/cobra')
eqtr_metacyc_map.set_index('equilibrator', inplace=True)
pathways = pd.read_csv('pathways.csv')
pathways.drop(columns=['Unnamed: 0'], inplace=True)
pathways = pathways.iloc[:, 0].to_list()

conc_of_mets = pd.read_csv('metabolites_w_con.csv')
for i in conc_of_mets.index:
    try:
        conc_of_mets.loc[i, 'InChl'] = conc_of_mets.loc[i, 'InChl'].replace('InChI=', '')
    except AttributeError:
        pass
conc_of_mets.drop_duplicates('InChl', keep='first', inplace=True)
conc_of_mets.dropna(inplace=True)
conc_of_mets.set_index('InChl', inplace=True)

# %% codecell
model = cobra.Model()
mets_wo_info = []
transport_rxns = []
multi_compartment_rxns = []
mets_wo_thermo_info = []
mets_w_concentration = []
rxns_wo_thermo_info = pd.DataFrame()
should_be_balanced_rxns = []
unbalanced_rxns = []
unmatched_pthwys = []
proton_unbalanced_rxns = []


n = 0
step = 0.1
for pthwy in pathways:
    n += 1
    if n / len(pathways) > step:
        step += 0.1
        print(str(n) + ' pathways parsed %3.2f completed' % (n / len(pathways) * 100))

    pthwy_id = '|' + pthwy + '|'
    if pthwy_id in all_pathways:
        for rxn_id in metacyc_db[pthwy_id].reaction_list:
            # List of compartments
            compartments = []
            for comp in metacyc_db.compartments_of_reaction(rxn_id):
                compartments.extend([comp.replace('|', '').replace('CCO-', '').lower()])

            if len(compartments) > 1 and (not metacyc_db.reaction_type(rxn_id) == '|TRANSPORT|'):
                multi_compartment_rxns.extend([rxn_id])
            elif metacyc_db.reaction_type(rxn_id) == '|TRANSPORT|':
                transport_rxns.extend([rxn_id])
            else:
                # Dictionary of reactants and products
                substrates = metacyc_db.reaction_reactants_and_products(rxn_id, pwy=pthwy_id)
                if not substrates:
                    substrates = [metacyc_db[rxn_id].left,
                                  metacyc_db[rxn_id].right]

                mets_in_reaction = dict(zip(['reactants', 'products'],
                                            substrates))

                # List of metabolites
                metabolites = mets_in_reaction['reactants'].copy()
                metabolites.extend(mets_in_reaction['products'])
                comp = '_' + compartments[0][0:2]

                for met in metabolites:
                    met_id = met + comp
                    if met_id in model.metabolites:
                        pass
                    else:
                        thermo_info = False
                        inchi_flag = False
                        db_flag = False
                        if metacyc_db[met].inchi_key:
                            inchi_flag = True
                            eqtr_met = cc.search_compound_by_inchi_key(metacyc_db[met].inchi_key.replace('InChIKey=', ''))
                            if eqtr_met:
                                eqtr_met = cc.get_compound_by_internal_id(eqtr_met[0].id)
                                thermo_info = True
                            else:
                                try:
                                    if metacyc_db[met].dblinks:
                                        db_flag = True
                                    else:
                                        raise dbLinksMissingError
                                except dbLinksMissingError:
                                    db_flag = False
                                    thermo_info = False

                        else:
                            inchi_flag = False
                            try:
                                if metacyc_db[met].dblinks:
                                    db_flag = True
                                else:
                                    raise dbLinksMissingError
                            except dbLinksMissingError:
                                db_flag = False
                                thermo_info = False

                        if db_flag:
                            filter_db = [db in metacyc_db[met].dblinks.keys() for db in eqtr_metacyc_map['metacyc']]
                            dblinks = {}
                            for db in eqtr_metacyc_map.loc[filter_db, 'metacyc'].index:
                                dblinks[db] = metacyc_db[met].dblinks[eqtr_metacyc_map.loc[db, 'metacyc']][0]

                            for db in dblinks.keys():
                                cid = db + ':' + dblinks[db]
                                eqtr_met = cc.get_compound(cid)
                                if eqtr_met:
                                    thermo_info = True
                                    break
                                else:
                                    thermo_info = False
                        elif not db_flag:
                            dblinks = {}

                        if thermo_info:
                            if metacyc_db[met].chemical_formula:
                                elemnts = metacyc_db[met].chemical_formula.copy()
                                formula = ''
                                for key in elemnts.keys():
                                    formula += key.replace('|', '') + str(elemnts[key][0])
                            else:
                                formula = eqtr_met.formula
                            charge = eqtr_met.net_charge
                            name = eqtr_met.get_common_name()
                            smiles = eqtr_met.smiles
                            dblinks = {}
                            database_old = []
                            for link in eqtr_met.identifiers:
                                database_new = link.registry.namespace
                                if database_new == database_old:
                                    pass
                                else:
                                    dblinks[database_new] = link.accession
                                database_old = database_new
                        elif not thermo_info:
                            if metacyc_db[met].chemical_formula:
                                elemnts = metacyc_db[met].chemical_formula.copy()
                                charge = []
                                name = metacyc_db[met].common_name
                                smiles = metacyc_db[met].smiles
                                formula = ''
                                for key in elemnts.keys():
                                    formula += key.replace('|', '') + str(elemnts[key][0])
                                mets_wo_thermo_info.extend([met_id])
                            else:
                                charge = []
                                name = metacyc_db[met].common_name
                                smiles = metacyc_db[met].smiles
                                formula = ''
                                mets_wo_info.extend([met_id])

                        if inchi_flag:
                            inchi = metacyc_db[met].inchi.replace('InChI=', '')
                            inchi_key = metacyc_db[met].inchi_key.replace('InChIKey=', '')
                        else:
                            inchi = ''
                            inchi_key = ''

                        if inchi in conc_of_mets.index:
                            mets_w_concentration.extend([met_id])
                            lb = conc_of_mets.loc[inchi, 'lb']
                            ub = conc_of_mets.loc[inchi, 'ub']
                            mean = conc_of_mets.loc[inchi, 'lb']
                        else:
                            lb = 1e-5
                            ub = 0.01
                            mean = 0

                        if not formula:
                            if smiles:
                                formula = str('C' + str(smiles.count('C')) +
                                              'H' + str(smiles.count('H')) +
                                              'O' + str(smiles.count('O')) +
                                              'N' + str(smiles.count('N')) +
                                              'S' + str(smiles.count('S'))
                                              )
                        if met in corrected_metabolites.keys():
                            formula = corrected_metabolites[met]['formula']

                        structure_atoms = metacyc_db[met].structure_atoms
                        structure_bonds = metacyc_db[met].structure_bonds
                        met_to_add = cobra.Metabolite(met_id,
                                                      formula=formula,
                                                      compartment=compartments[0],
                                                      name=name)
                        met_to_add.annotation = {'dblinks': dblinks,
                                                 'InChI': inchi,
                                                 'InChIKey': inchi_key,
                                                 'smiles': smiles,
                                                 'structure_atoms': structure_atoms,
                                                 'structure_bonds': structure_bonds,
                                                 'lb': lb,
                                                 'ub': ub,
                                                 'mean': mean}

                        model.add_metabolites(met_to_add)

                metabolites = [met + comp for met in metabolites]
                if rxn_id in model.reactions:
                    rxn = model.reactions.get_by_id(rxn_id)
                    rxn.subsystem = rxn.subsystem + ' ' + pthwy_id
                else:
                    # Add reaction
                    rxn_to_add = cobra.Reaction(rxn_id)
                    model.add_reaction(rxn_to_add)

                    # Add Stoichiometry
                    if rxn_id in corrected_stoichiometry.keys():
                        stoichiometry = corrected_stoichiometry[rxn_id]['stoichiometry']
                    else:
                        coefficients = []
                        for susbstrate in mets_in_reaction['reactants']:
                            coefficients.extend([-1])
                        for product in mets_in_reaction['products']:
                            coefficients.extend([1])
                        stoichiometry = dict(zip(metabolites,
                                                 coefficients))
                    rxn_to_add.add_metabolites(stoichiometry)

                    # Check if metabolites have thermodynamic info
                    mets_w_thermo_info = [met.id not in mets_wo_thermo_info and met.id not in mets_wo_info for met in rxn_to_add.metabolites]

                    # Contra-intuitive check-mass_balance returns false when reaction is mass_balanced
                    mass_balanced = rxn_to_add.check_mass_balance()
                    if not mass_balanced:
                        # if all metabolites have thermo info
                        if all(mets_w_thermo_info):
                            met_conc_lb = []
                            met_conc_ub = []

                            reactants_string = []
                            for reactant in rxn_to_add.reactants:
                                coeff = str(-stoichiometry[reactant.id])
                                try:
                                    db = list(reactant.annotation['dblinks'].keys())[0]
                                    db_id = reactant.annotation['dblinks'][db]
                                    db_link = db + ':' + db_id
                                    compound = cc.get_compound(db_link)
                                    if not compound:
                                        raise idNotFoundError
                                except idNotFoundError:
                                    found = []
                                    print(db_link + ' not found in compund cache, searching by next id')
                                    for db in reactant.annotation['dblinks'].keys():
                                        db_id = reactant.annotation['dblinks'][db]
                                        db_link = db + ':' + db_id
                                        compound = cc.get_compound(db_link)
                                        if compound:
                                            print('compund found')
                                            found = True
                                            break
                                        else:
                                            found = False
                                    if not found:
                                        print(db_link + ' not found by any id')

                                reactants_string.extend([coeff + ' ' + db_link])
                                met_conc_lb.extend([(db_link, reactant.annotation['lb'])])
                                met_conc_ub.extend([(db_link, reactant.annotation['ub'])])

                            products_string = []
                            for product in rxn_to_add.products:
                                try:
                                    db = list(product.annotation['dblinks'].keys())[0]
                                    db_id = product.annotation['dblinks'][db]
                                    db_link = db + ':' + db_id
                                    compound = cc.get_compound(db_link)
                                    if not compound:
                                        raise idNotFoundError
                                except idNotFoundError:
                                    found = []
                                    print(db_link + ' not found in compund cache, searching by next id')
                                    for db in product.annotation['dblinks'].keys():
                                        db_id = product.annotation['dblinks'][db]
                                        db_link = db + ':' + db_id
                                        compound = cc.get_compound(db_link)
                                        if compound:
                                            print('compund found')
                                            found = True
                                            break
                                        else:
                                            found = False
                                    if not found:
                                        print(db_link + ' not found by any id')

                                products_string.extend([coeff + ' ' + db_link])
                                met_conc_lb.extend([(db_link, product.annotation['ub'])])
                                met_conc_ub.extend([(db_link, product.annotation['lb'])])

                            # Calculate dG_low and dG_high
                            reaction_string = ' + '.join(reactants_string) + \
                                ' -> ' + ' + '.join(products_string)
                            eqtr_rxn = cc.parse_reaction_formula(reaction_string)

                            for cid, conc in met_conc_lb:
                                try:
                                    compound = cc.get_compound(cid)
                                    eqtr_rxn.set_phase(compound, 'aqueous')
                                except:
                                    met_conc_lb.remove((cid, conc))

                            for cid, conc in met_conc_ub:
                                try:
                                    compound = cc.get_compound(cid)
                                    eqtr_rxn.set_phase(compound, 'aqueous')
                                except:
                                    met_conc_ub.remove((cid, conc))

                            for cid, conc in met_conc_lb:
                                compound = cc.get_compound(cid)
                                if compound.formula == 'H':
                                    pass
                                else:
                                    abundance = eqtr.Q_(conc, "M")
                                    eqtr_rxn.set_abundance(compound, abundance)

                            try:
                                dG_lb = cc.dg_prime(eqtr_rxn)
                                dG_flag = True
                            except Exception as e:
                                print(e)
                                woti = pd.Series([met.id for met in rxn_to_add.metabolites])
                                rxns_wo_thermo_info.loc[rxn_id, 'metabolites'] = ' '.join(
                                    woti.loc[[not flag for flag in mets_w_thermo_info]])
                                dG_lb = []
                                dG_flag = False

                            for cid, conc in met_conc_ub:
                                compound = cc.get_compound(cid)
                                if compound.formula == 'H':
                                    pass
                                else:
                                    abundance = eqtr.Q_(conc, "M")
                                    eqtr_rxn.set_abundance(compound, abundance)
                            try:
                                dG_ub = cc.dg_prime(eqtr_rxn)
                            except:
                                pass

                            # Add reversibility constrains
                            if dG_flag:
                                if (dG_lb < 0) and (dG_lb < 0):
                                    lower_bound = 0
                                    upper_bound = 100
                                elif (dG_lb > 0) and (dG_lb > 0):
                                    lower_bound = -100
                                    upper_bound = 0
                                else:
                                    lower_bound = -100
                                    upper_bound = 100
                            else:
                                lower_bound = -100
                                upper_bound = 100
                        else:
                            woti = pd.Series([met.id for met in rxn_to_add.metabolites])
                            rxns_wo_thermo_info.loc[rxn_id, 'metabolites'] = ' '.join(
                                woti.loc[[not flag for flag in mets_w_thermo_info]])
                            dG_lb = []
                            dG_ub = []
                            lower_bound = -100
                            upper_bound = 100

                    else:
                        if all(mets_w_thermo_info):
                            pass
                        else:
                            woti = pd.Series([met.id for met in rxn_to_add.metabolites])
                            rxns_wo_thermo_info.loc[rxn_id, 'metabolites'] = ' '.join(
                                woti.loc[[not flag for flag in mets_w_thermo_info]])
                        unbalanced_elements = list(mass_balanced.keys())
                        if len(unbalanced_elements) == 1 and unbalanced_elements[0] == 'H':
                            proton_unbalanced_rxns.extend([rxn_id])
                        elif metacyc_db[rxn_id].reaction_balance_status == '|BALANCED|':
                            should_be_balanced_rxns.extend([rxn_id])
                        else:
                            unbalanced_rxns.extend([rxn_id])

                        dG_lb = []
                        dG_ub = []
                        lower_bound = -100
                        upper_bound = 100

                    rxn_to_add.lower_bound = lower_bound
                    rxn_to_add.upper_bound = upper_bound
                    genes_of_reaction = [s.replace('|', '')
                                         for s in metacyc_db.genes_of_reaction(rxn_id)]
                    if genes_of_reaction:
                        rxn_to_add.gene_reaction_rule = '( ' + ' or '.join(genes_of_reaction) + ' )'

                    rxn_to_add.subsystem = pthwy_id
                    atom_mappings = metacyc_db[rxn_id].atom_mappings
                    ec_number = metacyc_db[rxn_id].ec_number
                    dblinks = metacyc_db[rxn_id].dblinks
                    rxn_to_add.annotation = {'atom_mappings': atom_mappings,
                                             'ec_number': ec_number,
                                             'dblinks': dblinks,
                                             'dG_lb': dG_lb,
                                             'dG_ub': dG_ub}
    else:
        unmatched_pthwys.extend([pthwy_id])

print(str(len(unmatched_pthwys)) + ' pathways do not match')
print(str(len(mets_w_concentration)) + ' metabolites have concentration info')
print(str(len(mets_wo_thermo_info)) + ' metabolites don`t have thermo info')
print(str(len(mets_wo_info)) + ' metabolites have no info')
print(str(len(rxns_wo_thermo_info)) + ' reactions have missing thermo info')
print(str(len(transport_rxns)) + ' reactions are transport reactions')
print(str(len(multi_compartment_rxns)) + ' reactions are multicompartment reactions')
print(str(len(should_be_balanced_rxns)) + ' reactions are wrongly annotated')
print(str(len(proton_unbalanced_rxns)) + ' reactions are H+ unbalanced')
