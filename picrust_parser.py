from __future__ import print_function

import cobra
from cobra.flux_analysis.helpers import normalize_cutoff
from cobra.flux_analysis.fastcc import _find_sparse_mode, _flip_coefficients
import os
import pandas as pd

import pythoncyc as pcyc
import equilibrator_api as eqtr
os.chdir('/home/alexis/UAM/cobra')
from data_files.corrected_datasets import eqtr_metacyc_map, corrected_stoichiometry, corrected_metabolites, added_pthwys
eqtr_metacyc_map.set_index('equilibrator', inplace=True)

# Load metacyc
metacyc_db = pcyc.select_organism('meta')
all_pathways = metacyc_db.all_pathways()

# Prepare componenet contribution
cc = eqtr.ComponentContribution()
cc.p_h = eqtr.Q_(7.5)
cc.p_mg = eqtr.Q_(3.0)
cc.ionic_strength = eqtr.Q_("0.25M")
cc.temperature = eqtr.Q_("298.15K")

# Load metabolites with concentration
conc_of_mets = pd.read_csv('data_files/metabolites_w_con.csv')
for i in conc_of_mets.index:
    try:
        conc_of_mets.loc[i, 'InChl'] = conc_of_mets.loc[i, 'InChl'].replace('InChI=', '')
    except AttributeError:
        pass
conc_of_mets.drop_duplicates('InChl', keep='first', inplace=True)
conc_of_mets.dropna(inplace=True)
conc_of_mets.set_index('InChl', inplace=True)


class idNotFoundError(Exception):
    pass


class dbLinksMissingError(Exception):
    pass


class p_model():

    def __init__(self, path_to_pthwys='data_files/pathways.csv'):
        self.transport_rxns = []
        self.multi_compartment_rxns = []
        self.mets_w_concentration = []
        self.rxns_w_thermo_error = []
        self.should_be_balanced_rxns = []
        self.unbalanced_rxns = []
        self.unmatched_pthwys = []
        self.one_proton_unbalanced_rxns = []
        self.proton_unbalanced_rxns = []
        self.one_charge_unbalanced_rxns = []
        self.charge_unbalanced_rxns = []
        self.proton_and_charge_unbalanced_rxns = []
        self.corrected_hc_rxns = []
        self.compartments = []
        self.eqtr_unbalanced_rxns = []
        self.pathways = pd.read_csv(path_to_pthwys)
        self.pathways.drop(columns=['Unnamed: 0'], inplace=True)
        self.pathways.drop_duplicates(inplace = True)
        self.pathways = self.pathways.iloc[:, 0].to_list()
        self.pathways.extend(added_pthwys)

    def replace_redox_pairs(self, metabolites):
        c_redox_pairs = [met == '|NADH-P-OR-NOP|' or met == '|NADH-P-OR-NOP|' or met == '|Reduced-2Fe-2S-Ferredoxins|' or met ==
                         '|Oxidized-2Fe-2S-Ferredoxins|' or met == '|METHYLENE-THF-GLU-N|' or met == '|5-METHYL-THF-GLU-N|' or met ==
                         '|5-10-METHENYL-THF-GLU-N|' or met == '|THF-GLU-N|' or met ==
                         '|FORMYL-THF-GLU-N|' or met == '|N5-Formyl-THF-Glu-N|' or met ==
                         '|Reduced-Flavoproteins|' or met == '|Oxidized-Flavoproteins|' or met ==
                         '|DIHYDROFOLATE-GLU-N|_cy' for met in metabolites]
        if any(c_redox_pairs):
            new_metabolites = []
            for met in metabolites:
                if met == '|NADH-P-OR-NOP|':
                    new_metabolites.extend(['|NADH|'])
                elif met == '|NAD-P-OR-NOP|':
                    new_metabolites.extend(['|NAD|'])
                elif met == '|Reduced-2Fe-2S-Ferredoxins|':
                    new_metabolites.extend(['|Reduced-ferredoxins|'])
                elif met == '|Oxidized-2Fe-2S-Ferredoxins|':
                    new_metabolites.extend(['|Oxidized-ferredoxins|'])
                elif met == '|METHYLENE-THF-GLU-N|':
                    new_metabolites.extend(['|METHYLENE-THF|'])
                elif met == '|5-METHYL-THF-GLU-N|':
                    new_metabolites.extend(['|5-METHYL-THF|'])
                elif met == '|5-10-METHENYL-THF-GLU-N|':
                    new_metabolites.extend(['|5-10-METHENYL-THF|'])
                elif met == '|THF-GLU-N|':
                    new_metabolites.extend(['|THF|'])
                elif met == '|FORMYL-THF-GLU-N|':
                    new_metabolites.extend(['|10-FORMYL-THF|'])
                elif met == '|N5-Formyl-THF-Glu-N|':
                    new_metabolites.extend(['|5-FORMYL-THF|'])
                elif met == '|Reduced-Flavoproteins|':
                    new_metabolites.extend(['|FADH2|'])
                elif met == '|Oxidized-Flavoproteins|':
                    new_metabolites.extend(['|FAD|'])
                elif met == '|DIHYDROFOLATE-GLU-N|':
                    new_metabolites.extend(['|DIHYDROFOLATE|'])
                else:
                    new_metabolites.extend([met])
        else:
            new_metabolites = metabolites

        return new_metabolites

    def parse_metabolite(self, met, comp):
        met_id = met + comp
        thermo_info = False
        inchi_flag = False
        db_flag = False

        # Check if metabolite has thermodynamic information
        if metacyc_db[met].inchi_key:
            # First check by inchi key
            inchi_flag = True
            eqtr_met = cc.search_compound_by_inchi_key(metacyc_db[met].inchi_key.replace('InChIKey=', ''))
            if eqtr_met:
                eqtr_met = cc.get_compound_by_internal_id(eqtr_met[0].id)
                thermo_info = True
            # if it has inchi but compound can't be found, check if it has database links
            else:
                try:
                    if metacyc_db[met].dblinks:
                        db_flag = True
                    else:
                        raise dbLinksMissingError
                # If it has no database links, raise error
                except dbLinksMissingError:
                    db_flag = False
                    thermo_info = False
        # if met does not have inchi_key
        else:
            inchi_flag = False
            # look for databe links
            try:
                if metacyc_db[met].dblinks:
                    db_flag = True
                else:
                    raise dbLinksMissingError
            except dbLinksMissingError:
                db_flag = False
                thermo_info = False

        # if met does not have, or can't be found by inhci_key
        if db_flag:
            # parse db_links in metacyc to equilibrator compatible string
            filter_db = [db in metacyc_db[met].dblinks.keys() for db in eqtr_metacyc_map['metacyc']]
            dblinks = {}
            for db in eqtr_metacyc_map.loc[filter_db, 'metacyc'].index:
                dblinks[db] = metacyc_db[met].dblinks[eqtr_metacyc_map.loc[db, 'metacyc']][0]
            # Search for thermo_info with all db_links
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

        # if met has thermo_info
        if thermo_info:
            # Parse formula from metacyc
            if metacyc_db[met].chemical_formula:
                elemnts = metacyc_db[met].chemical_formula.copy()
                formula = ''
                for key in elemnts.keys():
                    formula += key.replace('|', '') + str(elemnts[key][0])
            # Or equilibrator
            else:
                formula = eqtr_met.formula
            # Parse rest of info from equilibrator
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
        # If met is not in equilibrator
        elif not thermo_info:
            # if formula is in metacyc
            if metacyc_db[met].chemical_formula:
                # Parse all info from metacyc
                elemnts = metacyc_db[met].chemical_formula.copy()
                charge = None
                name = metacyc_db[met].common_name
                smiles = metacyc_db[met].smiles
                formula = ''
                for key in elemnts.keys():
                    formula += key.replace('|', '') + str(elemnts[key][0])
            # if no formula is found
            else:
                # Then metabolite has no info
                charge = None
                name = metacyc_db[met].common_name
                smiles = metacyc_db[met].smiles
                formula = None

        # Parse inchi string and inchi key
        if inchi_flag:
            inchi = metacyc_db[met].inchi.replace('InChI=', '')
            inchi_key = metacyc_db[met].inchi_key.replace('InChIKey=', '')
        else:
            inchi = None
            inchi_key = None

        # look if met has concentration info by inchi string
        if inchi in conc_of_mets.index:
            self.mets_w_concentration.extend([met_id])
            lb = conc_of_mets.loc[inchi, 'lb']
            ub = conc_of_mets.loc[inchi, 'ub']
            mean = conc_of_mets.loc[inchi, 'lb']
        # if it does not add standard concenctration
        else:
            lb = 3e-6
            ub = 0.003
            mean = 0

        # parse atomic info
        structure_atoms = metacyc_db[met].structure_atoms
        structure_bonds = metacyc_db[met].structure_bonds
        # if up to this point met has no formula
        if not formula:
            # parse formula from atomic info
            if structure_atoms:
                formula = str('C' + str(str(structure_atoms).count('C')) +
                              'H' + str(str(structure_atoms).count('H')) +
                              'O' + str(str(structure_atoms).count('O')) +
                              'N' + str(str(structure_atoms).count('N')) +
                              'S' + str(str(structure_atoms).count('S')) +
                              'P' + str(str(structure_atoms).count('P'))
                              )
            # Or smiles
            elif smiles:
                formula = str('C' + str(smiles.count('C')) +
                              'H' + str(smiles.count('H')) +
                              'O' + str(smiles.count('O')) +
                              'N' + str(smiles.count('N')) +
                              'S' + str(smiles.count('S')) +
                              'P' + str(smiles.count('P'))
                              )
        # if metabolite is in the list of mets with corrected info
        if met in corrected_metabolites.keys():
            # set correct charge and formula
            met_to_correct = corrected_metabolites[met]
            if 'formula' in met_to_correct.keys():
                formula = met_to_correct['formula']
            if 'charge' in met_to_correct.keys():
                charge = met_to_correct['charge']

        metabolite_to_add = cobra.Metabolite(met_id,
                                             formula=formula,
                                             charge=charge,
                                             name=name
                                             )
        metabolite_to_add.annotation = {'dblinks': dblinks,
                                        'InChI': inchi,
                                        'InChIKey': inchi_key,
                                        'smiles': smiles,
                                        'structure_atoms': structure_atoms,
                                        'structure_bonds': structure_bonds,
                                        'lb': lb,
                                        'ub': ub,
                                        'mean': mean
                                        }
        return metabolite_to_add

    def parse_stoichiometry(self, model, rxn_id, mets_in_reaction, metabolites):
        # if rxn in corrected reactions
        if rxn_id in corrected_stoichiometry.keys():
            # Set stoichiometry as indicated
            stoichiometry = corrected_stoichiometry[rxn_id]['stoichiometry']
            # check if mets in corrected have been added to the model
            for met_id in stoichiometry.keys():
                if met_id not in model.metabolites:
                    comp = met_id[-3:]
                    if comp == '_cy':
                        compartment = 'cytosol'
                    elif comp == '_pe':
                        compartment = 'periplasm'
                    elif comp == '_ex':
                        compartment = 'extracellular'
                    elif comp == 'it':
                        compartment = 'intramembrane'
                    met_to_add = self.parse_metabolite(met_id[:-3],
                                                       comp)
                    met_to_add.compartment = compartment
                    model.add_metabolites(met_to_add)
        # if rxn is not in corrected dataset
        else:
            # Set all coefficients to 1
            coefficients = []
            for susbstrate in mets_in_reaction['reactants']:
                coefficients.extend([-1])
            for product in mets_in_reaction['products']:
                coefficients.extend([1])
            stoichiometry = dict(zip(metabolites,
                                     coefficients)
                                 )
        return stoichiometry

    def dblink_to_equilibrator_id(self, metabolite):
        try:
            databases = metabolite.annotation['dblinks']
            if databases:
                for db in databases.keys():
                    # get id
                    db_id = metabolite.annotation['dblinks'][db]
                    # write database:met_id
                    db_link = db + ':' + db_id
                    # Look for metabolite in equilibrator cache
                    compound = cc.get_compound(db_link)
                    if compound:
                        found = True
                        break
                    else:
                        # print(db_link + ' not found in compound cache, searching by next id')
                        found = False
            else:
                found = False
                db_link = ''
            if not found:
                raise idNotFoundError
        except idNotFoundError:
            if db_link:
                pass
                # print(db_link + ' not found by any id')
            else:
                pass
                # print(metabolite.id + ' has no database links')

        return db_link

    def cobra_to_eqtr_rxn(self, rxn_to_add):
        met_conc_lb = []
        met_conc_ub = []

        reactants_string = []
        for reactant in rxn_to_add.reactants:
            coeff = str(-rxn_to_add.get_coefficient(reactant))
            db_link = self.dblink_to_equilibrator_id(reactant)
            reactants_string.extend([coeff + ' ' + db_link])
            met_conc_lb.extend([(db_link, reactant.annotation['lb'])])
            met_conc_ub.extend([(db_link, reactant.annotation['ub'])])

        products_string = []
        for product in rxn_to_add.products:
            coeff = str(rxn_to_add.get_coefficient(product))
            db_link = self.dblink_to_equilibrator_id(product)
            products_string.extend([coeff + ' ' + db_link])
            # driving force of reaction is lowest whe products are high
            met_conc_lb.extend([(db_link, product.annotation['ub'])])
            # and highest when products are low
            met_conc_ub.extend([(db_link, product.annotation['lb'])])

        # write reaction string
        reaction_string = ' + '.join(reactants_string) + \
            ' -> ' + ' + '.join(products_string)

        return met_conc_lb, met_conc_ub, reaction_string

    def calculate_dg_bound(self, rxn_to_add, eqtr_rxn, concentration_bounds):
        # set concentration bounds
        for c_id, cc_bound in concentration_bounds:
            compound = cc.get_compound(c_id)
            try:
                # if compound in aqeuos phase
                eqtr_rxn.set_phase(compound, 'aqueous')
                if compound.formula == 'H':
                    pass
                else:
                    # set concentration
                    abundance = eqtr.Q_(cc_bound, "M")
                    eqtr_rxn.set_abundance(compound, abundance)
            except:
                # else, do nothing
                pass

        # Double check mass balance
        try:
            balanced = eqtr_rxn.is_balanced(raise_exception = True)
        except:
            pass
        else:
            if not balanced:
                if rxn_to_add.id not in self.eqtr_unbalanced_rxns:
                    self.eqtr_unbalanced_rxns.extend([rxn_to_add.id])

        # Calculate dg bound
        try:
            dG_bound = cc.dg_prime(eqtr_rxn)
            rev_index = cc.ln_reversibility_index(eqtr_rxn)
        # if error
        except Exception:
            # print(e)
            if rxn_to_add.id not in self.rxns_w_thermo_error:
                self.rxns_w_thermo_error.extend([rxn_to_add.id])
            dG_bound = []
            rev_index = []

        return dG_bound, rev_index

    def add_thermo_constraints(self, rxn_to_add):
        met_conc_lb, met_conc_ub, reaction_string = self.cobra_to_eqtr_rxn(rxn_to_add)
        try:
            eqtr_rxn = cc.parse_reaction_formula(reaction_string)
            dG_lb, ri_lb = self.calculate_dg_bound(rxn_to_add, eqtr_rxn, met_conc_lb)
            dG_ub, ri_ub = self.calculate_dg_bound(rxn_to_add, eqtr_rxn, met_conc_ub)
        except:
            if rxn_to_add.id not in self.rxns_w_thermo_error:
                self.rxns_w_thermo_error.extend([rxn_to_add.id])
            dG_lb = []
            dG_ub = []
            ri_lb = []
            ri_ub = []

        # Add reversibility constrains
        if all([dG_lb, dG_ub]):
            if (dG_lb.value - dG_lb.error) < 0 and (dG_ub.value + dG_ub.error) < 0:
                lower_bound = 0
                upper_bound = 100
            elif (dG_lb.value - dG_lb.error) > 0 and (dG_ub.value + dG_ub.error) > 0:
                lower_bound = -100
                upper_bound = 0
            else:
                lower_bound = -100
                upper_bound = 100
        else:
            lower_bound = -100
            upper_bound = 100

        return lower_bound, upper_bound, dG_lb, dG_ub, ri_lb, ri_ub

    def classify_unbalanced_rxn(self, rxn_to_add):
        mass_balanced = rxn_to_add.check_mass_balance()
        unbalanced_elements = list(mass_balanced.keys())

        if len(unbalanced_elements) == 1 and unbalanced_elements[0] == 'H':
            if abs(mass_balanced['H']) == 1:
                self.one_proton_unbalanced_rxns.extend([rxn_to_add.id])
            else:
                self.proton_unbalanced_rxns.extend([rxn_to_add.id])
        elif len(unbalanced_elements) == 1 and unbalanced_elements[0] == 'charge':
            if abs(mass_balanced['charge']) == 1:
                self.one_charge_unbalanced_rxns.extend([rxn_to_add.id])
            else:
                self.charge_unbalanced_rxns.extend([rxn_to_add.id])
        elif len(unbalanced_elements) == 2 and ('charge' in unbalanced_elements and 'H' in unbalanced_elements):
            self.proton_and_charge_unbalanced_rxns.extend([rxn_to_add.id])
        elif metacyc_db[rxn_to_add.id].reaction_balance_status == '|BALANCED|':
            self.should_be_balanced_rxns.extend([rxn_to_add.id])
        else:
            self.unbalanced_rxns.extend([rxn_to_add.id])

    def add_rxn_from_stoichiometry(self, model, pthwy_id, rxn_id, stoichimetry,
                                   reversible = True, direction = None, get_thermo_constraints = False):
        rxn_to_add = cobra.Reaction(rxn_id)
        model.add_reaction(rxn_to_add)

        for met_id in stoichimetry.keys():
            if met_id not in model.metabolites:
                comp = met_id[-3:]
                if comp == '_cy':
                    compartment = 'cytosol'
                elif comp == '_pe':
                    compartment = 'periplasm'
                elif comp == '_ex':
                    compartment = 'extracellular'
                elif comp == '_it':
                    compartment = 'intramembrane'
                else:
                    print(met_id + ' in ' + rxn_id + ' has no compartment specified')

                met_to_add = self.parse_metabolite(met_id[:-3], comp)
                met_to_add.compartment = compartment
                model.add_metabolites(met_to_add)
        rxn_to_add.add_metabolites(stoichimetry)

        if get_thermo_constraints:
            lower_bound, upper_bound, dG_lb, dG_ub, ri_lb, ri_ub = self.add_thermo_constraints(rxn_to_add)
        else:
            dG_lb = []
            dG_ub = []
            ri_lb = []
            ri_ub = []
            if reversible:
                lower_bound = -100
                upper_bound = 100
            elif not reversible:
                if direction == 'L2R':
                    lower_bound = 0
                    upper_bound = 100
                elif direction == 'R2L':
                    lower_bound = -100
                    upper_bound = 0
                else:
                    print('Reaction direction should de L2R or R2L')

        rxn_to_add.lower_bound = lower_bound
        rxn_to_add.upper_bound = upper_bound
        try:
            genes_of_reaction = [gene.replace('|', '')
                                 for gene in metacyc_db.genes_of_reaction(rxn_id[rxn_id.index('|'):])]
        except:
            genes_of_reaction = []

        if genes_of_reaction:
            rxn_to_add.gene_reaction_rule = '( ' + ' or '.join(genes_of_reaction) + ' )'

        rxn_to_add.subsystem = pthwy_id
        try:
            atom_mappings = metacyc_db[rxn_id[rxn_id.index('|'):]].atom_mappings
            ec_number = metacyc_db[rxn_id[rxn_id.index('|'):]].ec_number
            dblinks = metacyc_db[rxn_id[rxn_id.index('|'):]].dblinks
            rxn_to_add.annotation = {'atom_mappings': atom_mappings,
                                     'ec_number': ec_number,
                                     'dblinks': dblinks,
                                     'dG': [dG_lb, dG_ub],
                                     'rIndex': [ri_lb, ri_ub]}
        except:
            pass
        try:
            self.classify_unbalanced_rxn(rxn_to_add)
        except:
            pass

    def search_biomass_components(self, model, path_to_biomass = 'data_files/biomass.csv'):
        self.biomass_in_model = pd.DataFrame(
            columns = ['class', 'name', 'mmol/g', 'met_id']
            )
        biomas_db = pd.read_csv(path_to_biomass)
        biomas_db
        biomas_db.set_index(
            'metacyc_identifier', inplace = True
            )

        for met in biomas_db.index:
            for metabolite in model.metabolites:
                if '|' + met + '|' in metabolite.id:
                    self.biomass_in_model = self.biomass_in_model.append(
                        pd.DataFrame(
                            [[biomas_db.loc[met, 'class'],
                              biomas_db.loc[met, 'name'],
                              biomas_db.loc[met, 'mmol/g'],
                              metabolite.id]],
                            columns = ['class', 'name', 'mmol/g', 'met_id']
                            )
                        )
        self.biomass_in_model.reset_index(
            inplace = True
            )
        self.biomass_in_model.drop(
            columns=['index'], inplace=True
            )
        dup_filter = self.biomass_in_model['name'].duplicated()
        self.biomass_in_model.drop(
            self.biomass_in_model[dup_filter].index, inplace = True
            )

        return self.biomass_in_model

    def add_biomass_rxn(self, model, biomass_in_model, ATP_growth = 0.0063):
        BM_weight = 0
        # For Every Macromolecule
        for tp in biomass_in_model.loc[:, 'class'].unique():
            # Create Metabolite for Macromolecule (AAs, FAMES, etc.)
            met_to_add = cobra.Metabolite(
                tp, name = tp
            )
            model.add_metabolites(met_to_add)

            # Filter metabolites belonging to Macromolecule
            class_filt = biomass_in_model.loc[:, 'class'] == tp
            components = biomass_in_model.loc[class_filt, 'met_id'].to_list()
            coeffs = (-biomass_in_model.loc[class_filt, 'mmol/g'] / 1000).to_list()
            components.extend([tp])
            coeffs.extend([1.0])
            # Create Dictionary with Stoichiometry
            stoichiometry = dict(zip(components, coeffs))

            # Create reaction
            rxn_to_add = cobra.Reaction(
                'BM_' + tp
                )
            model.add_reaction(
                rxn_to_add
            )
            # Add Stoichiometry
            rxn_to_add.add_metabolites(
                stoichiometry
            )
            rxn_to_add.lower_bound = 0
            rxn_to_add.upper_bound = 1000

            # Calculate Macromolecule Formula from elemental balance
            element_balance = rxn_to_add.check_mass_balance()
            formula = ''
            for e in element_balance:
                if not e == 'charge':
                    if not e == 'F':
                        if not e == 'E':
                            formula += str(e) + str(format(-element_balance[e], '.30f'))
                elif e == 'charge':
                    met_to_add.charge = element_balance[e]
            met_to_add.formula = formula

            # Add macromolecule weight to total biomass weight
            BM_weight += met_to_add.formula_weight

        # Create biomass metabolite
        met_to_add = cobra.Metabolite(
            'biomass', name = "biomass"
            )
        model.add_metabolites(
            met_to_add
            )

        # Add stoichimetry for every macromolecule normalized to BM_weight so that
        # biomass molecular weight = 1g/mol
        biomas_coeff = {
            com: -1 / BM_weight for com in biomass_in_model.loc[:, 'class'].unique()
            }
        biomas_coeff['biomass'] = 1
        # Create Biomass reaction
        rxn_to_add = cobra.Reaction(
            'BIOMASS'
            )
        model.add_reaction(
            rxn_to_add
            )
        # Add Stoichiometry
        rxn_to_add.add_metabolites(
            biomas_coeff
            )
        rxn_to_add.upper_bound = 1000

        # Calculate Biomass formula from element balance
        element_balance = rxn_to_add.check_mass_balance()
        formula = ''
        for e in element_balance:
            if not e == 'charge':
                formula += str(e) + str(format(-element_balance[e], '.30f'))
        met_to_add.formula = formula
        # Add ATP requirement
        rxn_to_add.add_metabolites(
            {'|ATP|_cy': -ATP_growth, '|WATER|_cy': -ATP_growth, '|ADP|_cy': ATP_growth, '|Pi|_cy': ATP_growth}
            )
        # Add NAD(P)H and PROTON requirements
        element_balance = rxn_to_add.check_mass_balance()
        rxn_to_add.add_metabolites(
            {'|NADPH|_cy': -0.001, '|PROTON|_cy': -element_balance['charge'], '|NADP|_cy': 0.001, '|NADH|_cy': -0.002, '|NAD|_cy': 0.002}
            )

        return model

    def run_blocked_demand(self, model, demands_to_test):
        for i in tqdm(range(int(len(demands_to_test)))):
            comp = demands_to_test[i]
            with model as model:
                comp_to_test = model.add_boundary(model.metabolites.get_by_id(comp), type = 'demand')
                model.objective = comp_to_test
                self.solution = model.optimize()
                if self.solution.objective_value < model.tolerance:
                    blocked = cobra.flux_analysis.find_blocked_reactions(model)
                    self.blocked_demands.append([comp, self.solution.objective_value, self.solution.fluxes['|RXN-12219|']])
                    for rxn in blocked:
                        if rxn not in self.blocked_rxns:
                            self.blocked_rxns.extend([rxn])
                else:
                    self.available_demands.append([comp, self.solution.objective_value, self.solution.fluxes['|RXN-12219|']])

    def test_blocked_components(self, model, media, demands = [], sinks = []):
        self.blocked_demands = []
        self.available_demands = []
        self.blocked_rxns = []

        with model as model:
            for met in media:
                model.add_boundary(model.metabolites.get_by_id(met), type = 'exchange')
            model.medium = dict([[rxn.id, 100] for rxn in model.exchanges])
            print(model.medium)

            if sinks:
                for met in sinks:
                    model.add_boundary(model.metabolites.get_by_id(met), type = 'sink')

            if demands:
                demands_to_test = demands
                self.run_blocked_demand(model, demands_to_test)
            else:
                demands_to_test = self.biomass_in_model.loc[:, 'met_id'].to_list()
                self.run_blocked_demand(model, demands_to_test)
            return self.blocked_demands, self.available_demands, self.blocked_rxns

    def fastcc(self, model, flux_threshold = 1.0, zero_cutoff = None):
        zero_cutoff = normalize_cutoff(model, zero_cutoff)

        irreversible_rxns = [rxn for rxn in model.reactions if not rxn.reversibility]
        rxns_to_check = irreversible_rxns

        with model:
            rxns_to_keep = _find_sparse_mode(
                model, rxns_to_check, flux_threshold, zero_cutoff
                )

        rxns_to_check = list(set(model.reactions).difference(rxns_to_keep))

        while rxns_to_check:
            with model:
                new_rxns = _find_sparse_mode(
                    model, rxns_to_check, flux_threshold, zero_cutoff
                    )
                rxns_to_keep.extend(new_rxns)

                # this condition will be valid for all but the last iteration
                if list(set(rxns_to_check).intersection(rxns_to_keep)):
                    rxns_to_check = list(set(rxns_to_check).difference(rxns_to_keep))
                else:
                    rxns_to_flip = list(set(rxns_to_check).difference(irreversible_rxns))
                    _flip_coefficients(model, rxns_to_flip)
                    sol = model.optimize(min)
                    to_add_rxns = sol.fluxes.index[sol.fluxes.abs() > zero_cutoff].tolist()
                    rxns_to_keep.extend(
                        [model.reactions.get_by_id(rxn) for rxn in to_add_rxns]
                        )
                    # since this is the last iteration, it needs to break or else
                    # it will run forever since rxns_to_check won't be empty
                    break

        consistent_rxns = set(rxns_to_keep)
        # need the ids since Reaction objects are created fresh with model.copy()
        rxns_to_remove = [
            rxn.id for rxn in set(model.reactions).difference(consistent_rxns)
            ]

        return consistent_rxns, rxns_to_remove

    def rxns_of_metabolite(self, model, metabolite = '', in_solution = False):
        for rxn in model.reactions:
            met_in_rxn = [metabolite in met.id for met in rxn.metabolites]
            if any(met_in_rxn):
                if in_solution:
                    if rxn.id in self.solution.fluxes[((self.solution.fluxes > 0) | (self.solution.fluxes < 0))].index:
                        print(rxn.subsystem, rxn.id, rxn.build_reaction_string(), self.solution.fluxes[rxn.id])
                else:
                    print(rxn.subsystem, rxn.id, rxn.build_reaction_string())

    def rxns_in_pthwy(self, model, pathway = '', in_solution = False):
        for rxn in model.reactions:
            if in_solution:
                if pathway in rxn.subsystem and rxn.id in self.solution.fluxes[((self.solution.fluxes > 0) | (self.solution.fluxes < 0))].index:
                    print(rxn.id, rxn.build_reaction_string(), self.solution.fluxes[rxn.id])
            else:
                if pathway in rxn.subsystem:
                    print(rxn.id, rxn.build_reaction_string())
