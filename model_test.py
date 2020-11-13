import cobra
import picrust_parser as pparser
from data_files.corrected_datasets import corrected_revesibility
pc = pparser.p_model()

# %%codecell
del model
model = c_model.copy()
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

biomass_in_model = pc.search_biomass_components(model)
model = pc.add_biomass_rxn(model, biomass_in_model)
model.metabolites.get_by_id('biomass').formula_weight


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

with model as model:
    for substrate in media:
        met = model.metabolites.get_by_id(substrate)
        model.add_boundary(
            met, type = 'exchange'
            )
    medium = {rxn.id: 1000 for rxn in model.exchanges}
    medium['EX_|CH4|_ex'] = 100
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

    model.objective = {
        # model.reactions.get_by_id('DM_Intracellular_Metabolites'): 1.0
        model.reactions.get_by_id('BIOMASS'): 1.0
        }

    sol = cobra.flux_analysis.pfba(model)
    print(model.summary())
    print(model.metabolites.get_by_id('|OXYGEN-MOLECULE|_cy').summary())
