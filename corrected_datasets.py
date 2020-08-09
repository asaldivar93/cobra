import pandas as pd

eqtr_databases = {'bigg.metabolite',
                  'chebi',
                  'envipath',
                  'hmdb',
                  'kegg',
                  'metacyc.compound',
                  'metanetx.chemical',
                  'reactome',
                  'sabiork.compound',
                  'seed',
                  'synonyms'}

eqtr_metacyc_map = pd.DataFrame([('bigg.metabolite', '|BIGG|'),
                                 ('chebi', '|CHEBI|'),
                                 ('hmdb', '|HMDB|'),
                                 ('kegg', '|KEGG|'),
                                 ('metanetx.chemical', '|METANETX|'),
                                 ('seed', '|SEED|')], columns=['equilibrator', 'metacyc'])

corrected_metabolites = {'|Donor-H2|': {'formula': 'H2'}}

corrected_stoichiometry = {
    '|RXN-9364|': {'stoichiometry': {"|FARNESYL-PP|_cy": -1, "|DELTA3-ISOPENTENYL-PP|_cy": -10, "|CPD-9972|_cy": 1, "|PPI|_cy": 10}, 'balanced': []},
    '|SUCCINATE-SEMIALDEHYDE-DEHYDROGENASE-RXN|': {'stoichiometry': {'|NAD|_cy': -1, '|SUCC-S-ALD|_cy': -1, '|WATER|_cy': -1, '|NADH|_cy': 1, '|PROTON|_cy': 2, '|SUC|_cy': 1}, 'balanced': []},
    '|RXN-9138|': {'stoichiometry': {'|DELTA3-ISOPENTENYL-PP|_cy': -8, '|FARNESYL-PP|_cy': -1, '|CPD-9649|_cy': 1, '|PPI|_cy': 8}, 'balanced': []},
    '|TRANS-HEXAPRENYLTRANSTRANSFERASE-RXN|': {'stoichiometry': {'|DELTA3-ISOPENTENYL-PP|_cy': -4, '|FARNESYL-PP|_cy': -1, '|ALL-TRANS-HEPTAPRENYL-DIPHOSPHATE|_cy': 1, '|PPI|_cy': 4}, 'balanced': []},
    '|RXN-17122|': {'stoichiometry': {'|5-PHOSPHORIBOSYL-5-AMINOIMIDAZOLE|_cy': -1, '|Donor-H2|_cy': -1, '|S-ADENOSYLMETHIONINE|_cy': -1, '|AMMONIUM|_cy': 1, '|Acceptor|_cy': 1, '|CH33ADO|_cy': 1, '|CPD-18497|_cy': 1, '|FORMATE|_cy': 1, '|MET|_cy': 1, '|PROTON|_cy': 2, '|Pi|_cy': 1}, 'balanced': []},
    '|1.13.11.55-RXN|': {'stoichiometry': {'|Elemental-Sulfur|_cy': -4, '|OXYGEN-MOLECULE|_cy': -1, '|WATER|_cy': -4, '|HS|_cy': 2, '|PROTON|_cy': 4, '|SO3|_cy': 2}, 'balanced': []},
    '|ADENYLYLSULFATE-REDUCTASE-RXN|': {'stoichiometry': {'|AMP|_cy': -1, '|Acceptor|_cy': -1, '|PROTON|_cy': -2, '|SO3|_cy': -1, '|APS|_cy': 1, '|Donor-H2|_cy': 1}, 'balanced': []},
    '|SULFITE-DEHYDROGENASE-RXN|': {'stoichiometry': {'|Cytochromes-C-Oxidized|_cy': -2, '|SO3|_cy': -1, '|WATER|_cy': -1, '|Cytochromes-C-Reduced|_cy': 2, '|PROTON|_cy': 2, '|SULFATE|_cy': 1}, 'balanced': []},
    '|1.8.5.2-RXN|': {'stoichiometry': {'|Quinones|_cy': -1, '|S2O3|_cy': -2, '|CPD-14|_cy': 1, '|Reduced-Quinones|_cy': 1}, 'balanced': []},
    '|RXN-16432|': {'stoichiometry': {'|Acceptor|_cy': -1, '|CPD-14|_cy': -1, '|DsrE3A-L-cysteine|_cy': -2, '|Donor-H2|_cy': 1, '|DsrE3A-L-cysteine-S-thiosulfonates|_cy': 2}, 'balanced': []},
    '|RXN-17802|': {'stoichiometry': {'|Acceptor|_cy': -1, '|TusA-L-cysteine-S-thiosulfonates|_cy': -1, '|WATER|_cy': -3, '|Donor-H2|_cy': 1, '|PROTON|_cy': 3, '|SO3|_cy': 2, '|TusA-L-cysteine|_cy': 1}, 'balanced': []},
    '|CARDIOLIPSYN-RXN|': {'stoichiometry': {'|L-1-PHOSPHATIDYL-GLYCEROL|_cy': -2, '|CARDIOLIPIN|_cy': 1, '|GLYCEROL|_cy': 1}, 'balanced': []},
    '|PHOSPHASERSYN-RXN|': {'stoichiometry': {'|CDPDIACYLGLYCEROL|_cy': -1, '|SER|_cy': -1, '|CMP|_cy': 1, '|L-1-PHOSPHATIDYL-SERINE|_cy': 1, '|PROTON|_cy': 1}, 'balanced': True},
    '': {'stoichiometry': {}, 'balanced': []}
}
