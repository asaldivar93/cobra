media = ['|CH4|_ex', '|OXYGEN-MOLECULE|_ex', '|NITRATE|_ex', '|FE+2|_ex', '|Pi|_ex', '|SULFATE|_ex']

sinks = ['|WATER|_pe', '|WATER|_cy', '|Menaquinones|_it', '|ACP|_cy',
         '|CPD-1302|_cy', '|CPD-1301|_cy', '|Guanine34-in-tRNAs|_cy', '|CL-|_cy',
         '|Hpr-Histidine|_cy', '|biotin-L-lysine-in-BCCP-dimers|_cy',
         '|carboxybiotin-L-lysine-in-BCCP-dimers|_cy', '|S-ADENOSYL-4-METHYLTHIO-2-OXOBUTANOATE|_cy',
         '|Unsulfurated-Sulfur-Acceptors|_cy', '|Corrinoid-Adenosyltransferases|_cy',
         '|LysW-C-Terminal-L-Glutamate|_cy', '|CPD-17931|_pe', '|Cytochromes-C-Oxidized|_pe',
         '|CoI-Corrinoid-Fe-S-proteins|_cy', '|Cysteine-Desulfurase-L-cysteine|_cy',
         '|D-alanine-carrier-protein|_cy', '|DsrE3A-L-cysteine|_cy', '|Ox-Thioredoxin|_cy',
         '|Sulfur-Carrier-Proteins-ThiI|_cy', '|Thi-S|_cy', '|4Fe-4S+1|_cy'
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
                   '|CPD-20903|_cy', '|TOLUENE|_cy', '|Starch|_ex'
                   ]

true_DM = ['|NA+|_pe', '|DTDP-RHAMNOSE|_cy', '|ACP|_pe', '|FORMAMIDE|_cy',
           '|UNDECAPRENYL-DIPHOSPHATE|_pe', '|CPD-10640|_cy', '|UNKNOWN|_cy',
           '|N-ACETYL-D-GLUCOSAMINE|_cy', '|3-5-ADP|_cy', '|CPD-15999|_cy',
           '|ETHANOL-AMINE|_cy', '|CPD-1091|_cy', '|ALLYSINE|_cy', '|Alcohols|_cy',
           '|1-AMINO-PROPAN-2-OL|_cy', '|P3I|_cy', '|Cysteine-Desulfurase-L-cysteine|_cy',
           '|GLYCOLALDEHYDE|_cy', '|4Fe-4S+2|_cy', '|HYDROGEN-PEROXIDE|_cy'
           ]

possible_product = ['|NITROGEN-MOLECULE|_pe', '|HYDROGEN-MOLECULE|_cy', '|ETOH|_cy',
                    '|CPD-10755|_cy', '|NITROGEN-MOLECULE|_cy', '|Methylketones|_cy',
                    '|CPD-347|_cy'
                    ]

posible_biomass = ['|PROTOHEME|_cy', '|ECTOINE|_cy', '|ADENOSYLCOBALAMIN|_cy',
                   '|PALMITYL-COA|_cy', '|LAUROYLCOA-CPD|_cy', '|CPD-9955|_cy',
                   '|CPD-9957|_cy', '|CPD-9958|_cy', '|URATE|_cy', '|CPD-9247|_pe',
                   '|STEAROYL-COA|_cy', '|CPD-9245|_cy', '|CPD-11994|_cy', '|C3|_cy',
                   '|CPD-12336|_cy', '|tRNAs-with-queuine|_cy',
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
                   '|CPD-12129|_cy', '|CARDIOLIPIN|_cy'
                   ]
