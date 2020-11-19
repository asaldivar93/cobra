import pandas as pd
import p_tools

# %% codecell
pt = p_tools.picrust_tools()
pwy_strat = pt.sample_strat_pathway()
pwy_strat['CIR_19'].dropna(how = 'any')
# %% codecell
taxonomy = pd.read_csv('/home/alexis/UAM/qiime/Mar20/OTUs_99/silva_taxonomy.csv', sep ='\t')
taxonomy.set_index('sequence_identifier', inplace = True)
seq_counts = pd.read_csv('/home/alexis/UAM/qiime/Mar20/OTUs_99/otu-table-corrected.tsv', sep = '\t', index_col = 'OTU_ID')
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
