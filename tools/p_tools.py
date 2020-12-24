import pandas as pd
import os


class picrust_tools():
    def __init__(self, picrust_dir = '/home/alexis/UAM/picrust/Mar20/OTUs_99/', results_dir = 'pathways_out/'):
        os.chdir(picrust_dir)
        self.picrust_dir = picrust_dir
        self.results_dir = results_dir
        os.chdir(self.picrust_dir + self.results_dir)

    def change_p_dir(self, new_p_dir):
        self.picrust_dir = new_p_dir
        os.chdir(self.picrust_dir + self.results_dir)

    def change_r_dir(self, new_r_dir):
        self.results_dir = new_r_dir
        os.chdir(self.picrust_dir + self.results_dir)

    def sample_strat_pathway(self, pathways_path = 'path_abun_contrib.tsv.gz'):
        pthwy_strat_df = pd.read_csv(
            pathways_path, compression = 'gzip', sep = '\t'
            )

        samples = pthwy_strat_df.loc[:, 'sample'].unique()
        sample_strat = {}
        for sample in samples:
            sample_filter = pthwy_strat_df.loc[:, 'sample'] == sample
            pwy_sample_df = pthwy_strat_df.loc[sample_filter, :]

            #   We construct a tax_func sparse matrix where rows will be functions and columns will be taxons
            otus = pwy_sample_df['taxon'].unique()
            functions = pwy_sample_df['function'].unique()
            sample_strat[sample] = pd.DataFrame(
                index = functions, columns = otus
                )

            # iterate over all taxons
            for otu_id in otus:
                # Create a dataset with the information of all functions but only one taxon
                tax_filter = pwy_sample_df.loc[:, 'taxon'] == otu_id
                temp_df = pwy_sample_df.loc[tax_filter, ['function', 'taxon_rel_function_abun']]
                temp_df.set_index(
                    'function', inplace = True
                    )
                # iterate for all functions to which taxonid contributes
                for pwy_id in temp_df.index:
                    # replace the entry for funid in relabun_matrix with taxon_rel_function_abun of the filtered datased for said funid
                    sample_strat[sample].loc[pwy_id, otu_id] = temp_df.loc[pwy_id, 'taxon_rel_function_abun']
                # delete the filtered dataset
                del temp_df

        return sample_strat


def build_genus_database(samples, seq_counts, taxonomy):
    otus = pd.DataFrame()
    for sample in samples:
        sample_counts = seq_counts[sample].loc[seq_counts[sample] != 0].to_frame().rename(columns = {sample: 'seq_counts'}).reset_index()
        sample_counts.loc[:, 'sample'] = sample
        otus = otus.append(sample_counts)

    for seq in otus['OTU_ID'].unique():
        seq_filter = otus.loc[:, 'OTU_ID'] == seq
        otus.loc[seq_filter, 'Family'] = taxonomy.loc[seq, 'Family']
        otus.loc[seq_filter, 'Genus'] = taxonomy.loc[seq, 'Genus']

    taxa_database = pd.DataFrame(columns = ['id', 'seq_counts', 'sample', 'Genus', 'model_path'])
    for sample in samples:
        sample_filter = otus.loc[:, 'sample'] == sample
        new_otus = otus.loc[sample_filter, :].set_index('OTU_ID')
        for genus in new_otus.loc[:, 'Genus'].unique():
            if genus == 'uncultured' or genus == 'unknown':
                gen_filter = new_otus.loc[new_otus['Genus'] == genus].index
                families = new_otus.loc[gen_filter, :]
                for family in families['Family'].unique():
                    gen_filter = families.loc[families['Family'] == family].index
                    gen_id = family[:2] + family[-2:]
                    model_path = 'models/models_gapfilled/{}.xml'.format(gen_id)
                    taxa_database = taxa_database.append(pd.DataFrame(
                        [[gen_id,
                          families.loc[gen_filter, 'seq_counts'].sum(),
                          sample,
                          family,
                          model_path]],
                        columns = ['id', 'seq_counts', 'sample', 'Genus', 'model_path']
                        )
                    )

            else:
                gen_filter = new_otus['Genus'] == genus
                gen_id = genus[:2] + genus[-2:]
                model_path = 'models/models_gapfilled/{}.xml'.format(gen_id)
                taxa_database = taxa_database.append(pd.DataFrame(
                    [[gen_id,
                      new_otus.loc[gen_filter, 'seq_counts'].sum(),
                      sample,
                      genus,
                      model_path]],
                    columns = ['id', 'seq_counts', 'sample', 'Genus', 'model_path']
                    )
                )

    taxa_database.reset_index(inplace = True)
    for sample in samples:
        sample_filter = taxa_database.loc[:, 'sample'] == sample
        for k in taxa_database.loc[sample_filter, :].index:
            taxa_database.loc[k, 'relative_abundance'] = (taxa_database.loc[k, 'seq_counts'] / taxa_database.loc[sample_filter, 'seq_counts'].sum())

    return taxa_database
