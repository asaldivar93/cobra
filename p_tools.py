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


# %%codecell
