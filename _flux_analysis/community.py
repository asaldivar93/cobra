from __future__ import absolute_import
import os
import pandas as pd
import cobra
import warnings
import plotly.express as px
import seaborn as sns
import matplotlib.pyplot as plot
import random
warnings.filterwarnings("ignore")
os.chdir('/home/alexis/UAM/cobra/')

from numpy import (sign,
                   log10,
                   arange,
                   zeros,
                   diag,
                   cov,
                   dot,
                   mean)
from numpy.random import (uniform,
                          binomial)
from numpy.linalg import pinv
from skbio.stats.composition import clr
from scipy.stats import (multivariate_normal,
                         bernoulli)
from datetime import datetime
from cobra.flux_analysis import (micom,
                                 micom_tradeoff)
from cobra.flux_analysis.micom import (_build_micom_constraints,
                                       _build_micom_objective)
from tools.cobra_tools import add_boundaries_for_simulation
from tools.p_tools import build_genus_database
from tools.plotting import plot_heatmap
from tqdm import trange


class community():
    def __init__(
        self,
        name,
        home_path='/home/alexis/UAM/cobra/',
        taxonomy_path='data_files/silva_taxonomy.csv',
        counts_path='data_files/otu-table-corrected.tsv',
        results_path='results/'
    ):
        os.chdir(home_path)
        self.taxonomy = pd.read_csv(taxonomy_path, sep='\t')
        self.taxonomy.set_index('sequence_identifier', inplace = True)
        self.taxonomy.fillna('unknown', inplace = True)
        self.seq_counts = pd.read_csv(counts_path, sep = '\t', index_col = 'OTU_ID')
        self.seq_counts.loc['0281ce539a44c8b1b059ae85c15e3fab', 'My_20'] = 0  # Outlier Mesoaciditoga
        self.results_path = results_path

    def set_samples(self, sample_dict):
        self.samples = list(sample_dict.keys())
        self.taxa_database = build_genus_database(self.samples, self.seq_counts, self.taxonomy)

    def load_models(self, sample_dict, to_close, close_ex=True, micom=True):
        self.models = dict()
        self.micom = micom
        for sample in self.samples:
            community_model = cobra.Model()
            community_model.solver = 'gurobi'
            sample_filter = self.taxa_database.loc[:, 'sample'] == sample
            sample_db = self.taxa_database[sample_filter].set_index('id')
            for k in self.taxa_database.loc[sample_filter].index:
                new_model = cobra.io.read_sbml_model(self.taxa_database.loc[k, 'model_path'])
                community_model.merge(new_model)

            self.permanent_close = ['Ubiquinols', 'CO_A', 'FMN', '10_FORMYL_THF', 'FAD', 'GLUTATHIONE', '2_OCTAPRENYLPHENOL', 'METHYLENE_THF', 'PYRIDOXAL_PHOSPHATE', 'S_ADENOSYLMETHIONINE', 'THF', 'Menaquinols', 'PROTOHEME', 'Glycogens', 'Poly_Hydroxybutyrate']

            for ex in self.permanent_close:
                for rxn in community_model.reactions:
                    if 'EX' in rxn.id and ex in rxn.id:
                        rxn.lower_bound = 0
                        rxn.upper_bound = 0
                        community_model.remove_reactions([rxn])

            for taxa_k in sample_db.index:
                for rxn in community_model.reactions:
                    if 'EX' in rxn.id and taxa_k in rxn.id:
                        proton = ['PROTON' in met.id for met in rxn.metabolites]
                        if any(proton):
                            rxn.add_metabolites({'_PROTON_pe_' + taxa_k: -1, '_PROTON_cy_' + taxa_k: 1})

            for rxn in community_model.reactions:
                if 'BIOMASS' in rxn.id:
                    for met in rxn.metabolites:
                        if '_ATP_' in met.id:
                            current_atpdm = -rxn.get_coefficient(met.id)
                            new_atpdm = sample_dict[sample]
                            to_add = current_atpdm - new_atpdm
                            rxn.add_metabolites({met: to_add})
            self.to_close = to_close
            if close_ex:
                for rxn_id in self.to_close:
                    rxn = community_model.reactions.get_by_id(rxn_id)
                    rxn.lower_bound = -100
                    rxn.upper_bound = 0

            if micom:
                micom_model = _build_micom_constraints(
                    community_model, sample_db
                )

            self.models[sample] = micom_model
        del community_model

    def update_exchanges(self, sample, actives):
        for rxn_id in self.to_close:
            if actives[rxn_id] == 1:
                rxn = self.models[sample].reactions.get_by_id(rxn_id)
                rxn.lower_bound = -100
                rxn.upper_bound = 100
            elif actives[rxn_id] == 0:
                rxn = self.models[sample].reactions.get_by_id(rxn_id)
                rxn.lower_bound = -100
                rxn.upper_bound = 0

    def plot_micom_tradeoff(self, media, template='none'):
        non_zero = pd.DataFrame(
            columns = ['Sample', 'Trade off', 'Growth Rate']
        )
        for sample in self.samples:
            sample_filter = self.taxa_database.loc[:, 'sample'] == sample
            sample_db = self.taxa_database[sample_filter].set_index('id')
            with self.models[sample] as model:
                add_boundaries_for_simulation(
                    model, media=media, carbon_uptake=0.00625
                    )
                nonzero = micom_tradeoff(
                    model, sample_db
                    )
            nonzero.loc[:, 'Sample'] = sample
            non_zero = non_zero.append(nonzero)

        non_zero.to_csv(self.results_path + 'micom_tradeoff.csv', sep='/t')
        fig = px.box(
            non_zero,
            y='Growth Rate',
            x='Trade off',
            points='all',
            color='Sample',
            template=template
        )
        fig.write_html(self.results_path + 'plots/micom_tradeoff.html')
        fig.write_image(self.results_path + 'plots/micom_tradeoff.png')
        return fig

    def plot_yields(self, media, carbon_uptake, metabolites, trade_off, template='none'):
        yields_df = pd.DataFrame(
            columns = ['Sample', 'Reaction', 'Metabolite', 'Yield']
        )
        for sample in self.samples:
            sample_filter = self.taxa_database.loc[:, 'sample'] == sample
            sample_db = self.taxa_database[sample_filter].set_index('id')
            with self.models[sample] as model:
                add_boundaries_for_simulation(model, media=media, carbon_uptake=carbon_uptake)
                if self.micom:
                    micom_solution = micom(
                        model=[],
                        micom_model=model,
                        trade_off=trade_off,
                        taxa_db=sample_db
                    )[3]
                else:
                    micom_solution = micom(
                        model=model,
                        trade_off=trade_off,
                        taxa_db=sample_db
                    )[3]

            biomass = self._get_biomass_yield(sample, micom_solution)

            yields_df = yields_df.append(
                pd.DataFrame(
                    [[sample, 'Biomass', 'Biomass', biomass]],
                    columns = ['Sample', 'Reaction', 'Metabolite', 'Yield']
                )
            )
            for rxn_id in metabolites.keys():
                met_yield = round(abs(micom_solution[rxn_id] / carbon_uptake), 2)
                yields_df = yields_df.append(
                    pd.DataFrame(
                        [[sample, rxn_id, metabolites[rxn_id]['name'], met_yield]],
                        columns = ['Sample', 'Reaction', 'Metabolite', 'Yield']
                    )
                )

        fig = px.bar(
            yields_df,
            x='Metabolite',
            y='Yield',
            color='Sample',
            barmode='group',
            text='Yield',
            template=template
        )

        try:
            fig.write_html(self.results_path + 'plots/' + self.name + '/yield_plots.html')
            fig.write_image(self.results_path + 'plots/' + self.name + '/yield_plots.svg')
        except FileNotFoundError:
            os.mkdir(self.results_path + 'plots/' + self.name)
            fig.write_html(self.results_path + 'plots/' + self.name + '/yield_plots.html')
            fig.write_image(self.results_path + 'plots/' + self.name + '/yield_plots.svg')
        return fig, yields_df

    def plot_metoh_ex(self, media, carbon_uptake, template='none'):
        metoh_ex_df = pd.DataFrame(
            columns=['Sample', 'Trade off', 'Exchanged Metanol']
        )
        for sample in self.samples:
            sample_filter = self.taxa_database.loc[:, 'sample'] == sample
            sample_db = self.taxa_database[sample_filter].set_index('id')
            for trade_off in arange(0.1, 1, 0.1):
                with self.models[sample] as model:
                    add_boundaries_for_simulation(model, media=media, carbon_uptake=carbon_uptake)
                    if self.micom:
                        micom_solution = micom(
                            model=[],
                            micom_model=model,
                            trade_off=trade_off,
                            taxa_db=sample_db
                        )[3]
                    else:
                        micom_solution = micom(
                            model=model,
                            trade_off=trade_off,
                            taxa_db=sample_db
                        )[3]

                metoh_ex = self._get_metohex(self.models[sample], micom_solution)
                metoh_ex_df = metoh_ex_df.append(
                    pd.DataFrame(
                        [[sample, trade_off, metoh_ex]],
                        columns=['Sample', 'Trade off', 'Exchanged Metanol']
                    )
                )

        fig = px.line(
            metoh_ex_df,
            x='Trade off',
            y='Exchanged Metanol',
            color='Sample',
            template=template
        )

        try:
            fig.write_html(self.results_path + 'plots/' + self.name + '/ex_metoh.html')
            fig.write_image(self.results_path + 'plots/' + self.name + '/ex_metoh.svg')
        except FileNotFoundError:
            os.mkdir(self.results_path + 'plots/' + self.name)
            fig.write_html(self.results_path + 'plots/' + self.name + '/ex_metoh.html')
            fig.write_image(self.results_path + 'plots/' + self.name + '/ex_metoh.svg')

        try:
            metoh_ex_df.to_csv(self.results_path + 'exchanges/' + self.name + '/ex_metoh.html')
        except FileNotFoundError:
            os.mkdir(self.results_path + 'exchanges/' + self.name)
            metoh_ex_df.to_csv(self.results_path + 'exchanges/' + self.name + '/ex_metoh.html')

        return fig

    def plot_ex_fluxes(self, media, carbon_uptake, trade_off):
        for sample in self.samples:
            external_metabolites = self._parse_external_metabolites(sample)
            sample_filter = self.taxa_database.loc[:, 'sample'] == sample
            sample_db = self.taxa_database[sample_filter].set_index('id')
            with self.models[sample] as model:
                add_boundaries_for_simulation(
                    model, media=media, carbon_uptake=carbon_uptake
                    )
                if self.micom:
                    micom_solution = micom(
                        model=[],
                        micom_model=model,
                        trade_off=trade_off,
                        taxa_db=sample_db
                    )[3]
                else:
                    micom_solution = micom(
                        model=model,
                        trade_off=trade_off,
                        taxa_db=sample_db
                    )[3]

            import_fluxes, export_fluxes = self._exfluxes_by_cmol(
                sample, external_metabolites, sample_db, micom_solution
                )[:-1]

            min_I = 7.787481e-12
            min_II = 7.787481e-12
            import_fluxes.fillna(0, inplace=True)
            import_fluxes = import_fluxes + 2 * min_I
            export_fluxes.fillna(0, inplace=True)
            export_fluxes = export_fluxes + 2 * min_II

            import_clr = import_fluxes.copy()
            export_clr = export_fluxes.copy()
            for m in import_clr.index:
                import_clr.loc[m, :] = clr(import_clr.loc[m, :])
            for m in export_clr.index:
                export_clr.loc[m, :] = clr(export_clr.loc[m, :])

            import_fluxes = log10(import_fluxes)
            export_fluxes = log10(export_fluxes)
            exp = import_fluxes - export_fluxes
            import_fluxes.reset_index(inplace=True)
            export_fluxes.reset_index(inplace=True)

            try:
                import_clr.to_csv(self.results_path + 'exchanges/' + self.name + '/import_clr_{}.csv'.format(sample))
                export_clr.to_csv(self.results_path + 'exchanges/' + self.name + '/export_clr_{}.csv'.format(sample))
                exp.to_csv(self.results_path + 'exchanges/' + self.name + '/exchanges_log10_{}.csv'.format(sample))
            except FileNotFoundError:
                os.mkdir(self.results_path + 'exchanges/' + self.name)
                import_clr.to_csv(self.results_path + 'exchanges/' + self.name + '/import_clr_{}.csv'.format(sample))
                export_clr.to_csv(self.results_path + 'exchanges/' + self.name + '/export_clr_{}.csv'.format(sample))
                exp.to_csv(self.results_path + 'exchanges/' + self.name + '/exchanges_log10_{}.csv'.format(sample))

            fig1 = plot_heatmap(
                exp, cmap = sns.diverging_palette(271, 130, s = 80, l = 47, sep = 20, as_cmap = True)
                )

            fig2 = plot_heatmap(
                import_clr, cmap = sns.diverging_palette(271, 130, s = 80, l = 47, sep = 20, as_cmap = True)
                )

            fig3 = plot_heatmap(
                export_clr, cmap = sns.diverging_palette(271, 130, s = 80, l = 47, sep = 20, as_cmap = True)
                )
            try:
                fig1.savefig(self.results_path + 'plots/' + self.name + '/export_clr_{}.svg'.format(sample))
                fig2.savefig(self.results_path + 'plots/' + self.name + '/exchanges_log10_{}.svg'.format(sample))
                fig3.savefig(self.results_path + 'plots/' + self.name + '/import_clr_{}.svg'.format(sample))
            except FileNotFoundError:
                os.mkdir(self.results_path + 'exchanges/' + self.name)
                fig1.savefig(self.results_path + 'plots/' + self.name + '/export_clr_{}.svg'.format(sample))
                fig2.savefig(self.results_path + 'plots/' + self.name + '/exchanges_log10_{}.svg'.format(sample))
                fig3.savefig(self.results_path + 'plots/' + self.name + '/import_clr_{}.svg'.format(sample))

    def plot_interaction_matrix(self, size=(9, 4)):
        self.interaction_matrix = dict()
        for sample in self.samples:
            external_fluxes = pd.read_csv(self.results_path + 'exchanges/' + self.name + '/exchanges_{}.csv'.format(sample), index_col=0)
            external_fluxes.drop(index='CH4', inplace=True)
            self.interaction_matrix[sample] = self._calculate_interaction_matrix(sample, external_fluxes)
            try:
                self.interaction_matrix[sample].to_csv(self.results_path + 'interaction_matrix/' + self.name + '/ic_{}.csv'.format(sample))
            except FileNotFoundError:
                os.mkdir(self.results_path + 'interaction_matrix/' + self.name)
                self.interaction_matrix[sample].to_csv(self.results_path + 'interaction_matrix/' + self.name + '/ic_{}.csv'.format(sample))

        fig, (ax1, ax2) = plot.subplots(1, 2)
        sns.heatmap(
            self.interaction_matrix['CIR_19'],
            cmap = sns.diverging_palette(41, 200, s=100, l=45, sep=1, center='light', as_cmap=True),
            vmin=-1,
            vmax=1,
            ax=ax1,
            xticklabels=True,
            cbar=False)
        sns.heatmap(
            self.interaction_matrix['My_20'],
            cmap = sns.diverging_palette(41, 200, s=100, l=45, sep=1, center='light', as_cmap=True),
            vmin=-1,
            vmax=1,
            xticklabels=True,
            ax=ax2)
        fig.set_size_inches(size, forward=True)

        try:
            fig.savefig(self.results_path + 'plots/' + self.name + '/ic_matrix.svg')
        except FileNotFoundError:
            os.mkdir(self.results_path + 'plots/' + self.name)
            fig.savefig(self.results_path + 'plots/' + self.name + '/ic_matrix.svg')

        return fig

    def _exfluxes_by_cmol(self, sample, external_metabolites, sample_db, solution):
        import_fluxes = pd.DataFrame(index = external_metabolites, columns = sample_db.index)
        export_fluxes = pd.DataFrame(index = external_metabolites, columns = sample_db.index)
        external_fluxes = pd.DataFrame(index = external_metabolites, columns = sample_db.index)
        for rxn in self.models[sample].reactions:
            if '_EX' in rxn.id:
                metabolite = rxn.id[4:-5]
                genus_k = rxn.id[-4:]
                abundance_k = sample_db.loc[genus_k, 'relative_abundance']
                if metabolite in import_fluxes.index:
                    if solution[rxn.id] > 0:
                        import_fluxes.loc[metabolite, genus_k] = 1000 * solution[rxn.id] * abundance_k * list(rxn.metabolites)[0].cmol_by_mol
                    elif solution[rxn.id] < 0:
                        export_fluxes.loc[metabolite, genus_k] = -1000 * solution[rxn.id] * abundance_k * list(rxn.metabolites)[0].cmol_by_mol
                    external_fluxes.loc[metabolite, genus_k] = 1000 * solution[rxn.id] * abundance_k * list(rxn.metabolites)[0].cmol_by_mol
        external_fluxes.fillna(0, inplace = True)
        try:
            external_fluxes.to_csv(self.results_path + 'exchanges/' + self.name + '/exchanges_{}.csv'.format(sample))
        except FileNotFoundError:
            os.mkdir(self.results_path + 'exchanges/' + self.name)
            external_fluxes.to_csv(self.results_path + 'exchanges/' + self.name + '/exchanges_{}.csv'.format(sample))

        return import_fluxes, export_fluxes, external_fluxes

    def _parse_external_metabolites(self, sample):
        ex_mets = []
        # Filter external reactions
        for rxn in self.models[sample].reactions:
            if 'EX' in rxn.id:
                # Save only the name of themetabolite
                ex_mets.append(rxn.id[4:-5])
        # Remove  sitem inputs (CH4, O2, NO3, etc.)
        external_metabolites = list(pd.Series(ex_mets).unique()[0:2])
        b = list(pd.Series(ex_mets).unique()[13:])
        external_metabolites.extend(b)

        return external_metabolites

    def _calculate_interaction_matrix(self, sample, external_fluxes):
        interaction_matrix = pd.DataFrame(index = external_fluxes.columns, columns = external_fluxes.columns, dtype = 'float')
        for k in interaction_matrix.index:
            imports = external_fluxes[k] > 0
            total_imports_k = external_fluxes.loc[imports, k].sum()
            for j in interaction_matrix.columns:
                flux_from_j = -external_fluxes.loc[imports, j].sum()
                total_external = external_fluxes.loc[imports].drop(columns = [j, k], axis=1).sum()
                imports_others = -total_external[total_external > 0].sum()
                exports_others = -total_external[total_external < 0].sum()
                if k == j:
                    IC_j_in_k = 0
                elif sign(flux_from_j) == sign(total_imports_k):
                    IC_j_in_k = (total_imports_k - (total_imports_k - flux_from_j)) / (flux_from_j + exports_others)
                else:
                    # IC_j_in_k = -(exports_others - (exports_others + flux_from_j)) / (exports_others + imports_others)
                    IC_j_in_k = flux_from_j / (exports_others + imports_others)
                interaction_matrix.loc[k, j] = IC_j_in_k

        return interaction_matrix

    def _get_metohex(self, model, solution):
        summary = model.metabolites._METOH_pe_Meis.summary(solution=solution)
        try:
            metoh_ex = summary.consuming_flux.loc['_EX_METOH_Meis', 'percent']
        except KeyError as e:
            print(e)
            metoh_ex = 0

        return metoh_ex

    def _get_biomass_yield(self, sample, solution):
        biomass_fluxes = solution.fluxes.index.str.contains('biomass')
        biomass = 0
        for tax in solution.fluxes.loc[biomass_fluxes].index:
            biomass += solution[tax] * self.models[sample].metabolites.get_by_id(tax.replace('DM_', '')).cmol_by_mol

        return abs(round(biomass / solution['EX__CH4_ou'], 3))

    def _melt_ic_matrix(self, ic_matrix):
        columns = ['Giver', 'Reciver', 'Interaction Coefficient']
        ic_long = pd.DataFrame(
            columns=columns
        )
        for giver in ic_matrix.columns:
            for reciver in ic_matrix.index:
                ic_long = ic_long.append(
                    pd.DataFrame(
                        [[giver, reciver, ic_matrix.loc[reciver, giver]]],
                        columns=columns
                    )
                )

        return ic_long

    def fit_model(
        self,
        sample,
        media,
        carbon_uptake,
        trade_off,
        yield_mets,
        n_sim=100,
        warmup=50,
        permute_init=False,
        dir=[]
    ):

        external_metabolites = self._parse_external_metabolites(sample)

        sample_filter = self.taxa_database.loc[:, 'sample'] == sample
        sample_db = self.taxa_database[sample_filter].set_index('id')
        micom_model = _build_micom_constraints(
            self.models[sample], sample_db
        )
        micom_objective = _build_micom_objective(
            micom_model, sample_db
        )
        add_boundaries_for_simulation(
            micom_model, media=media, carbon_uptake=carbon_uptake
            )

        gh = gibbs(
            self.results_path,
            micom_objective,
            trade_off,
            sample_db,
            yield_mets,
            dir=dir)

        n_ex = len(self.to_close)  # number of external metabolites to fit
        rxns_to_fit = pd.Series(self.to_close)  # reactions to fit

        if dir:
            theta_sampler = pd.read_csv(self.results_path + dir + '/theta_samples.csv', index_col='Unnamed: 0')
            columns = pd.Series(theta_sampler.columns.to_list(), dtype='int64')
            theta_sampler.columns = columns
            iter_k = theta_sampler.columns[-1]
        else:
            iter_k = 0
            # Initialize binary matrix
            # set initial vector to be all zeros
            theta_sampler = pd.DataFrame(
                zeros((n_ex, n_sim)), index=self.to_close, columns=range(n_sim), dtype='int64'
                )
            # Set initial vector
            if permute_init:
                for i in range(n_ex):
                    # initial theta set to be uniformly distributed
                    theta_sampler.iloc[i, 0] = int(binomial(size=1, n=1, p=uniform(0, 1))[0])

        for k in trange(iter_k, n_sim - 1):
            theta_old = theta_sampler.loc[:, k].copy()
            n_changes = random.randint(1, n_ex)  # number of changes to make in iter i
            index = random.sample(range(0, n_ex), n_changes)  # indexes of the reactions to be changed
            rxns_in_iter = rxns_to_fit.iloc[index]  # reactions to change in this interation

            # copy old proposition
            theta_new = theta_old.copy()
            # generate new propositions for rxns in rxns_in_iter
            if k <= warmup:
                # metropolized gibbs
                for rxn in rxns_in_iter:
                    m_i = 1 - theta_old[rxn]
                    theta_new.loc[rxn] = int(binomial(size=1, n=1, p=m_i))
            else:
                # adaptive gibbs
                W = pd.DataFrame(
                    pinv(cov(theta_sampler.loc[:, :k])), index=rxns_to_fit, columns=rxns_to_fit
                )  # Covariance matrix of theta i in k iterations
                for rxn in rxns_in_iter:
                    mu = theta_sampler.loc[rxn, :k].mean()
                    m_i = gh._m_theta_(rxn, mu, W, theta_old)
                    theta_new.loc[rxn] = int(binomial(size=1, n=1, p=m_i))

            # estimate conditional probabilities
            p_data_new = gh._prob_of_simulation(sample, micom_model, theta_new)[0]
            p_data_old, sol_old = gh._prob_of_simulation(sample, micom_model, theta_old)

            if k <= warmup:
                # metropolized gibbs
                r = p_data_new / p_data_old  # calculate acceptance rate
                # update solution
                if r >= uniform(0, 1):
                    theta_sampler.loc[:, k + 1] = theta_new.copy()
                else:
                    theta_sampler.loc[:, k + 1] = theta_old.copy()
            else:
                for rxn in rxns_in_iter:
                    mu_old = theta_sampler.loc[rxn, :k].mean()
                    mu_new = list(theta_sampler.loc[rxn, :k])
                    mu_new.append(theta_new.loc[rxn])
                    mu_new = mean(mu_new)
                    m_old = gh._m_theta_(rxn, mu_old, W, theta_old)
                    m_new = gh._m_theta_(rxn, mu_new, W, theta_new)
                    p_theta_new = bernoulli.pmf(theta_new.loc[rxn], m_old)  # p of new given old
                    p_theta_old = bernoulli.pmf(theta_old.loc[rxn], m_new)  # p of old given new
                    # calculate acceptance rate
                    r = (p_data_new * p_theta_old) / (p_data_old * p_theta_new)
                    # update solution
                    if not r >= uniform(0, 1):
                        theta_new[rxn] = theta_old[rxn]
                theta_sampler.loc[:, k + 1] = theta_new.copy()

            # summarize external fluxes by c_mol
            external_fluxes = self._exfluxes_by_cmol(
                sample, external_metabolites, sample_db, sol_old
                )[-1]
            # Get variables
            interaction_matrix = self._calculate_interaction_matrix(
                sample, external_fluxes
                )
            ic_long = self._melt_ic_matrix(
                interaction_matrix
                )
            metoh_ex = self._get_metohex(
                micom_model, sol_old
            )
            biomass = self._get_biomass_yield(
                sample, sol_old
                )
            yields_hat = {rxn: round(abs(sol_old[rxn] / sol_old['EX__CH4_ou']), 3) for rxn in yield_mets.keys()}
            n_open = theta_sampler.loc[:, k].sum()

            gh._update_files(
                sample, n_open, ic_long, metoh_ex, biomass, yields_hat
                )
            theta_sampler.loc[:, :k].to_csv(self.results_path + gh.dir + '/theta_samples.csv')


class gibbs():
    def __init__(
        self,
        results_path,
        micom_objective,
        trade_off,
        sample_db,
        yield_mets,
        dir
    ):
        self.results_path = results_path
        self.micom_objective = micom_objective
        self.trade_off = trade_off
        self.sample_db = sample_db
        self.yield_mets = yield_mets

        if not dir:
            now = datetime.now()
            dt_string = now.strftime("%d-%m-%Y_%H-%M-%S")
            self.dir = 'fit_' + dt_string
            os.mkdir(self.results_path + self.dir)
            self.ic_filename = self.dir + '/interactions.csv'
            with open(self.results_path + self.ic_filename, 'w+') as file:
                pd.DataFrame(columns=['Giver', 'Reciver', 'Interaction Coefficient', 'Sample', 'N Open']).to_csv(file, index=False)

            self.metoh_filename = self.dir + '/ex_metoh.csv.csv'
            with open(self.results_path + self.metoh_filename, 'w+') as file:
                pd.DataFrame(columns=['Sample', 'Exchanged Metanol', 'N Open']).to_csv(file, index=False)

            self.yields_filename = self.dir + '/yields.csv'
            with open(self.results_path + self.yields_filename, 'w+') as file:
                pd.DataFrame(columns=['Sample', 'Metabolite', 'Yield', 'N Open']).to_csv(file, index=False)
        else:
            self.dir = dir
            self.ic_filename = self.dir + '/interactions.csv'
            self.metoh_filename = self.dir + '/ex_metoh.csv.csv'
            self.yields_filename = self.dir + '/yields.csv'

    def _read_means(self, sample):
        with open(self.results_path + self.yields_filename, 'r') as file:
            df = pd.read_csv(file)
            df = df.query("Sample==@sample")
            mu_new = dict()
            for rxn_id in self.yield_mets.keys():
                met = self.yield_mets[rxn_id]['name']
                mu = df.query("Metabolite==@met")['Yield'].mean()
                mu_new[rxn_id] = mu
            iter_n = len(df.query("Metabolite==@met")['Yield']) - 1

        return [mu_new, iter_n]

    def _prob_of_simulation(self, sample, micom_model, theta_vector):
        with micom_model as model:
            for i in theta_vector.index:
                if theta_vector.loc[i] == 1:
                    rxn = micom_model.reactions.get_by_id(i)
                    rxn.lower_bound = -100
                    rxn.upper_bound = 100
                else:
                    rxn = micom_model.reactions.get_by_id(i)
                    rxn.lower_bound = -100
                    rxn.upper_bound = 0

            micom_model, linear_solution, quadratic_solution, micom_solution, taxa_growth_rates = micom(
                micom_model=model,
                micom_objective=self.micom_objective,
                trade_off=self.trade_off,
                taxa_db=self.sample_db,
                model=[]
                )

        new_y_hat = {rxn_id: abs(micom_solution[rxn_id] / micom_solution['EX__CH4_ou']) for rxn_id in self.yield_mets.keys()}
        y_obsv = [self.yield_mets[rxn_id]['yields'][sample] for rxn_id in self.yield_mets.keys()]
        s_obsv = [self.yield_mets[rxn_id]['std'][sample]**2 for rxn_id in self.yield_mets.keys()]
        mvn = multivariate_normal(mean=y_obsv, cov=diag(s_obsv))
        cdf_upper = mvn.cdf([new_y_hat[rxn_id] + 0.5 * self.yield_mets[rxn_id]['std'][sample] for rxn_id in new_y_hat.keys()])
        cdf_lower = mvn.cdf([new_y_hat[rxn_id] - 0.5 * self.yield_mets[rxn_id]['std'][sample] for rxn_id in new_y_hat.keys()])
        probability = cdf_upper - cdf_lower

        return probability, micom_solution

    def _m_theta_(self, rxn, mean, W, theta_vector):
        delta = 0.1
        f_theta = mean - dot(W.loc[:, rxn].drop(rxn), theta_vector.drop(rxn)) / W.loc[rxn, rxn]
        m_i = min(max(f_theta, delta), 1 - delta)

        return m_i

    def _update_files(self, sample, n_open, ic_long, metoh_ex, biomass, yields_hat):
        ic_long.loc[:, 'Sample'] = sample
        ic_long.loc[:, 'N Open'] = n_open
        with open(self.results_path + self.ic_filename, 'a+') as file:
            ic_long.to_csv(file, index=False, header=False)

        with open(self.results_path + self.metoh_filename, 'a+') as file:
            pd.DataFrame(
                [[sample, metoh_ex, n_open]],
                columns=['Sample', 'Exchanged Metanol', 'N Open']
                ).to_csv(file, index=False, header=False)

        with open(self.results_path + self.yields_filename, 'a+') as file:
            pd.DataFrame(
                [[sample, 'Biomass', biomass, n_open]],
                columns=['Sample', 'Metabolite', 'Yield', 'N Open']
                ).to_csv(file, index=False, header=False)

        for rxn_id in yields_hat.keys():
            with open(self.results_path + self.yields_filename, 'a+') as file:
                pd.DataFrame(
                    [[sample, self.yield_mets[rxn_id]['name'], yields_hat[rxn_id], n_open]],
                    columns=['Sample', 'Metabolite', 'Yield', 'N Open']
                    ).to_csv(file, index=False, header=False)
