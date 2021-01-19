


def simulate_randomwalk(self, sample, media, carbon_uptake, trade_off, yield_mets, n_sim=100):
    ex_mets = []
    for rxn in self.models[self.samples[0]].reactions:
        if 'EX' in rxn.id:
            ex_mets.append(rxn.id[4:-5])
    external_metabolites = list(pd.Series(ex_mets).unique()[0:2])
    b = list(pd.Series(ex_mets).unique()[13:])
    external_metabolites.extend(b)
    exp_means = dict(EX__OXYGEN_MOLECULE_ou = float(1.31), _DM_CARBON_DIOXIDE_ex = float(0.64))
    with open(self.results_path + 'yields.csv', 'r') as file:
        df = pd.read_csv(file)
        df = df.query("Sample==@sample")
        o_mean = df.query("Metabolite=='O2'")['Yield'].mean()
        c_mean = df.query("Metabolite=='CO2'")['Yield'].mean()
        o_std = df.query("Metabolite=='O2'")['Yield'].std()
        c_std = df.query("Metabolite=='CO2'")['Yield'].std()
        means = dict(EX__OXYGEN_MOLECULE_ou=o_mean, _DM_CARBON_DIOXIDE_ex=c_mean)
        stds = dict(EX__OXYGEN_MOLECULE_ou=o_std, _DM_CARBON_DIOXIDE_ex=c_std)
        k = len(df.query("Metabolite=='O2'")['Yield'])

    n_ex = len(self.to_close)
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

    n = 0
    while n < n_sim:
        max_del = random.randint(0, n_ex)
        index = random.sample(range(0, n_ex), max_del)

        with micom_model as model:
            for met in pd.Series(self.to_close).iloc[index]:
                rxn = model.reactions.get_by_id(met)
                rxn.lower_bound = -100
                rxn.upper_bound = 100

            micom_model, linear_solution, quadratic_solution, micom_solution, taxa_growth_rates = micom(
                micom_model=micom_model,
                micom_objective=micom_objective,
                trade_off=trade_off,
                taxa_db=sample_db,
                model=[]
                )

        joint_last = 1
        joint_new = 1
        for rxn_id in yield_mets.keys():
            met_yield = round(abs(micom_solution[rxn_id] / carbon_uptake), 2)
            p_met_last = norm.pdf(
                met_yield, loc=means[rxn_id], scale=stds[rxn_id]
                )
            joint_last = joint_last * p_met_last
            # p_met = norm.pdf(
            #     met_yield, loc=means[rxn_id] + (met_yield - means[rxn_id]) / k, scale=stds[rxn_id]
            #     )
            p_met = norm.pdf(
                met_yield, loc=exp_means[rxn_id], scale=0.1
                )
            joint_new = joint_new * p_met

        r = joint_new / joint_last
        if r > uniform(0, 1):
            n += 1
            print(n)
            external_fluxes = pd.DataFrame(index = external_metabolites, columns = sample_db.index)
            for rxn in micom_model.reactions:
                if '_EX' in rxn.id:
                    metabolite = rxn.id[4:-5]
                    genus_k = rxn.id[-4:]
                    abundance_k = sample_db.loc[genus_k, 'relative_abundance']
                    if metabolite in external_fluxes.index:
                        external_fluxes.loc[metabolite, genus_k] = 1000 * micom_solution[rxn.id] * abundance_k * list(rxn.metabolites)[0].cmol_by_mol

            external_fluxes.drop(index='CH4', inplace=True)
            external_fluxes.fillna(0, inplace = True)

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
            ic_long = self.ic_matrix_melt(interaction_matrix)
            ic_long.loc[:, 'Sample'] = sample
            ic_long.loc[:, 'N Blocked'] = max_del
            with open(self.results_path + 'simulations.csv', 'a+') as simulations:
                ic_long.to_csv(simulations, index=False, header=False)

            summary = micom_model.metabolites._METOH_pe_Meis.summary(solution=micom_solution)
            try:
                metoh_ex = summary.consuming_flux.loc['_EX_METOH_Meis', 'percent']
            except KeyError as e:
                print(e)
                metoh_ex = 0

            with open(self.results_path + 'metoh_ex.csv', 'a+') as file:
                pd.DataFrame(
                    [[sample, metoh_ex, max_del]],
                    columns=['Sample', 'Exchanged Metanol', 'N Blocked']
                    ).to_csv(file, index=False, header=False)

            biomass_fluxes = micom_solution.fluxes.index.str.contains('biomass')
            biomass = 0
            for tax in micom_solution.fluxes.loc[biomass_fluxes].index:
                biomass += micom_solution[tax] * self.models[sample].metabolites.get_by_id(tax.replace('DM_', '')).cmol_by_mol
            biomass = round(biomass / carbon_uptake, 2)

            with open(self.results_path + 'yields.csv', 'a+') as file:
                pd.DataFrame(
                    [[sample, 'Biomass', biomass, max_del]],
                    columns=['Sample', 'Metabolite', 'Yield', 'N Blocked']
                    ).to_csv(file, index=False, header=False)

            for rxn_id in yield_mets.keys():
                met_yield = round(abs(micom_solution[rxn_id] / carbon_uptake), 2)
                with open(self.results_path + 'yields.csv', 'a+') as file:
                    pd.DataFrame(
                        [[sample, yield_mets[rxn_id], met_yield, max_del]],
                        columns=['Sample', 'Metabolite', 'Yield', 'N Blocked']
                        ).to_csv(file, index=False, header=False)

            with open(self.results_path + 'yields.csv', 'r') as file:
                df = pd.read_csv(file)
                df = df.query("Sample==@sample")
                o_mean = df.query("Metabolite=='O2'")['Yield'].mean()
                c_mean = df.query("Metabolite=='CO2'")['Yield'].mean()
                o_std = df.query("Metabolite=='O2'")['Yield'].std()
                c_std = df.query("Metabolite=='CO2'")['Yield'].std()
                means = dict(EX__OXYGEN_MOLECULE_ou=o_mean, _DM_CARBON_DIOXIDE_ex=c_mean)
                stds = dict(EX__OXYGEN_MOLECULE_ou=o_std, _DM_CARBON_DIOXIDE_ex=c_std)
                k = len(df.query("Metabolite=='O2'")['Yield'])






def simulate_blocked_ex(self, sample, media, carbon_uptake, trade_off, yield_mets, n_sim=25):
    ex_mets = []
    for rxn in self.models[self.samples[0]].reactions:
        if 'EX' in rxn.id:
            ex_mets.append(rxn.id[4:-5])
    external_metabolites = list(pd.Series(ex_mets).unique()[0:2])
    b = list(pd.Series(ex_mets).unique()[13:])
    external_metabolites.extend(b)

    columns = ['Giver', 'Reciver', 'Interaction Coefficient', 'Sample', 'N Blocked']
    if not os.path.exists(self.results_path + 'simulations.csv'):
        with open(self.results_path + 'simulations.csv', 'w+') as file:
            pd.DataFrame(columns=columns).to_csv(file, index=False)
    if not os.path.exists(self.results_path + 'metoh_ex.csv'):
        with open(self.results_path + 'metoh_ex.csv', 'w+') as file:
            pd.DataFrame(columns=['Sample', 'Exchanged Metanol', 'N Blocked']).to_csv(file, index=False)
    if not os.path.exists(self.results_path + 'yields.csv'):
        with open(self.results_path + 'yields.csv', 'w+') as file:
            pd.DataFrame(columns=['Sample', 'Metabolite', 'Yield', 'N Blocked']).to_csv(file, index=False)

    n_ex = len(self.to_close)
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

    for n in trange(n_sim):
        max_del = random.randint(0, n_ex)
        index = random.sample(range(0, n_ex), max_del)

        with micom_model as model:
            for met in pd.Series(self.to_close).iloc[index]:
                rxn = model.reactions.get_by_id(met)
                rxn.lower_bound = -100
                rxn.upper_bound = 100

            micom_model, linear_solution, quadratic_solution, micom_solution, taxa_growth_rates = micom(
                micom_model=micom_model,
                micom_objective=micom_objective,
                trade_off=trade_off,
                taxa_db=sample_db,
                model=[]
                )

        external_fluxes = pd.DataFrame(index = external_metabolites, columns = sample_db.index)
        for rxn in micom_model.reactions:
            if '_EX' in rxn.id:
                metabolite = rxn.id[4:-5]
                genus_k = rxn.id[-4:]
                abundance_k = sample_db.loc[genus_k, 'relative_abundance']
                if metabolite in external_fluxes.index:
                    external_fluxes.loc[metabolite, genus_k] = 1000 * micom_solution[rxn.id] * abundance_k * list(rxn.metabolites)[0].cmol_by_mol

        external_fluxes.drop(index='CH4', inplace=True)
        external_fluxes.fillna(0, inplace = True)

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
        ic_long = self.ic_matrix_melt(interaction_matrix)
        ic_long.loc[:, 'Sample'] = sample
        ic_long.loc[:, 'N Blocked'] = max_del
        with open(self.results_path + 'simulations.csv', 'a+') as simulations:
            ic_long.to_csv(simulations, index=False, header=False)

        summary = micom_model.metabolites._METOH_pe_Meis.summary(solution=micom_solution)
        try:
            metoh_ex = summary.consuming_flux.loc['_EX_METOH_Meis', 'percent']
        except KeyError as e:
            print(e)
            metoh_ex = 0

        with open(self.results_path + 'metoh_ex.csv', 'a+') as file:
            pd.DataFrame(
                [[sample, metoh_ex, max_del]],
                columns=['Sample', 'Exchanged Metanol', 'N Blocked']
                ).to_csv(file, index=False, header=False)

        biomass_fluxes = micom_solution.fluxes.index.str.contains('biomass')
        biomass = 0
        for tax in micom_solution.fluxes.loc[biomass_fluxes].index:
            biomass += micom_solution[tax] * self.models[sample].metabolites.get_by_id(tax.replace('DM_', '')).cmol_by_mol
        biomass = round(biomass / carbon_uptake, 2)

        with open(self.results_path + 'yields.csv', 'a+') as file:
            pd.DataFrame(
                [[sample, 'Biomass', biomass, max_del]],
                columns=['Sample', 'Metabolite', 'Yield', 'N Blocked']
                ).to_csv(file, index=False, header=False)

        for rxn_id in yield_mets.keys():
            met_yield = round(abs(micom_solution[rxn_id] / carbon_uptake), 2)
            with open(self.results_path + 'yields.csv', 'a+') as file:
                pd.DataFrame(
                    [[sample, yield_mets[rxn_id], met_yield, max_del]],
                    columns=['Sample', 'Metabolite', 'Yield', 'N Blocked']
                    ).to_csv(file, index=False, header=False)

# %%codecell
n_sim = 100
n_ex = len(to_close)  # number of external metabolites to fit
rxns_to_fit = pd.Series(to_close)  # reactions to fit
# Initialize binary matrix
theta_sampler = pd.DataFrame(
    zeros((n_ex, n_sim)), index=to_close, columns=range(n_sim)
    )
# Set initial vector
for i in range(n_ex):
    # initial theta set to be uniformly distributed
    theta_sampler.iloc[i, 0] = float(binomial(size=1, n=1, p=uniform(0, 1))[0])

for k in range(n_sim):
    theta_old = theta_sampler.loc[:, k]
    for i in rxns_to_fit:
        print(theta_sampler.loc[i,k])
        # rxn = micom_model.reactions.get_by_id(i)
        # if theta_sampler.loc[i, k] == 1:
        #     rxn.lower_bound = -100
        #     rxn.upper_bound = 100
        # else:
        #     rxn.lower_bound = -100
        #     rxn.upper_bound = 100

index
