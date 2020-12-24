import pandas as pd


def search_biomass_components(model, path_to_biomass = 'data_files/biomass.csv'):
    biomass_in_model = pd.DataFrame(
        columns = ['class', 'name', 'mmol/g', 'met_id']
        )
    biomas_db = pd.read_csv(path_to_biomass)
    biomas_db
    biomas_db.set_index(
        'metacyc_identifier', inplace = True
    )

    for biomass_component in biomas_db.index:
        for metabolite in model.metabolites:
            bc = biomass_component.replace('|', '').replace('.', '').replace('-', '_').replace('+', '_')
            if bc + '_cy' in metabolite.id or bc + '_pe' in metabolite.id or bc + '_it' in metabolite.id:
                biomass_in_model = biomass_in_model.append(
                    pd.DataFrame(
                        [[biomas_db.loc[biomass_component, 'class'],
                          biomas_db.loc[biomass_component, 'name'],
                          biomas_db.loc[biomass_component, 'mmol/g'],
                          metabolite.id]],
                        columns = ['class', 'name', 'mmol/g', 'met_id']
                    )
                )
    biomass_in_model.reset_index(
        inplace = True
    )
    biomass_in_model.drop(
        columns=['index'], inplace=True
    )
    dup_filter = biomass_in_model['name'].duplicated()
    biomass_in_model.drop(
        biomass_in_model[dup_filter].index, inplace = True
    )

    return biomass_in_model


def rxns_of_metabolite(model, metabolite = '', in_solution = False):
    for rxn in model.reactions:
        met_in_rxn = [metabolite in met.id for met in rxn.metabolites]
        if any(met_in_rxn):
            print(rxn.subsystem, rxn.id, rxn.build_reaction_string())


def add_boundaries_for_simulation(
    model,
    carbon_uptake = 0.00625,
    media = [],
    demands = [],
    sinks = []
):
    if media:
        for substrate in media:
            met = model.metabolites.get_by_id(substrate)
            model.add_boundary(
                met, type = 'exchange'
            )
        medium = {rxn.id: 1000 for rxn in model.exchanges}
        if 'EX__CH4_ou' in medium.keys():
            medium['EX__CH4_ou'] = carbon_uptake

        model.medium = medium

    if demands:
        for met in demands:
            try:
                dm = model.metabolites.get_by_id(met)
            except:
                pass
            else:
                model.add_boundary(
                    dm, type = 'demand'
                )
    if sinks:
        for met in sinks:
            try:
                sk = model.metabolites.get_by_id(met)
            except:
                pass
            else:
                model.add_boundary(
                    sk, type = 'sink'
                )

    o2_ex = model.reactions.get_by_id('_ou_O2')
    o2_ex.upper_bound = 1000
    o2_ex.lower_bound = 0
    for rxn in model.reactions:
        if '_EX_O2' in rxn.id:
            # print(rxn.id)
            rxn.upper_bound = 1000
        elif '_DM_ATP_cy' in rxn.id:
            rxn.lower_bound = 3.5 / 1000
            rxn.upper_bound = 1000
        elif '_BIOMASS' in rxn.id:
            rxn.lower_bound = 0
            rxn.upper_bound = 2000

    for met in model.metabolites:
        if '_biomass' in met.id:
            model.add_boundary(
                met, type = 'demand', ub = 2000
                )


def unused_plotrelative_tradeoff():
    media = ['_CH4_ou', '_OXYGEN_MOLECULE_ou', '_NITRATE_ou', '_FE_2_ou',
             '_Pi_ou', '_SULFATE_ou', '_NA__ou', '_MG_2_ou', '_CO_2_ou',
             '_CL__ou']
    sample_db = taxa_database[sample].set_index('id')
    with open('data_files/fit_test.csv', 'w+') as file:
        pd.DataFrame(
            columns = ['pivot',
                       'target',
                       'relative_abundance',
                       'fitness',
                       'relative_fitness',
                       'relative_k',
                       'import_flux',
                       'export_flux']).to_csv(file, index = False)
    sample_db.index[-3:]

    for k in ['Meis']:
        with community_model as model:
            add_boundaries_for_simulation(model, media=media, carbon_uptake=0.00625)
            new = sample_db.copy()
            others = new['seq_counts'].loc[~(new['seq_counts'].index == k)].sum()
            print(others)
            if k == 'Chae':
                range = arange(0.7, 1, 0.05)
            else:
                range = arange(0.05, 1, 0.05)
            for r in range:
                new_counts = (r * others) / (1 - r)
                new['seq_counts'][k] = new_counts
                new['relative_abundance'] = new['seq_counts'] / new['seq_counts'].sum()
                relative_abundance = new['relative_abundance'][k]
                micom_model, linear_solution, quadratic_solution, micom_solution, taxa_growth_rates = micom(
                    model, 0.8, new
                )
                print(relative_abundance)
                if micom_solution.status == 'optimal':
                    average_fitness = 0
                    for gen in taxa_growth_rates:
                        gen_k = gen.id[-4:]
                        abundance_k = new.loc[gen_k, 'relative_abundance']
                        average_fitness += taxa_growth_rates[gen]
                        if gen_k == k:
                            fitness_k = taxa_growth_rates[gen]

                    for gen in taxa_growth_rates:
                        gen_k = gen.id[-4:]
                        abundance_k = new.loc[gen_k, 'relative_abundance']
                        fitness = taxa_growth_rates[gen]
                        relative_fitness = fitness / average_fitness
                        relative_k = fitness / fitness_k
                        pivot_genus = sample_db.loc[k, 'Genus']
                        target_genus = sample_db.loc[gen_k, 'Genus']
                        import_flux = 0
                        export_flux = 0
                        for rxn in micom_model.reactions:
                            ind_in_id = [ind in rxn.id for ind in ['EX', gen_k]]
                            if all(ind_in_id):
                                if micom_solution[rxn.id] > 0:
                                    import_flux += micom_solution[rxn.id] * list(rxn.metabolites)[0].cmol_by_mol
                                elif micom_solution[rxn.id] < 0:
                                    export_flux += micom_solution[rxn.id] * list(rxn.metabolites)[0].cmol_by_mol
                        import_flux = import_flux * abundance_k * 1000
                        export_flux = export_flux * abundance_k * 1000

                        with open('data_files/fit_test.csv', 'a+') as fitness_database:
                            pd.DataFrame(
                                [[pivot_genus,
                                  target_genus,
                                  relative_abundance,
                                  fitness,
                                  relative_fitness,
                                  relative_k,
                                  import_flux,
                                  export_flux]],
                                columns = ['pivot',
                                           'target',
                                           'relative_abundance',
                                           'fitness',
                                           'relative_fitness',
                                           'relative_k',
                                           'import_flux',
                                           'export_flux']).to_csv(fitness_database, header=False, index=False)
                else:
                    for gen in taxa_growth_rates:
                        gen_k = gen.id[-4:]
                        fitness = 0
                        relative_fitness = 0
                        pivot_genus = sample_db.loc[k, 'Genus']
                        target_genus = sample_db.loc[gen_k, 'Genus']
                        import_flux = 0
                        export_flux = 0
                        relative_k = 0
                        with open('data_files/fit_test.csv', 'a+') as fitness_database:
                            pd.DataFrame(
                                [[pivot_genus,
                                  target_genus,
                                  relative_abundance,
                                  fitness,
                                  relative_fitness,
                                  relative_k,
                                  import_flux,
                                  export_flux]],
                                columns = ['pivot',
                                           'target',
                                           'relative_abundance',
                                           'fitness',
                                           'relative_fitness',
                                           'relative_k',
                                           'import_flux',
                                           'export_flux']).to_csv(fitness_database, header=False, index=False)

    # %%codecell
    fit_test = pd.read_csv('data_files/fit_test.csv')
    fit_test.head(17)

    fit_test.tail(17)

    fit_test['export_flux'] = -fit_test['export_flux']
    fitness_database = pd.read_csv('data_files/fitness_database.csv')
    pivot = fitness_database.query("pivot=='Methylocystis'")
    px.line(pivot, x = 'relative_abundance', y = 'relative_fitness', line_group = 'target', color = 'target', log_y = True)
    px.line(fit_test, x = 'relative_abundance', y = 'import_flux', line_group = 'target', color = 'target', log_y = True)
    px.line(fit_test, x = 'relative_abundance', y = 'fitness', line_group = 'target', color = 'target', log_y = True)
