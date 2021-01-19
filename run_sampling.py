from cobra.flux_analysis import community

home_path = '/home/alexis/UAM/cobra/'
taxonomy_path = 'data_files/silva_taxonomy.csv'
counts_path = 'data_files/otu-table-corrected.tsv'
results_path = 'results/'

com = community.community(
    home_path=home_path, taxonomy_path=taxonomy_path, counts_path=counts_path, results_path=results_path
    )
media = ['_CH4_ou', '_OXYGEN_MOLECULE_ou', '_NITRATE_ou', '_FE_2_ou',
         '_Pi_ou', '_SULFATE_ou', '_NA__ou', '_MG_2_ou', '_CO_2_ou',
         '_CL__ou']
to_close = ['_EX_ETOH_Meis', '_EX_BUTANEDIOL_Meis', '_EX_ACET_Meis', '_EX_FORMATE_Meis',
            '_EX_PUTRESCINE_Meis', '_EX_PROPANE_1_2_DIOL_Meis', '_EX_GLT_Meis',
            '_EX_SUC_Meis', '_EX_FUM_Meis', '_EX_FRUCTOSE_6P_Meis', '_EX_CIT_Meis',
            '_EX_RIBOSE_5P_Meis', '_EX_MAL_Meis', '_EX_SUCROSE_Meis', '_EX_GLYCERATE_Meis',
            '_EX_FRUCTOSE_16_DIPHOSPHATE_Meis', '_EX_PYRUVATE_Meis', '_EX_GAP_Meis',
            '_EX_2_KETO_3_DEOXY_6_P_GLUCONATE_Meis', '_EX_CPD_2961_Meis',
            '_EX_PHOSPHO_ENOL_PYRUVATE_Meis', '_EX_G3P_Meis', '_EX_L_ALPHA_ALANINE_Meis',
            '_EX_SER_Meis', '_EX_THR_Meis', '_EX_GLN_Meis', '_EX_L_ASPARTATE_Meis']
yield_mets = dict(
    EX__OXYGEN_MOLECULE_ou = dict(
        name='O2', yields=dict(CIR_19=1.31, My_20=1.46), std=dict(CIR_19=0.1, My_20=0.1)
        ),
    _DM_CARBON_DIOXIDE_ex = dict(
        name='CO2', yields=dict(CIR_19=0.64, My_20=0.69), std=dict(CIR_19=0.1, My_20=0.1)
        )
    )

samples = dict(CIR_19 = 0.0063, My_20 = 0.00913)
com.set_samples(samples)
com.load_models(samples, to_close=to_close, close_ex=True)

com.fit_model(
    sample='CIR_19',
    media=media,
    carbon_uptake=0.00625,
    trade_off=0.9,
    yield_mets=yield_mets,
    n_sim=600,
    warmup=250,
    permute_init=True,
    dir = 'fit_16-01-2021_17-13-07'
)
