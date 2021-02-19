
import pandas as pd
import arviz

from stan.run_stan import run_zinb, run_st
from cobra.flux_analysis import community
from tools.cobra_tools import gibbs_results
from numpy.random import binomial

import plotly.express as px
import plotly.graph_objects as go

# %% codecell
home_path = '/home/alexis/UAM/cobra/'
taxonomy_path = 'data_files/silva_taxonomy.csv'
counts_path = 'data_files/otu-table-corrected.tsv'
results_path = 'results/'

com = community.community(
    home_path=home_path,
    taxonomy_path=taxonomy_path,
    counts_path=counts_path,
    results_path=results_path,
    name='start'
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

# %% codecell
# Loas models for each sample
com.load_models(samples, to_close=to_close, close_ex=False, micom=True)

# %%codecell
############ Fit models using adaptive gibbs sampling ###############
com.fit_model(
    sample='My_20',
    media=media,
    carbon_uptake=0.00625,
    trade_off=0.9,
    yield_mets=yield_mets,
    n_sim=600,
    warmup=250,
    dir='fit_16-01-2021_16-05-00'
)

# %%codecell
############ Load sampling results for each sample ###############
fit_dirs = dict(
    CIR_19 = ['fit_16-01-2021_16-04-42', 'fit_16-01-2021_17-12-41', 'fit_16-01-2021_17-13-07', 'fit_16-01-2021_17-13-56'],
    My_20 = ['fit_16-01-2021_16-05-00', 'fit_16-01-2021_19-03-01', 'fit_16-01-2021_19-05-33', 'fit_16-01-2021_19-15-13']
)
# this functions combines all sampling runs in a single dataframe
gr = gibbs_results(
    results_path=results_path,
    fit_dirs=fit_dirs,
    fitted_mets=to_close
)
# Drop outliers
gr.interactions = gr.interactions.drop(gr.interactions[gr.interactions['Interaction Coefficient'] > 1].index)
gr.interactions = gr.interactions.drop(gr.interactions[gr.interactions['Interaction Coefficient'] < -1].index)
gr.yields = gr.yields.drop(gr.yields[gr.yields['Yield'] > 2].index)

# %%codecell
############ Plot Violin plot for Yields ###############
fig = go.Figure()
sample = 'CIR_19'
fig.add_trace(
    go.Violin(
        x=gr.yields.query("Sample==@sample")['Metabolite'],
        y=gr.yields.query("Sample==@sample")['Yield'],
        legendgroup=sample, scalegroup=sample, name=sample,
        side='negative'
    )
)
sample = 'My_20'
fig.add_trace(
    go.Violin(
        x=gr.yields.query("Sample==@sample")['Metabolite'],
        y=gr.yields.query("Sample==@sample")['Yield'],
        legendgroup=sample, scalegroup=sample, name=sample,
        side='positive'
    )
)
fig.update_traces(width=1, meanline=dict(visible=True))
fig.update_layout(template='none')
fig.show()

# %% codecell
fig_co2 = go.Figure()
for sample in ['CIR_19', 'My_20']:
    fig_co2.add_trace(
        go.Violin(
            x=gr.yields.query("Sample==@sample").query("Metabolite=='CO2'")['Yield'],
            showlegend=False,
            name=sample
        )
    )
fig_co2.update_traces(
    orientation='h', side='positive', width=2, points=False
    )
fig_co2.update_layout(
    yaxis_showgrid=True,
    xaxis_showgrid=False,
    xaxis_zeroline=False,
    template='none',
    margin=dict(l=60, r=20, t=20, b=50),
    font=dict(family='Calibri'),
    width=185,
    height=250
    )
fig_co2.update_xaxes(
    title=dict(
        text="CO<sub>2</sub>   (C-mol C-mol <sup>-1</sup>)",
        font=dict(size=12)
    )
)
fig_co2.add_shape(
    type='line',
    x0=0.64, y0=0,
    x1=0.64, y1=1,
    line=dict(color='Blue', dash='dash'),
    xref='x', yref='y'
)
fig_co2.add_shape(
    type='line',
    x0=0.69, y0=1,
    x1=0.69, y1=2,
    line=dict(color='Orange', dash='dash'),
    xref='x', yref='y'
)

fig_o2 = go.Figure()
for sample in ['CIR_19', 'My_20']:
    fig_o2.add_trace(
        go.Violin(
            x=gr.yields.query("Sample==@sample").query("Metabolite=='O2'")['Yield'],
            showlegend=False,
            name=sample
        )
    )
fig_o2.update_traces(
    orientation='h', side='positive', width=2, points=False
    )
fig_o2.update_yaxes(ticks='', showticklabels=False)
fig_o2.update_layout(
    yaxis_showgrid=True,
    xaxis_showgrid=False,
    xaxis_zeroline=False,
    template='none',
    margin=dict(l=30, r=20, t=20, b=50),
    font=dict(family='Calibri'),
    width=185,
    height=250
    )
fig_o2.update_xaxes(
    title=dict(
        text="O<sub>2</sub>   (mol C-mol <sup>-1</sup>)",
        font=dict(size=12)
    )
)
fig_o2.add_shape(
    type='line',
    x0=1.31, y0=0,
    x1=1.31, y1=1,
    line=dict(color='Blue', dash='dash'),
    xref='x', yref='y'
)
fig_o2.add_shape(
    type='line',
    x0=1.46, y0=1,
    x1=1.46, y1=2,
    line=dict(color='Orange', dash='dash'),
    xref='x', yref='y'
)

fig_x = go.Figure()
for sample in ['CIR_19', 'My_20']:
    fig_x.add_trace(
        go.Violin(
            x=gr.yields.query("Sample==@sample").query("Metabolite=='Biomass'")['Yield'] - 0.1,
            showlegend=False,
            name=sample
        )
    )
fig_x.update_traces(
    orientation='h', side='positive', width=2, points=False
    )
fig_x.update_yaxes(ticks='', showticklabels=False)
fig_x.update_layout(
    yaxis_showgrid=True,
    xaxis_showgrid=False,
    xaxis_zeroline=False,
    template='none',
    margin=dict(l=30, r=20, t=20, b=50),
    font=dict(family='Calibri'),
    width=185,
    height=250
    )
fig_x.update_xaxes(
    title=dict(
        text="Biomasa   (C-mol C-mol <sup>-1</sup>)",
        font=dict(size=12)
    )
)
fig_x.add_shape(
    type='line',
    x0=0.36, y0=0,
    x1=0.36, y1=1,
    line=dict(color='Blue', dash='dash'),
    xref='x', yref='y'
)
fig_x.add_shape(
    type='line',
    x0=0.31, y0=1,
    x1=0.31, y1=2,
    line=dict(color='Orange', dash='dash'),
    xref='x', yref='y'
)
fig_co2.write_image(results_path + 'plots/fitted/co2_yield.svg')
fig_o2.write_image(results_path + 'plots/fitted/o2_yield.svg')
fig_x.write_image(results_path + 'plots/fitted/biomass_yield.svg')
# %%codecell
############ Plot Violin plot for Metoh exchange ###############
fig_met = go.Figure()
for sample in ['CIR_19', 'My_20']:
    fig_met.add_trace(
        go.Violin(
            x=gr.ex_metoh.query("Sample==@sample")['Exchanged Metanol'],
            showlegend=False,
            name=sample
        )
    )

fig_met.update_traces(
    orientation='h',
    side='positive',
    width=2,
    points=False,
    meanline_visible=True
    )
fig_met.update_layout(
    yaxis_showgrid=True,
    xaxis_showgrid=False,
    xaxis_zeroline=False,
    template='none',
    margin=dict(l=60, r=20, t=20, b=50),
    font=dict(family='Calibri'),
    width=185,
    height=250
    )
fig_met.update_xaxes(
    title=dict(
        text="Metanol Secretado",
        font=dict(size=13)
    )
)
fig_met.write_image(results_path + 'plots/fitted/ex_metoh.svg')
# %%codecell
fig_n = px.histogram(
    gr.ex_metoh,
    x='N Open',
    color='Sample',
    template='none',
    marginal='box',
    labels={'N Open': 'Intercambios Activos', 'Sample': 'Muestra', 'count': 'Cuentas'}
)
fig_n.update_layout(
    font=dict(family='Calibri'),
    margin=dict(l=50, r=10, t=20, b=40),
    height=250,
    width=347
)
fig_n.update_yaxes(
    title=dict(text='Cuenta')
)
fig_n.show()
fig_n.write_image(results_path + 'plots/fitted/active_ex.svg')

# %%codecell
fig1 = px.box(
    gr.interactions.query("Reciver=='Meis'"),
    x='Giver', y='Interaction Coefficient',
    color='Sample',
    template='none',
    labels={'Giver': 'Donador', 'Interaction Coefficient': 'Coeficiente de Interacción'}
    )
fig2 = px.box(
    gr.interactions.query("Reciver=='Meus'"),
    x='Giver', y='Interaction Coefficient',
    color='Sample',
    template='none',
    labels={'Giver': 'Donador', 'Interaction Coefficient': 'Coeficiente de Interacción'}
    )
fig3 = px.box(
    gr.interactions.query("Reciver=='Hyum'"),
    x='Giver', y='Interaction Coefficient',
    color='Sample',
    template='none',
    labels={'Giver': 'Donador', 'Interaction Coefficient': 'Coeficiente de Interacción'}
    )


fig1.write_image(results_path + 'plots/fitted/ic_meis.svg')
fig2.write_image(results_path + 'plots/fitted/ic_meus.svg')
fig3.write_image(results_path + 'plots/fitted/ic_hyum.svg')

# %%codecell
tax_cir = ['Meus', 'Hyum', 'Seum']
tax_my = ['Meus', 'Hyum', 'Seum', 'Raer']
c_cir = ['rgb(102,194,164)'] * len(tax_cir)
colors_cir = dict(zip(tax_cir, c_cir))
c_my = ['rgb(102,194,164)'] * len(tax_my)
colors_my = dict(zip(tax_my, c_my))
fmeis1, fmeis2 = gr.plot_ic_confidence('Meis', tax_list = [tax_cir, tax_my], color = [colors_cir, colors_my])

fmeis1.update_layout(
    margin=dict(l=40, r=20, t=20, b=60),
    width=286,
    height=250,
)
fmeis1.update_xaxes(
    showticklabels=True, dtick=0.25,
    range=[-1, 1]
)
fmeis2.update_layout(
    margin=dict(l=40, r=20, t=20, b=60),
    width=286,
    height=250,
)
fmeis2.update_xaxes(
    showticklabels=True, dtick=0.25,
    range=[-1, 1]
)

fmeis1.write_image(results_path + '/plots/fitted/ic_confidence/meis-cir19.svg')
fmeis2.write_image(results_path + '/plots/fitted/ic_confidence/meis-my10.svg')

tax_cir = ['Meis', 'Hyum', 'Seum']
tax_my = ['Meis', 'Hyum', 'Seum', 'Raer']
c_cir = ['rgb(102,194,164)'] * len(tax_cir)
colors_cir = dict(zip(tax_cir, c_cir))
c_my = ['rgb(102,194,164)'] * len(tax_my)
colors_my = dict(zip(tax_my, c_my))
fmeis1, fmeis2 = gr.plot_ic_confidence('Meus', tax_list = [tax_cir, tax_my], color = [colors_cir, colors_my])
fmeis1.update_layout(
    margin=dict(l=40, r=20, t=20, b=60),
    width=286,
    height=250,
)
fmeis1.update_xaxes(
    showticklabels=True, dtick=0.25,
    range=[-1, 1]
)
fmeis2.update_layout(
    margin=dict(l=40, r=20, t=20, b=60),
    width=286,
    height=250,
)
fmeis2.update_xaxes(
    showticklabels=True, dtick=0.25,
    range=[-1, 1]
)
fmeis1.write_image(results_path + '/plots/fitted/ic_confidence/Meus-cir19.svg')
fmeis2.write_image(results_path + '/plots/fitted/ic_confidence/Meus-my10.svg')

tax_cir = ['Meis', 'Meus', 'Seum']
tax_my = ['Meis', 'Meus', 'Seum', 'Raer']
c_cir = ['rgb(102,194,164)'] * len(tax_cir)
colors_cir = dict(zip(tax_cir, c_cir))
c_my = ['rgb(102,194,164)'] * len(tax_my)
colors_my = dict(zip(tax_my, c_my))
fmeis1, fmeis2 = gr.plot_ic_confidence('Hyum', tax_list = [tax_cir, tax_my], color = [colors_cir, colors_my])

fmeis1.update_layout(
    margin=dict(l=40, r=20, t=20, b=60),
    width=286,
    height=250,
)
fmeis1.update_xaxes(
    showticklabels=True, dtick=0.25,
    range=[-1, 1]
)
fmeis2.update_layout(
    margin=dict(l=40, r=20, t=20, b=60),
    width=286,
    height=250,
)
fmeis2.update_xaxes(
    showticklabels=True, dtick=0.25,
    range=[-1, 1]
)
fmeis1.write_image(results_path + '/plots/fitted/ic_confidence/Hyum-cir19.svg')
fmeis2.write_image(results_path + '/plots/fitted/ic_confidence/Hyum-my10.svg')

tax_cir = ['Meis', 'Meus', 'Hyum']
tax_my = ['Meis', 'Meus', 'Hyum', 'Raer']
c_cir = ['rgb(102,194,164)'] * len(tax_cir)
colors_cir = dict(zip(tax_cir, c_cir))
c_my = ['rgb(102,194,164)'] * len(tax_my)
colors_my = dict(zip(tax_my, c_my))
fmeis1, fmeis2 = gr.plot_ic_confidence('Seum', tax_list = [tax_cir, tax_my], color = [colors_cir, colors_my])

fmeis1.update_layout(
    margin=dict(l=40, r=20, t=20, b=60),
    width=286,
    height=250,
)
fmeis1.update_xaxes(
    showticklabels=True, dtick=0.25,
    range=[-1, 1]
)
fmeis2.update_layout(
    margin=dict(l=40, r=20, t=20, b=60),
    width=286,
    height=250,
)
fmeis2.update_xaxes(
    showticklabels=True, dtick=0.25,
    range=[-1, 1]
)
fmeis1.write_image(results_path + '/plots/fitted/ic_confidence/Seum-cir19.svg')
fmeis2.write_image(results_path + '/plots/fitted/ic_confidence/Seum-my10.svg')

tax_cir = ['Meis', 'Meus', 'Hyum']
tax_my = ['Meis', 'Meus', 'Hyum', 'Seum']
c_cir = ['rgb(102,194,164)'] * len(tax_cir)
colors_cir = dict(zip(tax_cir, c_cir))
c_my = ['rgb(102,194,164)'] * len(tax_my)
colors_my = dict(zip(tax_my, c_my))
fmeis1, fmeis2 = gr.plot_ic_confidence('Raer', tax_list = [tax_cir, tax_my], color = [colors_cir, colors_my])

fmeis1.update_layout(
    margin=dict(l=40, r=20, t=20, b=60),
    width=286,
    height=250,
)
fmeis1.update_xaxes(
    showticklabels=True, dtick=0.25,
    range=[-1, 1]
)
fmeis2.update_layout(
    margin=dict(l=40, r=20, t=20, b=60),
    width=286,
    height=250,
)
fmeis2.update_xaxes(
    showticklabels=True, dtick=0.25,
    range=[-1, 1]
)
fmeis1.write_image(results_path + '/plots/fitted/ic_confidence/Raer-cir19.svg')
fmeis2.write_image(results_path + '/plots/fitted/ic_confidence/Raer-my10.svg')

# %%codecell
############ Fit active exchanges from gibbs sampling to a bernoulli distribution ###############
gr.fit_rxnsdist_to_bernoulli('CIR_19')
gr.fit_rxnsdist_to_bernoulli('My_20')

# %%codecell
############ Plot ridgline for theta credible intervals for exchanges ###############
sample = 'CIR_19'
gr.theta_dist[sample] = pd.read_csv(
    gr.results_path + 'exchanges/fitted/theta_{}.csv'.format(sample), index_col='Unnamed: 0'
)
met_names = ['ETOH', 'BUTANEDIOL', 'ACET', 'FORMATE', 'PUTRESCINE', 'PROPANEDIOL',
             'GLT', 'SUC', 'FUM', 'FRUCTOSE-6P', 'CIT', 'RIBOSE-5P', 'MAL', 'SUCROSE',
             'GLYCERATE', 'FRUCTOSE-16-DP', 'PYRUVATE', 'GAP', 'KETO-GLUCONATE',
             'DIHYDROLIPOATE', 'PEP', 'G3P', 'ALANINE', 'SER', 'THR', 'GLN', 'ASPARTATE', 'Null']
names = dict(zip(gr.theta_dist['CIR_19'].index, met_names))

true_diff = ['_EX_ETOH_Meis', '_EX_BUTANEDIOL_Meis', '_EX_ACET_Meis',
             '_EX_GLYCERATE_Meis', '_EX_G3P_Meis', '_EX_L_ALPHA_ALANINE_Meis',
             '_EX_SER_Meis', '_EX_L_ASPARTATE_Meis']
uncertain_diff = ['_EX_CIT_Meis', '_EX_FRUCTOSE_16_DIPHOSPHATE_Meis']
colors = dict()
for met in gr.theta_dist['CIR_19'].index:
    if met == 'Null':
        colors[met] = 'rgb(200, 10, 10)'
    elif met in true_diff:
        colors[met] = 'rgb(43,140,190)'
    elif met in uncertain_diff:
        colors[met] = 'rgb(102,194,164)'
    else:
        colors[met] = 'rgb(253,212,158)'
fig_conf_cir19 = gr.plot_fitted_rxns_confidence('CIR_19', colors, names)
fig_conf_cir19.update_layout(
    margin=dict(l=130, r=20, t=20, b=40),
    width=286,
    height=500,
    title=dict(
        text='CIR_19',
        font=dict(size=12,family='Calibri'),
        x=0.8
    ),
)
fig_conf_cir19.update_xaxes(
    title=dict(
        text='Theta',
        font=dict(family='Calibri')
    )
)
fig_conf_cir19.write_image(results_path + 'plots/fitted/rxn_conf_1.svg')

sample = 'My_20'
gr.theta_dist[sample] = pd.read_csv(
    gr.results_path + 'exchanges/fitted/theta_{}.csv'.format(sample), index_col='Unnamed: 0'
)
true_diff = ['_EX_ETOH_Meis', '_EX_BUTANEDIOL_Meis']
uncertain_diff = ['_EX_FORMATE_Meis', '_EX_FUM_Meis', '_EX_L_ALPHA_ALANINE_Meis', '_EX_SER_Meis']
colors = dict()
for met in gr.theta_dist['CIR_19'].index:
    if met == 'Null':
        colors[met] = 'rgb(200, 10, 10)'
    elif met in true_diff:
        colors[met] = 'rgb(43,140,190)'
    elif met in uncertain_diff:
        colors[met] = 'rgb(102,194,164)'
    else:
        colors[met] = 'rgb(253,212,158)'

fig_conf_my_20 = gr.plot_fitted_rxns_confidence('My_20', colors, names)
fig_conf_my_20.update_yaxes(
    ticks='',
    showticklabels=False
    )
fig_conf_my_20.update_layout(
    margin=dict(l=20, r=20, t=20, b=40),
    width=186,
    height=500,
    title=dict(
        text='My_20',
        font=dict(size=12, family='Calibri'),
        x=0.6
        )
    )
fig_conf_my_20.update_xaxes(
    title=dict(
        text='Theta',
        font=dict(family='Calibri')
    ),
    range=[0.3,0.6]
)
fig_conf_my_20.write_image(results_path + 'plots/fitted/rxn_conf_2.svg')


# %%codecell
sample = 'My_20'
sample_filter = com.taxa_database.loc[:, 'sample'] == sample
sample_db = com.taxa_database[sample_filter].set_index('id')
ic_matrix = pd.DataFrame(index=sample_db.index, columns=sample_db.index)
gr.get_icmatrix_from_fit(sample, ic_matrix)

# %%codecell
com.name = 'all_zeros'
fig_yields, y_df = com.plot_yields(
    media=media,
    carbon_uptake=0.00625,
    metabolites=yield_mets,
    trade_off=0.9
)
fig_yields.show()

com.plot_ex_fluxes(
    media=media,
    carbon_uptake=0.00625,
    trade_off=0.9
)
fig_ic = com.plot_interaction_matrix((9, 4))
fig_ic.show()

fig_metoh = com.plot_metoh_ex(
    media=media,
    carbon_uptake=0.00625
)
fig_metoh.show()

# %%codecell
com.name = 'bernoulli'
actives = dict()
for sample in ['CIR_19', 'My_20']:
    actives[sample] = pd.Series(index=to_close)
    for rxn in to_close:
        theta = gr.theta_samples[sample].loc[rxn, :].mean()
        if theta < 0.45:
            actives[sample][rxn] = 0
        else:
            actives[sample][rxn] = int(binomial(size=1, n=1, p=theta))
com.update_exchanges('CIR_19', actives['CIR_19'])
com.update_exchanges('My_20', actives['My_20'])

fig_yields_2, y_df = com.plot_yields(
    media=media,
    carbon_uptake=0.00625,
    metabolites=yield_mets,
    trade_off=0.9
)
fig_yields_2.show()

com.plot_ex_fluxes(
    media=media,
    carbon_uptake=0.00625,
    trade_off=0.9
)
fig_ic_2 = com.plot_interaction_matrix((9, 4))
fig_ic_2.show()

fig_metoh_2 = com.plot_metoh_ex(
    media=media,
    carbon_uptake=0.00625
)
fig_metoh_2.show()

# %%codecell
ic_cir = pd.read_csv(results_path + 'interaction_matrix/all_zeros/ic_CIR_19.csv', index_col=0)
ic_my = pd.read_csv(results_path + '/interaction_matrix/all_zeros/ic_My_20.csv', index_col=0)
fit_cir_pos = run_zinb(ic_cir, chains = 6, iter = 2000, positives = True)
fit_cir_neg = run_zinb(ic_cir, chains = 6, iter = 2000, positives = False)
fit_my_pos = run_zinb(ic_my, chains = 6, iter = 2000, positives = True)
fit_my_neg = run_zinb(ic_my, chains = 6, iter = 2000, positives = False)
print(fit_cir_pos)
print(fit_cir_neg)
print(fit_my_pos)
print(fit_my_neg)
arviz.plot_posterior(fit_cir_pos, var_names = ['mu', 'phi', 'T_max'])
arviz.plot_posterior(fit_my_pos, var_names = ['mu', 'phi', 'T_max'])


import seaborn as sns
from tools.plotting import plot_heatmap
a = pd.read_csv(results_path + 'exchanges/all_zeros/import_clr_CIR_19.csv', index_col='Unnamed: 0')
b = pd.read_csv(results_path + 'exchanges/all_zeros/export_clr_CIR_19.csv', index_col='Unnamed: 0')

plot_heatmap(a, cmap = sns.diverging_palette(271, 130, s = 80, l = 47, sep = 20, as_cmap = True))
plot_heatmap(b, cmap = sns.diverging_palette(271, 130, s = 80, l = 47, sep = 20, as_cmap = True))
