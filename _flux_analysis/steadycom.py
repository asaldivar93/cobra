from __future__ import absolute_import

from cobra.flux_analysis.helpers import normalize_cutoff
from optlang.symbolics import Zero
from copy import deepcopy
from numpy import arange


def _build_aggregate_constraints(model, taxa_db):
    for metabolite in model.metabolites:
        if '_ex' in metabolite.id:
            constraint = model.constraints.get(metabolite.id)
            variables = constraint.variables
            coefficients = constraint.get_linear_coefficients(variables)
            for var in coefficients.keys():
                coeff = coefficients[var]
                if 'ou_' in var.name or 'DM_' in var.name or '_ou' in var.name:
                    pass
                else:
                    try:
                        if '_reverse' in var.name:
                            pos = var.name.find('_reverse')
                            taxa_k = var.name[(pos - 4):pos]
                        else:
                            taxa_k = var.name[-4:]
                        abundance_k = taxa_db.loc[taxa_k, 'relative_abundance']
                        coefficients[var] = coeff * abundance_k
                    except:
                        print(var.name)
                constraint.set_linear_coefficients(coefficients)

            for rxn in metabolite.reactions:
                if 'ou_' in rxn.id or 'DM_' in rxn.id or '_ou' in var.name:
                    pass
                else:
                    taxa_k = rxn.id[-4:]
                    abundance_k = taxa_db.loc[taxa_k, 'relative_abundance']
                    rxn.lower_bound = rxn.lower_bound * abundance_k
                    rxn.upper_bound = rxn.upper_bound * abundance_k

    # for rxn in model.reactions:
    #     if 'biomass' in rxn.id:
    #         taxa_k = rxn.id[-4:]
    #         abundance_k = taxa_db.loc[taxa_k, 'relative_abundance']
    #         rxn.lower_bound = rxn.lower_bound * abundance_k
    #         rxn.upper_bound = rxn.upper_bound * abundance_k

    model.solver.update()


def add_steadycom_constraints(model, taxa_db):
    problem = model.problem

    growth_rate = problem.Variable(
        'growth_rate',
        lb=0,
    )
    model.add_cons_vars(growth_rate)
    model.solver.update()
    print('Created growth rate variable')

    print('Adding abundance variable and growth rate constraint for k taxa')
    growth_cons = []
    for taxa_k in taxa_db.index:
        growth_rxn_k = model.reactions.get_by_id(
            '_BIOMASS_{}'.format(taxa_k)
        )
        abundance_k = taxa_db.loc[taxa_k, 'relative_abundance']
        growth_k = problem.Constraint(
            1.0 * growth_rxn_k.forward_variable - 1.0 * growth_rate * abundance_k,
            name='growth_{}'.format(taxa_k),
            lb=0,
            ub=0
        )
        growth_cons.append(growth_k)

    model.add_cons_vars(growth_cons)
    model.solver.update()

    _build_aggregate_constraints(model, taxa_db)


def SteadyCom(model, max_biomass=1.0, growth_rate_0=0.5, max_iter=1000):
    problem = model.problem
    zero_cutoff = normalize_cutoff(model, None)
    objective = Zero

    for variable in model.variables:
        if 'abundance' in variable.name:
            objective += variable

    model.objective = problem.Objective(objective, direction='max')
    growth_rate = model.variables.get('growth_rate')
    growth_rate.lb = growth_rate_0
    growth_rate.ub = growth_rate_0

    growth_rate_ub = []
    growth_rate_lb = []
    aggregate_biomass = model.slim_optimize()
    print('finding growth_rate bounds')
    iter = 0
    while not all([growth_rate_ub, growth_rate_lb]) and iter <= max_iter:
        ag_b_last = aggregate_biomass
        if ag_b_last > max_biomass:
            growth_rate_ub = growth_rate_0
            growth_rate_0 = (1 + 0.01) * growth_rate_0
        else:
            growth_rate_lb = growth_rate_0
            growth_rate_0 = (1 - 0.01) * growth_rate_0

        growth_rate.lb = growth_rate_0
        growth_rate.ub = growth_rate_0
        aggregate_biomass = model.slim_optimize()
        print(growth_rate_ub, growth_rate_lb, aggregate_biomass)
        iter += 1

    print(growth_rate_ub, growth_rate_lb)
    # error = 1000
    # iter = 0
    # while error > zero_cutoff or iter < max_iter:
    #     growth_rate_m = (growth_rate_ub + growth_rate_lb) / 2
    #
    #     growth_rate.lb = growth_rate_ub
    #     growth_rate.ub = growth_rate_ub
    #     f_ub = model.slim_optimize() - 1
    #
    #     growth_rate.lb = growth_rate_lb
    #     growth_rate.ub = growth_rate_lb
    #     f_lb = model.slim_optimize() - 1
    #
    #     growth_rate.lb = growth_rate_m
    #     growth_rate.ub = growth_rate_m
    #     f_m = model.slim_optimize() - 1
    #
    #     if f_ub * f_m < 0:
    #         growth_rate_lb = growth_rate_m
    #     elif f_lb * f_m < 0:
    #         growth_rate_ub = growth_rate_m
    #
    #     growth_rate_m_new = (growth_rate_ub + growth_rate_lb) / 2
    #     error = abs(growth_rate_m_new - growth_rate_m) / growth_rate_m_new
    #     iter += 1
    #
    # growth_rate.lb = growth_rate_m_new
    # growth_rate.ub = growth_rate_m_new
    #
    # solution = model.optimize()
    #
    # return solution, growth_rate_m_new
