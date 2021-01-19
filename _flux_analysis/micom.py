from __future__ import absolute_import

from cobra.flux_analysis.helpers import normalize_cutoff
from optlang.symbolics import Zero
from copy import deepcopy
from numpy import arange

import pandas as pd


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
                    rxn.lower_bound = rxn.lower_bound * abundance_k * 100
                    rxn.upper_bound = rxn.upper_bound * abundance_k * 100
    model.solver.update()


def _build_micom_constraints(model, taxa_db):
    micom_model = deepcopy(model)
    _build_aggregate_constraints(micom_model, taxa_db)

    return micom_model


def _micom_linear_optimization(micom_model, linear_objective):
    micom_model.objective = linear_objective
    linear_solution = micom_model.optimize(objective_sense = 'maximize')
    max_growth_rate = linear_solution.objective_value

    return max_growth_rate, linear_solution


def _build_micom_objective(model, taxa_db):
    linear_objective = {}
    quadratic_objective = Zero
    for taxa_k in taxa_db.index:
        for reaction in model.reactions:
            if 'BIOMASS' in reaction.id and taxa_k in reaction.id:
                abundance_k = taxa_db.loc[taxa_k, 'relative_abundance']
                linear_objective[reaction] = abundance_k * 100
                quadratic_objective += 1.0 * reaction.flux_expression ** 2

    micom_objective = [linear_objective, quadratic_objective]
    return micom_objective


def _micom_quadratic_optimization(model, max_growth_rate, trade_off, quadratic_objective, linear_objective):
    problem = model.problem
    community_growth_rate = Zero
    with model as quadratic_model:
        for reaction in linear_objective.keys():
            community_growth_rate += reaction.flux_expression * linear_objective[reaction]

        max_growth_tradeoff = problem.Constraint(
            community_growth_rate - trade_off * max_growth_rate,
            name = 'trade_off',
            lb = 0
        )
        quadratic_model.add_cons_vars(max_growth_tradeoff)
        quadratic_model.solver.update()

        objective = problem.Objective(
            quadratic_objective,
            direction = 'min'
        )
        quadratic_model.objective = objective
        quadratic_solution = quadratic_model.optimize()

    taxa_growth_rates = {}
    for reaction in linear_objective.keys():
        taxa_growth_rates[reaction] = quadratic_solution.fluxes[reaction.id]

    return quadratic_solution, taxa_growth_rates


def micom(model, trade_off, taxa_db, micom_model=[], micom_objective = []):
    if not micom_model:
        micom_model = _build_micom_constraints(
            model, taxa_db
        )
    if not micom_objective:
        micom_objective = _build_micom_objective(
            micom_model, taxa_db
        )

    linear_objective, quadratic_objective = micom_objective

    max_growth_rate, linear_solution = _micom_linear_optimization(
        micom_model, linear_objective
    )

    problem = micom_model.problem
    quadratic_solution, taxa_growth_rates = _micom_quadratic_optimization(
        micom_model, max_growth_rate, trade_off, quadratic_objective, linear_objective
    )

    taxa_growth_constraints = []
    for reaction in taxa_growth_rates.keys():
        growth_rate_k = problem.Constraint(
            1.0 * reaction.flux_expression - taxa_growth_rates[reaction],
            ub = 0
        )
        taxa_growth_constraints.append(growth_rate_k)

    micom_model.add_cons_vars(taxa_growth_constraints)
    micom_model.solver.update()
    max_growth_rate, micom_solution = _micom_linear_optimization(micom_model, linear_objective)
    if micom_solution.status == 'optimal':
        for reaction in taxa_growth_rates.keys():
            taxa_growth_rates[reaction] = micom_solution.fluxes[reaction.id]

    return micom_model, linear_solution, quadratic_solution, micom_solution, taxa_growth_rates


def micom_tradeoff(model, taxa_db):
    zero_cutoff = normalize_cutoff(model, None)
    nonzero_biomass = pd.DataFrame(
        columns = ['Trade off', 'Growth Rate']
    )
    for trade_off in arange(0.1, 1, 0.1):
        with model as m:
            micom_model, linear_solution, quadratic_solution, micom_solution, taxa_growth_rates = micom(
                m, trade_off, taxa_db
            )

        if micom_solution.status == 'optimal':
            for reaction in taxa_growth_rates.keys():
                taxa_growth_rates[reaction] = micom_solution.fluxes[reaction.id]

            non_zero = sum([taxa_growth_rates[k] > zero_cutoff for k in taxa_growth_rates.keys()])
        else:
            non_zero = 0

        for k in taxa_growth_rates.keys():
            nonzero_biomass = nonzero_biomass.append(
                pd.DataFrame(
                    [[trade_off, taxa_growth_rates[k]]],
                    columns = ['Trade off', 'Growth Rate']
                )
            )
        print(trade_off, micom_solution.status, non_zero)
    return nonzero_biomass
