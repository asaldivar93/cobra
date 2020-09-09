
from cobra.flux_analysis.helpers import normalize_cutoff
from optlang.symbolics import Zero

for met in model.metabolites:
    if 'ACETALD' in met.id:
        print(met.id)

# %% codecell
leaks = find_leaks(model)

mode = find_leak_mode(model, leaks)

# %% codecell


def find_leaks(model):
    with model:
        prob = model.problem

        for rxn in model.reactions:
            if rxn.boundary:
                model.remove_reactions(rxn)

        for const in model.constraints:
            model.remove_cons_vars(const)

        obj_vars = []
        for met in model.metabolites:
            met_var = prob.Variable("aux_{}".format(met.id),
                                    lb = 0
                                    )

            const = prob.Constraint(Zero,
                                    name = met.id,
                                    lb = 0,
                                    ub = 0
                                    )
            model.add_cons_vars([met_var, const])

            rxn_coeffs = []
            for rxn in met.reactions:
                coeff = rxn.metabolites[met]
                rxn_coeffs.append([rxn.forward_variable, coeff])
                rxn_coeffs.append([rxn.reverse_variable, -coeff])

            rxn_coeffs.append([met_var, -1])
            rxn_coeffs = dict(rxn_coeffs)
            model.constraints.get(met.id).set_linear_coefficients(rxn_coeffs)

            obj_vars.extend([met_var])

        model.objective = prob.Objective(Zero, direction = 'max')
        model.objective.set_linear_coefficients({o: 1 for o in obj_vars})
        model.optimize()

        leaks = []
        for var in model.variables:
            if 'aux' in var.name:
                if var.primal > 0:
                    leaks.extend([var.name.replace('aux_', '')])

    return leaks


def find_leak_mode(model, leaks = [], zero_cutoff = None):
    """ Solve the following LP problem

    minimize sum_r z_r
    s.t. S*v - y = 0
         v**2 - z = 0
         z >= 0
         y >= 1 for met in leaks
         y >= 0 for met not in leaks
    """
    zero_cutoff = normalize_cutoff(model, zero_cutoff)

    with model:
        prob = model.problem
        obj_vars = []
        rxn_vars_and_cons = []

        for rxn in model.reactions:
            if rxn.boundary:
                model.remove_reactions(rxn)

        for const in model.constraints:
            model.remove_cons_vars(const)

        for met in model.metabolites:
            met_var = prob.Variable("aux_{}".format(met.id),
                                    lb = 0
                                    )

            const = prob.Constraint(Zero,
                                    name = met.id,
                                    lb = 0,
                                    ub = 0
                                    )
            model.add_cons_vars([met_var, const])

            rxn_coeffs = []
            for rxn in met.reactions:
                coeff = rxn.metabolites[met]
                rxn_coeffs.append([rxn.forward_variable, coeff])
                rxn_coeffs.append([rxn.reverse_variable, -coeff])

            rxn_coeffs.append([met_var, -1])
            rxn_coeffs = dict(rxn_coeffs)
            model.constraints.get(met.id).set_linear_coefficients(rxn_coeffs)

        for rxn in model.reactions:
            rxn_var = prob.Variable("rxn_{}".format(rxn.id),
                                    lb = 0,
                                    ub = 0
                                    )

            const = prob.Constraint((rxn.forward_variable + rxn.reverse_variable)**2 - rxn_var,
                                    name = "rxn_{}".format(rxn.id),
                                    lb = zero_cutoff
                                    )
            rxn_vars_and_cons.extend([rxn_var, const])
            obj_vars.extend([rxn_var])
        model.add_cons_vars(rxn_vars_and_cons)

        model.objective = prob.Objective(Zero, direction = 'max')
        model.objective.set_linear_coefficients({o: 1 for o in obj_vars})

        leak_modes = {}
        for leak in leaks:
            rxns_in_mode = []

            met_var = model.variables.get("aux_{}".format(leak))
            met_var.lb = 1
            model.optimize()
            met_var.lb = 0

            rxns_in_mode = [[rxn.id, rxn.flux] for rxn in model.reactions if abs(rxn.flux) >= zero_cutoff]
            leak_modes[leak] = rxns_in_mode

    return leak_modes
