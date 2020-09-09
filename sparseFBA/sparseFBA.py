from cobra.flux_analysis.helpers import normalize_cutoff
from optlang.symbolics import Zero

model = c_model.copy()
# %% codecell
leaks = find_leaks(model)
leaks
leak_modes = find_leak_mode(model, leaks)

# %% codecell


def find_leaks(model):
    w_model = model.copy()
    prob = w_model.problem

    for rxn in w_model.reactions:
        if rxn.boundary:
            w_model.remove_reactions(rxn)

    for const in w_model.constraints:
        w_model.remove_cons_vars(const)

    obj_vars = []
    for met in w_model.metabolites:
        met_var = prob.Variable("aux_{}".format(met.id),
                                lb = 0
                                )

        const = prob.Constraint(Zero,
                                name = met.id,
                                lb = 0,
                                ub = 0
                                )
        w_model.add_cons_vars([met_var, const])

        rxn_coeffs = []
        for rxn in met.reactions:
            coeff = rxn.metabolites[met]
            rxn_coeffs.append([rxn.forward_variable, coeff])
            rxn_coeffs.append([rxn.reverse_variable, -coeff])

        rxn_coeffs.append([met_var, -1])
        rxn_coeffs = dict(rxn_coeffs)
        w_model.constraints.get(met.id).set_linear_coefficients(rxn_coeffs)

        obj_vars.extend([met_var])

    w_model.objective = prob.Objective(Zero, direction = 'max')
    w_model.objective.set_linear_coefficients({o: 1 for o in obj_vars})
    w_model.optimize()

    leaks = []
    for var in w_model.variables:
        if 'aux' in var.name:
            if var.primal > 0:
                leaks.extend([var.name.replace('aux_', '')])

    return leaks


def find_leak_mode(model, leaks = [], zero_cutoff = None):
    """ Solve the following LP problem

    minimize sum_r z_r**2
    s.t. S*v - y = 0
         v - z = 0
         z unbounded
         y >= 1 for met in leaks
         y >= 0 for met not in leaks
    """
    zero_cutoff = normalize_cutoff(model, zero_cutoff)
    w_model = model.copy()
    w_model.solver = 'cplex'
    prob = w_model.problem
    obj_vars = []
    objective = Zero
    rxn_vars_and_cons = []

    for rxn in w_model.reactions:
        if rxn.boundary:
            w_model.remove_reactions(rxn)

    for const in w_model.constraints:
        w_model.remove_cons_vars(const)

    for met in w_model.metabolites:
        met_var = prob.Variable("aux_{}".format(met.id),
                                lb = 0
                                )

        const = prob.Constraint(Zero,
                                name = met.id,
                                lb = 0,
                                ub = 0
                                )
        w_model.add_cons_vars([met_var, const])

        rxn_coeffs = []
        for rxn in met.reactions:
            coeff = rxn.metabolites[met]
            rxn_coeffs.append([rxn.forward_variable, coeff])
            rxn_coeffs.append([rxn.reverse_variable, -coeff])

        rxn_coeffs.append([met_var, -1])
        rxn_coeffs = dict(rxn_coeffs)
        w_model.constraints.get(met.id).set_linear_coefficients(rxn_coeffs)

    for rxn in w_model.reactions:
        rxn_var = prob.Variable("rxn_{}".format(rxn.id),
                                lb = 0,
                                ub = 0
                                )

        const = prob.Constraint(rxn.forward_variable + rxn.reverse_variable - rxn_var,
                                name = "rxn_{}".format(rxn.id),
                                lb = zero_cutoff
                                )
        rxn_vars_and_cons.extend([rxn_var, const])
        obj_vars.extend([rxn_var])
        objective += rxn_var**2

    w_model.add_cons_vars(rxn_vars_and_cons)
    # since solvers only accept Linear or Quadratic objectives
    w_model.objective = prob.Objective(objective, direction = 'min')

    # w_model.objective.set_linear_coefficients({o: 1 for o in obj_vars})
    print('finding leak modes')
    leak_modes = {}
    for leak in leaks:
        rxns_in_mode = []

        met_var = w_model.variables.get("aux_{}".format(leak))
        met_var.lb = 1
        sol = w_model.optimize()
        print(leak, sol.status)
        met_var.lb = 0

        rxns_in_mode = [[rxn.id, sol.fluxes[rxn.id]] for rxn in w_model.reactions if abs(sol.fluxes[rxn.id]) >= 2 * zero_cutoff]
        leak_modes[leak] = rxns_in_mode

    return leak_modes
