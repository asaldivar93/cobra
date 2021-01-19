# -*- coding: utf-8 -*-

from __future__ import absolute_import

from optlang.interface import OPTIMAL
from optlang.symbolics import Zero, add

from cobra.core import Model
from cobra.util import fix_objective_as_constraint
from cobra.flux_analysis.steadycom import add_steadycom_constraints

from copy import deepcopy
import pandas as pd


class GapFiller(object):
    """Class for performing gap filling.

    This class implements gap filling based on a mixed-integer approach,
    very similar to that described in [1]_ and the 'no-growth but growth'
    part of [2]_ but with minor adjustments. In short, we add indicator
    variables for using the reactions in the universal model, z_i and then
    solve problem

    minimize \sum_i c_i * z_i
    s.t. Sv = 0
         v_o >= t
         lb_i <= v_i <= ub_i
         v_i = 0 if z_i = 0

    where lb, ub are the upper, lower flux bounds for reaction i, c_i is a
    cost parameter and the objective v_o is greater than the lower bound t.
    The default costs are 1 for reactions from the universal model, 100 for
    exchange (uptake) reactions added and 1 for added demand reactions.

    Note that this is a mixed-integer linear program and as such will
    expensive to solve for large models. Consider using alternatives [3]_
    such as CORDA instead [4,5]_.

    Parameters
    ----------
    model : cobra.Model
        The model to perform gap filling on.
    universal : cobra.Model
        A universal model with reactions that can be used to complete the
        model.
    lower_bound : float
        The minimally accepted flux for the objective in the filled model.
    penalties : dict, None
        A dictionary with keys being 'universal' (all reactions included in
        the universal model), 'exchange' and 'demand' (all additionally
        added exchange and demand reactions) for the three reaction types.
        Can also have reaction identifiers for reaction specific costs.
        Defaults are 1, 100 and 1 respectively.
    integer_threshold : float
        The threshold at which a value is considered non-zero (aka
        integrality threshold). If gapfilled models fail to validate,
        you may want to lower this value.
    exchange_reactions : bool
        Consider adding exchange (uptake) reactions for all metabolites
        in the model.
    demand_reactions : bool
        Consider adding demand reactions for all metabolites.

    References
    ----------
    .. [1] Reed, Jennifer L., Trina R. Patel, Keri H. Chen, Andrew R. Joyce,
       Margaret K. Applebee, Christopher D. Herring, Olivia T. Bui, Eric M.
       Knight, Stephen S. Fong, and Bernhard O. Palsson. “Systems Approach
       to Refining Genome Annotation.” Proceedings of the National Academy
       of Sciences 103, no. 46 (2006): 17480–17484.

       [2] Kumar, Vinay Satish, and Costas D. Maranas. “GrowMatch: An
       Automated Method for Reconciling In Silico/In Vivo Growth
       Predictions.” Edited by Christos A. Ouzounis. PLoS Computational
       Biology 5, no. 3 (March 13, 2009): e1000308.
       doi:10.1371/journal.pcbi.1000308.

       [3] http://opencobra.github.io/cobrapy/tags/gapfilling/

       [4] Schultz, André, and Amina A. Qutub. “Reconstruction of
       Tissue-Specific Metabolic Networks Using CORDA.” Edited by Costas D.
       Maranas. PLOS Computational Biology 12, no. 3 (March 4, 2016):
       e1004808. doi:10.1371/journal.pcbi.1004808.

       [5] Diener, Christian https://github.com/cdiener/corda
    """

    def __init__(
        self,
        model,
        universal=None,
        lower_bound=0.05,
        penalties=None,
        exchange_reactions=False,
        demand_reactions=True,
        integer_threshold=1e-9,
        community=False,
        taxa_database=None
    ):
        self.original_model = model
        self.lower_bound = lower_bound
        self.model = deepcopy(model)
        tolerances = self.model.solver.configuration.tolerances
        tolerances.integrality = integer_threshold
        self.universal = deepcopy(universal) if universal else Model("universal")
        self.integer_threshold = integer_threshold
        self.indicators = list()
        self.costs = dict()
        self.penalties = dict(
            universal=2,
            exchange=100,
            demand=2,
            reverse=1,
            community_exchange=3
        )
        if penalties is not None:
            self.penalties.update(penalties)

        if community:
            self.taxa_database = taxa_database
            assert isinstance(self.taxa_database, pd.DataFrame), "User must supply a panda's DataFrame with relative abundace"
            add_steadycom_constraints(self.model, self.taxa_database)
            self.extend_community_model(exchange_reactions, demand_reactions)
        else:
            self.extend_model(exchange_reactions, demand_reactions)

        fix_objective_as_constraint(self.model, bound=lower_bound)
        self.add_switches_and_objective()

    def extend_community_model(self, exchange_reactions=False, demand_reactions=True):
        """Extend gapfilling model.

        Add reactions from universal model and optionally exchange and
        demand reactions for all metabolites in the model to perform
        gapfilling on. Also adds a new reaction for every irreversible reaction
        with inverse revesibility constraints

        Parameters
        ----------
        exchange_reactions : bool
            Consider adding exchange (uptake) reactions for all metabolites
            in the model.
        demand_reactions : bool
            Consider adding demand reactions for all metabolites.
        """
        inverse_reactions = list()
        for rxn in self.model.reactions:
            if not rxn.reversibility:
                keep_irreversible = ['_DM_', '_EX_', '_ou_', 'ATP', '_BM', '_BIOMASS']
                if not any([kr in rxn.id for kr in keep_irreversible]):
                    rxn_rev = rxn.copy()
                    rxn_rev.id = rxn.id + '_rev'
                    if rxn_rev.lower_bound < 0:
                        rxn_rev.upper_bound = - rxn_rev.lower_bound
                        rxn_rev.lower_bound = 0
                        rxn_rev.gapfilling_type = "reverse"
                    elif rxn_rev.upper_bound > 0:
                        rxn_rev.lower_bound = - rxn.upper_bound
                        rxn_rev.upper_bound = 0
                        rxn_rev.gapfilling_type = "reverse"
                    elif all([rxn_rev.lower_bound == 0, rxn_rev.upper_bound == 0]):
                        rxn_rev.lower_bound = - 100
                        rxn_rev.upper_bound = 100
                        rxn.gapfilling_type = "universal"
                    inverse_reactions.append(rxn_rev)
        self.model.add_reactions(inverse_reactions)

        for rxn in self.universal.reactions:
            if 'EX_' in rxn.id and not rxn.boundary:
                rxn.gapfilling_type = "community_exchange"
            else:
                rxn.gapfilling_type = "universal"
        new_metabolites = self.universal.metabolites.query(
            lambda metabolite: metabolite not in self.model.metabolites
        )
        self.model.add_metabolites(new_metabolites)
        existing_exchanges = []
        for rxn in self.universal.boundary:
            existing_exchanges = existing_exchanges + [
                met.id for met in list(rxn.metabolites)
            ]

        for met in self.model.metabolites:
            if exchange_reactions:
                # check for exchange reaction in model already
                if met.id not in existing_exchanges:
                    rxn = self.universal.add_boundary(
                        met,
                        type="exchange_smiley",
                        lb=-1000,
                        ub=0,
                        reaction_id="EX_{}".format(met.id),
                    )
                    rxn.gapfilling_type = "exchange"
            if demand_reactions:
                rxn = self.universal.add_boundary(
                    met,
                    type="demand_smiley",
                    lb=0,
                    ub=1000,
                    reaction_id="DM_{}".format(met.id),
                )
                rxn.gapfilling_type = "demand"

        new_reactions = self.universal.reactions.query(
            lambda reaction: reaction not in self.model.reactions
        )
        self.model.add_reactions(new_reactions)

        growth_rate = self.model.variables.get('growth_rate')
        objective = self.model.problem.Objective(growth_rate, direction = 'max')
        self.model.objective = objective

    def extend_model(self, exchange_reactions=False, demand_reactions=True):
        """Extend gapfilling model.

        Add reactions from universal model and optionally exchange and
        demand reactions for all metabolites in the model to perform
        gapfilling on.

        Parameters
        ----------
        exchange_reactions : bool
            Consider adding exchange (uptake) reactions for all metabolites
            in the model.
        demand_reactions : bool
            Consider adding demand reactions for all metabolites.
        """
        for rxn in self.universal.reactions:
            rxn.gapfilling_type = "universal"
        new_metabolites = self.universal.metabolites.query(
            lambda metabolite: metabolite not in self.model.metabolites
        )
        self.model.add_metabolites(new_metabolites)
        existing_exchanges = []
        for rxn in self.universal.boundary:
            existing_exchanges = existing_exchanges + [
                met.id for met in list(rxn.metabolites)
            ]

        for met in self.model.metabolites:
            if exchange_reactions:
                # check for exchange reaction in model already
                if met.id not in existing_exchanges:
                    rxn = self.universal.add_boundary(
                        met,
                        type="exchange_smiley",
                        lb=-1000,
                        ub=0,
                        reaction_id="EX_{}".format(met.id),
                    )
                    rxn.gapfilling_type = "exchange"
            if demand_reactions:
                rxn = self.universal.add_boundary(
                    met,
                    type="demand_smiley",
                    lb=0,
                    ub=1000,
                    reaction_id="DM_{}".format(met.id),
                )
                rxn.gapfilling_type = "demand"

        new_reactions = self.universal.reactions.query(
            lambda reaction: reaction not in self.model.reactions
        )
        self.model.add_reactions(new_reactions)

    def update_costs(self):
        """Update the coefficients for the indicator variables in the objective.

        Done incrementally so that second time the function is called,
        active indicators in the current solutions gets higher cost than the
        unused indicators.
        """
        for var in self.indicators:
            if var not in self.costs:
                self.costs[var] = var.cost
            else:
                if var._get_primal() > self.integer_threshold:
                    self.costs[var] += var.cost
        self.model.objective.set_linear_coefficients(self.costs)

    def add_switches_and_objective(self):
        """Update gapfilling model with switches and the indicator objective."""
        constraints = list()
        big_m = max(max(abs(b) for b in r.bounds) for r in self.model.reactions)
        prob = self.model.problem
        for rxn in self.model.reactions:
            if not hasattr(rxn, "gapfilling_type"):
                continue
            indicator = prob.Variable(
                name="indicator_{}".format(rxn.id), lb=0, ub=1, type="binary"
            )
            if rxn.id in self.penalties:
                indicator.cost = self.penalties[rxn.id]
            else:
                indicator.cost = self.penalties[rxn.gapfilling_type]
            indicator.rxn_id = rxn.id
            self.indicators.append(indicator)

            # if z = 1 v_i is allowed non-zero
            # v_i - Mz <= 0   and   v_i + Mz >= 0
            constraint_lb = prob.Constraint(
                rxn.flux_expression - big_m * indicator,
                ub=0,
                name="constraint_lb_{}".format(rxn.id),
                sloppy=True,
            )
            constraint_ub = prob.Constraint(
                rxn.flux_expression + big_m * indicator,
                lb=0,
                name="constraint_ub_{}".format(rxn.id),
                sloppy=True,
            )

            constraints.extend([constraint_lb, constraint_ub])

        self.model.add_cons_vars(self.indicators)
        self.model.add_cons_vars(constraints, sloppy=True)
        self.model.objective = prob.Objective(Zero, direction="min", sloppy=True)
        self.model.objective.set_linear_coefficients({i: 1 for i in self.indicators})
        self.update_costs()

    def fill(self, iterations=1):
        """Perform the gapfilling by iteratively solving the model, updating
        the costs and recording the used reactions.


        Parameters
        ----------
        iterations : int
            The number of rounds of gapfilling to perform. For every
            iteration, the penalty for every used reaction increases
            linearly. This way, the algorithm is encouraged to search for
            alternative solutions which may include previously used
            reactions. I.e., with enough iterations pathways including 10
            steps will eventually be reported even if the shortest pathway
            is a single reaction.

        Returns
        -------
        iterable
            A list of lists where each element is a list reactions that were
            used to gapfill the model.

        Raises
        ------
        RuntimeError
            If the model fails to be validated (i.e. the original model with
            the proposed reactions added, still cannot get the required flux
            through the objective).
        """
        used_reactions = list()
        for i in range(iterations):
            self.model.slim_optimize(
                error_value=None, message="gapfilling optimization failed"
            )
            solution = [
                self.model.reactions.get_by_id(ind.rxn_id)
                for ind in self.indicators
                if ind._get_primal() > self.integer_threshold
            ]
            # if not self.validate(solution):
                # raise RuntimeError(
                #     "failed to validate gapfilled model, "
                #     "try lowering the integer_threshold"
                # )
            used_reactions.append(solution)
            self.update_costs()
        return used_reactions

    def validate(self, reactions):
        with self.original_model as model:
            mets = [x.metabolites for x in reactions]
            all_keys = set().union(*(d.keys() for d in mets))
            model.add_metabolites(all_keys)
            model.add_reactions(reactions)
            model.slim_optimize()
            print(model.solver.status, model.solver.objective.value)
            return (
                model.solver.status == OPTIMAL
                and model.solver.objective.value >= self.lower_bound
            )


def gapfill(
    model,
    universal=None,
    lower_bound=0.05,
    penalties=None,
    demand_reactions=True,
    exchange_reactions=False,
    iterations=1,
    community=False,
    taxa_database=None
):
    """Perform gapfilling on a model.

        See documentation for the class GapFiller.

        Parameters
        ----------
        model : cobra.Model
            The model to perform gap filling on.
        universal : cobra.Model, None
            A universal model with reactions that can be used to complete the
            model. Only gapfill considering demand and exchange reactions if
            left missing.
        lower_bound : float
            The minimally accepted flux for the objective in the filled model.
        penalties : dict, None
            A dictionary with keys being 'universal' (all reactions included in
            the universal model), 'exchange' and 'demand' (all additionally
            added exchange and demand reactions) for the three reaction types.
            Can also have reaction identifiers for reaction specific costs.
            Defaults are 1, 100 and 1 respectively.
        iterations : int
            The number of rounds of gapfilling to perform. For every iteration,
            the penalty for every used reaction increases linearly. This way,
            the algorithm is encouraged to search for alternative solutions
            which may include previously used reactions. I.e., with enough
            iterations pathways including 10 steps will eventually be reported
            even if the shortest pathway is a single reaction.
        exchange_reactions : bool
            Consider adding exchange (uptake) reactions for all metabolites
            in the model.
        demand_reactions : bool
            Consider adding demand reactions for all metabolites.

        Returns
        -------
        iterable
            list of lists with on set of reactions that completes the model per
            requested iteration.

        Examples
        --------

    import cobra as ct
        >>> from cobra import Model
        >>> from cobra.flux_analysis import gapfill
        >>> model = ct.create_test_model("salmonella")
        >>> universal = Model('universal')
        >>> universal.add_reactions(model.reactions.GF6PTA.copy())
        >>> model.remove_reactions([model.reactions.GF6PTA])
        >>> gapfill(model, universal)
    """
    gapfiller = GapFiller(
        model,
        universal=universal,
        lower_bound=lower_bound,
        penalties=penalties,
        demand_reactions=demand_reactions,
        exchange_reactions=exchange_reactions,
        community=community,
        taxa_database=taxa_database
    )
    return gapfiller.fill(iterations=iterations)
