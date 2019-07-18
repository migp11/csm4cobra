import warnings
from cobra.core import Reaction
from cobra.core import Metabolite
from cobra.flux_analysis.deletion import find_gene_knockout_reactions


def set_medium(model, medium, skip_absent=True, inplace=False):
    if not inplace:
        model = model.copy()

    def set_active_bound(reaction, bound):
        if reaction.reactants:
            reaction.lower_bound = -bound
        elif reaction.products:
            reaction.upper_bound = bound

    # Set the given media bounds
    media_rxns = list()
    for rxn_id, bound in medium.items():
        if rxn_id not in model.reactions and skip_absent:
            warnings.warn("Exchange flux %s not found, skippied" % rxn_id)
            continue
        rxn = model.reactions.get_by_id(rxn_id)
        media_rxns.append(rxn)
        set_active_bound(rxn, bound)

    boundary_rxns = set([r for r in model.exchanges if r.id.startswith('EX')])
    media_rxns = set(media_rxns)

    # Turn off reactions not present in media
    for rxn in (boundary_rxns - media_rxns):
        set_active_bound(rxn, 0)

    if not inplace:
        return model


def add_reaction_dict(model, reaction_id, reaction_dict, compartment='c', lb=0.0, ub=1000, replace=False, inplace=True):
    
    if not inplace:
        model = model.copy()

    metabolites_to_add = {}
    for met_id, coeff in reaction_dict.items():
        if len(compartment) > 0:
            met_id += "_" + compartment
        
        try:
            metabolite = model.metabolites.get_by_id(met_id)
        except:
            warnings.warn("Metabolites %s no present in model, added as new" % met_id)
            metabolite = Metabolite(met_id)

        metabolites_to_add[metabolite] = coeff

    if reaction_id in model.reactions:
        if replace:
            model.remove_reactions([reaction_id])
        else:
            return
    
    reaction = Reaction(reaction_id)
    reaction.lower_bound = lb
    reaction.upper_bound = ub
    reaction.add_metabolites(metabolites_to_add)

    model.add_reaction(reaction)
    
    if not inplace:
        return model


def add_atpm_reaction(model, inplace=True):
    metabolites_dict = {"atp":-1,"h2o":-1,"adp":1,"pi":1,"h":1}
    rxn_id = 'ATPM'
    add_reaction_dict(model, rxn_id, metabolites_dict, inplace=inplace)

