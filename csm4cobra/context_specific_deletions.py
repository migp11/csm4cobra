#  the function format_gid, safe_eval_gpr and reaction_confidence were extracted
#  from corda.util.py modules
#
#  https://github.com/resendislab/corda/
#
#  Copyright 2016 Christian Diener <mail[at]cdiener.com>
#
#  MIT license. See LICENSE for more information.
#
################################################################################################

import re
import pandas as pd
from cobra.core.gene import parse_gpr
from ast import Name, And, Or, BoolOp, Expression
from cobra.flux_analysis.deletion import find_gene_knockout_reactions


def format_gid(gid):
    """Internal function to strip transcript dot-notation from IDs."""
    return re.sub(r"\.\d*", "", gid)


def safe_eval_gpr(expr, conf_genes):
    """Internal function to evaluate a gene-protein rule in an
    injection-safe manner (hopefully).
    """
    if isinstance(expr, Expression):
        return safe_eval_gpr(expr.body, conf_genes)
    elif isinstance(expr, Name):
        fgid = format_gid(expr.id)
        if fgid not in conf_genes:
            return 0
        return conf_genes[fgid]
    elif isinstance(expr, BoolOp):
        op = expr.op
        if isinstance(op, Or):
            return max(safe_eval_gpr(i, conf_genes) for i in expr.values)
        elif isinstance(op, And):
            return min(safe_eval_gpr(i, conf_genes) for i in expr.values)
        else:
            raise TypeError("unsupported operation " + op.__class__.__name__)
    elif expr is None:
        return 0
    else:
        raise TypeError("unsupported operation  " + repr(expr))


def reaction_confidence(rule, conf_genes):
    """Calculates the confidence for the reaction based on a gene-reaction
    rule.
    Args:
        rule (str): A gene-reaction rule. For instance "A and B".
        conf_genes (dict): A str->int map denoting the mapping of gene IDs
            to expression confidence values. Allowed confidence values are -1
            (absent/do not include), 0 (unknown), 1 (low confidence),
            2 (medium confidence) and 3 (high confidence).
    """
    ast_rule, _ = parse_gpr(rule)
    return safe_eval_gpr(ast_rule, conf_genes)


################################################################################################

def get_gene_knockout_reactions(model):
    genes_ko_reactions = {}
    for g in model.genes:
        reactions = find_gene_knockout_reactions(model, [g])
        if len(reactions) == 0:
            continue
        genes_ko_reactions[g.id] = reactions

    return genes_ko_reactions


def contextualized_ko_effect(reaction, conf_genes, ko_confidence=-1):
    """ Given a reaction and a dictionary with its associated genes confidences this function
    compute the change in the reaction confidences when each of

    Parameters
    ----------

    reaction: A cobra.core.Reaction object
    conf_genes: (dict): A mapping of all genes ids in `reaction`
        to an integer denoting the gene expression confidence.
        Allowed confidence values are -1 (absent/do not include),
        0 (unknown), 1 (low confidence), 2 (medium confidence)
        and 3 (high confidence). See `reaction_confidence` for a way
        to construct this dictionary.
    ko_confidence: the confidence value set to a knocked gene

    Return
    ----------
    dict
        A dictionary os str->int mapping the genes IDs in the reaction
        to the reaction confidence value obtained when each gene is knocked
        out (ie. its confidences is set to -1)
    """

    cs_ko_effect = dict()
    cs_ko_effect['wild_type'] = reaction_confidence(reaction.gene_reaction_rule, conf_genes)
    for g in reaction.genes:
        # Storing the original confidences value for the gene
        gene_conf = conf_genes[g.id]
        # If the gene confidence is already the lowest, the it will not change the reaction confidence
        if gene_conf == ko_confidence:
            cs_ko_effect[g.id] = cs_ko_effect['wild_type']
        else:
            # Setting the for the given gene confidence to -1 to simulate its ko
            conf_genes[g.id] = ko_confidence
            # recompute the reaction confidence for the new context
            cs_ko_effect[g.id] = reaction_confidence(reaction.gene_reaction_rule, conf_genes)

        # Setting back the original gene confidence value
        conf_genes[g.id] = gene_conf

    return cs_ko_effect


def get_cs_gene_ko_reactions(model, conf_genes, threshold=0.5):
    """ This function identify the set of reactions inactivated by each gene
    when the expression confidences are considered, this is, any gene whose
    inactivation results in a drop of the reaction confidence.

    Parameters
    ----------
    model: cobra.core.Model
    conf_genes: a python dict with gene ids as keys and expression or confidences as values.
    threshold: a float: 0 <= threshold <= 1 tha will be used as a fraction fraction
        above which a reaction is considered inactive when a gene is knockout

    Returns
    -------
    dict()
        A dictionary that map each gene into the reactions inactivated by its knockout.

    Notes
    -----
    A Reaction confidence is computed combining its GPR and the expression
    confidences of the involved gene. To asses the impact of knocking a gene,
    its confidence is set to the lowest values (-1) and the reaction confidences
    is recalculated. If the new confidence value drop below a given threshold
    the gene knockout is considered to be sufficient to inactivate the reaction.
    """
    cs_gene_knockout_reactions = {}
    for rxn in model.reactions:
        # Skipping reactions without coding genes (orphans, exchanges, etc)
        if len(rxn.genes) == 0:
            continue

        rxn_conf_genes = {g.id: conf_genes[g.id] if g.id in conf_genes else 0 for g in rxn.genes}

        cs_ko_effect = contextualized_ko_effect(rxn, rxn_conf_genes)

        wt_conf = cs_ko_effect['wild_type']
        if wt_conf == 0:
            continue

        del cs_ko_effect['wild_type']

        for gene, ko_conf in cs_ko_effect.items():
            # Here the reaction confidence is compared to the confidences
            # previous to knocking any of its associated genes. If the
            # confidence after the gene knockout drop above the threshold,
            # then the knockout of that gene is considered to inactivate the reaction
            if ko_conf/wt_conf < threshold:

                if gene not in cs_gene_knockout_reactions:
                    cs_gene_knockout_reactions[gene] = []

                cs_gene_knockout_reactions[gene].append(rxn)

    return cs_gene_knockout_reactions


def get_all_gene_ko_reactions(model, conf_genes, threshold=0.5):
    """ This function identify all the genes that inactivate one or more reactions
        The set include genes tha inactivate reactions because of the GPR as well
        genes that inactivate reactions when the expression confidences of all the
        genes are considered.

        Parameters
        ----------
        model: cobra.core.Model
        conf_genes: a python dict with gene ids as keys and expression or confidences as values.
        threshold: a float: 0 <= threshold <= 1 tha will be used as a fraction fraction
            above which a reaction is considered inactive when a gene is knockout

        Returns
        -------
        dict()
            A dictionary that map each gene into the reactions inactivated by its knockout.

        Notes
        -----
        A Reaction confidence is computed combining its GPR and the expression
        confidences of the involved gene. To asses the impact of knocking a gene,
        its confidence is set to the lowest values (-1) and the reaction confidences
        is recalculated. If the new confidence value drop below a given threshold
        the gene knockout is considered to be sufficient to inactivate the reaction.
        See get_cs_gene_ko_reactions
        """
    all_gene_ko_reactions = get_cs_gene_ko_reactions(model, conf_genes, threshold=threshold)

    gene_ko_reactions = get_gene_knockout_reactions(model)

    for g in gene_ko_reactions:
        if g not in all_gene_ko_reactions:
            all_gene_ko_reactions[g] = gene_ko_reactions[g]
        else:
            for rxn in gene_ko_reactions[g]:
                if rxn in all_gene_ko_reactions[g]:
                    continue
                all_gene_ko_reactions[g].append(rxn)

    return all_gene_ko_reactions


def context_specific_ko(model, genes_confidences, threshold=0.5):

    gene_ko_reactions = get_gene_knockout_reactions(model)
    all_genes_ko_reactions = get_all_gene_ko_reactions(model, genes_confidences, threshold=threshold)

    # Set of genes identified by contextualized gpr evaluation
    cs_gene_ko = set(all_genes_ko_reactions.keys()) - set(gene_ko_reactions.keys())

    index = []
    data = []

    wt_solution = model.optimize()

    for gene in model.genes:
        index.append(gene)
        rxn_list = all_genes_ko_reactions[gene.id]
        # The gene do not deactivate any reaction
        if len(rxn_list) == 0:
            data.append((wt_solution.objective_value, wt_solution.status, False))

        else:

            bounds_dict = {}
            for rxn in rxn_list:
                bounds_dict[rxn.id] = rxn.bounds
                rxn.bounds = (0, 0)

            solution = model.optimize()
            is_cs_ko = gene in cs_gene_ko
            data.append((solution.objective_value, solution.status, is_cs_ko))

            for rxn in rxn_list:
                rxn.bounds = bounds_dict[rxn.id]

    df_results = pd.DataFrame(data=data, index=index, columns=['growth', 'status', 'is_cs_ko'])
    df_results.index.name = 'gene_id'

    return df_results
