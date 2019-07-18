#!/usr/bin/env python3
import argparse
import os

import numpy as np
import pandas as pd
from cobra.io import read_sbml_model
from cobra.util import solvers

from csm4cobra.io import read_json
from csm4cobra.manipulate import set_medium
from csm4cobra.context_specific_deletions import get_all_gene_ko_reactions
from csm4cobra.context_specific_deletions import get_gene_knockout_reactions
from csm4cobra.context_specific_deletions import context_specific_ko


EXPERIMENTS = ['cs_gene_ko', 'double_gene_ko']

SOLVERS = list(solvers.keys())
if 'cglpk' in SOLVERS:
    SOLVERS.remove('cglpk')
    SOLVERS.append('glpk')


def create_parser():
    parser = argparse.ArgumentParser(description='Run an in-silico experiment genome-scale metabolic model.')
    parser.add_argument('sbml_fname', action="store", help='SBML file to use a the model reference')
    parser.add_argument('csv_confidences', action="store", help='CSV file storing the gene confidences')
    parser.add_argument('experiment', action="store", choices=EXPERIMENTS,
                         help='Choose the in-silico experiment to be performed')
    
    parser.add_argument('--ceres', action="store", dest="csv_ceres", default=None,
                           help='CSV file storing gene ceres score')

    parser.add_argument('--media', action="store", dest="json_exchanges", default=None,
                            help='JSON file storing the exchange bounds')

    parser.add_argument('--solver', action="store", dest="solver", choices=SOLVERS,
                         default='glpk', help='LP solver to perform optimizations')

    parser.add_argument('--out', action="store", dest="output_folder", default=".",
                         help='Output folder to store the builded CSM')

    parser.add_argument('--column-name', action="store", dest="col_name", default="confidence",
                        help='Column name where the RPKM values are stored')

    return parser


def generate_output_fname(output_folder, model_id, experiment,
                          solver, output_format='tsv'):
    fname = "_".join([model_id, experiment,solver])
    fname = ".".join([fname, output_format])
    fname = os.path.join(output_folder,fname)
    return fname


def run_context_specific_ko(model, conf_genes, threshold=2,
                            objectives=('biomass_reaction', 'ATPM')):
    result_dict = {}
    model_genes = {g.id for g in model.genes}
    objectives = [model.reactions.get_by_id(r) for r in objectives]
    for r in model.reactions:
        if r.objective_coefficient == 0:
            continue
        r.objective_coefficient = 0

    for obj_rxn in objectives:
        result_dict[obj_rxn.id] = {}
        obj_rxn.objective_coefficient = 1
        solution = model.optimize()
        result_dict[obj_rxn.id]['wild_type'] = solution.objective_value
        obj_rxn.objective_coefficient = 0

    gene_ko_reactions = get_gene_knockout_reactions(model)

    all_genes_ko_reactions = get_all_gene_ko_reactions(model, conf_genes, threshold=threshold)

    # Set of genes identified by contextualized gpr evaluation
    cs_gene_ko = set(all_genes_ko_reactions.keys()) - set(gene_ko_reactions.keys())
    for gene, rxn_list in all_genes_ko_reactions.items():
        bounds_dict = {}
        for rxn in rxn_list:
            bounds_dict[rxn.id] = rxn.bounds
            rxn.bounds = (0, 0)

        for obj_rxn in objectives:
            obj_rxn.objective_coefficient = 1
            solution = model.optimize()
            result_dict[obj_rxn.id][gene] = solution.objective_value
            obj_rxn.objective_coefficient = 0

        for rxn in rxn_list:
            rxn.bounds = bounds_dict[rxn.id]

    result_dict['confidence'] = {}
    result_dict['in_model'] = {}
    result_dict['is_cs_ko'] = {}
    result_dict['inactivate_reactions'] = {}
    genes_inactivate_reactions = set(all_genes_ko_reactions.keys())
    for gene, conf in conf_genes.items():
        result_dict['confidence'][gene] = conf
        result_dict['is_cs_ko'][gene] = False
        if gene in genes_inactivate_reactions:
            result_dict['inactivate_reactions'][gene] = len(all_genes_ko_reactions[gene])
            if gene in cs_gene_ko:
                result_dict['is_cs_ko'][gene] = True
        else:
            result_dict['inactivate_reactions'][gene] = 0
            result_dict['biomass_reaction'][gene] = result_dict['biomass_reaction']['wild_type']
            result_dict['ATPM'][gene] = result_dict['ATPM']['wild_type']

        if gene in model_genes:
            result_dict['in_model'][gene] = True
        else:
            result_dict['in_model'][gene] = False

    df_results = pd.DataFrame(result_dict)
    df_results.index.name = 'gene_id'

    return df_results

def run_insilico_ko(model, genes_ko_reactions, objectives=('biomass_reaction', 'ATPM')):

    result_dict = {'confidence': {},
                   'in_model': {},
                   'is_cs_ko': {},
                   'inactivate_reactions': {}
                   }

    model_genes = {g.id for g in model.genes}

    objectives = [model.reactions.get_by_id(r) for r in objectives]
    for r in model.reactions:
        if r.objective_coefficient == 0:
            continue
        r.objective_coefficient = 0

    result_dict['is_cs_ko']['wild_type'] = False
    result_dict['in_model']['wild_type'] = False
    for obj_rxn in objectives:
        result_dict[obj_rxn.id] = {}
        obj_rxn.objective_coefficient = 1
        solution = model.optimize()
        result_dict[obj_rxn.id]['wild_type'] = solution.f
        obj_rxn.objective_coefficient = 0

    for gene, rxn_list in genes_ko_reactions.items():
        bounds_dict = {}
        for rxn in rxn_list:
            bounds_dict[rxn.id] = rxn.bounds
            rxn.bounds = (0, 0)

        for obj_rxn in objectives:
            obj_rxn.objective_coefficient = 1
            solution = model.optimize()
            result_dict[obj_rxn.id][gene] = solution.f
            obj_rxn.objective_coefficient = 0

        for rxn in rxn_list:
            rxn.bounds = bounds_dict[rxn.id]

    genes_inactivate_reactions = set(genes_ko_reactions.keys())
    for gene in model_genes:
        result_dict['confidence'] = 0
        result_dict['is_cs_ko'][gene] = False
        result_dict['in_model'][gene] = True

        if gene in genes_inactivate_reactions:
            result_dict['inactivate_reactions'][gene] = len(genes_ko_reactions[gene])
        else:
            result_dict['inactivate_reactions'][gene] = 0
            result_dict['biomass_reaction'][gene] = result_dict['biomass_reaction']['wild_type']
            result_dict['ATPM'][gene] = result_dict['ATPM']['wild_type']
            result_dict['is_cs_ko']['wild_type'] = False

    df_results = pd.DataFrame(result_dict)
    df_results.index.name = 'gene_id'

    return df_results


def run_experiment(model, experiment, conf_genes_dict, ):
    df_results = None

    if experiment == 'cs_gene_ko':
        df_results = context_specific_ko(model, conf_genes_dict)

    if experiment == 'double_gene_ko':
        raise NotImplementedError('Function not implemented')

    return df_results


def main():

    parser = create_parser()
    args = parser.parse_args()

    assert os.path.isfile(args.sbml_fname)
    assert os.path.isdir(args.output_folder)

    # Reading SBML genome-scale model
    print("Reading SBML Model from %s:" % args.sbml_fname, end=" ")
    model = read_sbml_model(args.sbml_fname)
    model.solver = args.solver
    print("OK!")

    print("Reading gene confidences: %s:" % args.csv_confidences, end=" ")
    # df_genes_conf = read_genes_table(args.csv_conf)
    # conf_genes_dict = df_genes_conf[args.col_name].to_dict()
    df_conf_genes = pd.read_csv(args.csv_confidences, sep='\t', index_col='gene_id')
    conf_genes_dict = df_conf_genes.confidence.to_dict()
    print("OK!")

    if args.json_exchanges:
        assert os.path.isfile(args.json_exchanges)
        print("Reading exchange fluxes bounds: %s:" % args.json_exchanges, end=" ")
        media_dict = read_json(args.json_exchanges)
        print("OK!")

        print("Setting exchange fluxes bounds:", end=" ")
        set_medium(model, media_dict, inplace=True)
        print("OK!")

    print("Running in-silico experiment: \"%s\":" % args.experiment)
    df_results = run_experiment(model, args.experiment, conf_genes_dict)
    print("OK!")

    if args.csv_ceres:
        assert os.path.isfile(args.csv_ceres)
        print("Reading genes CERES: %s:" % args.csv_ceres, end=" ")
        df_ceres = pd.read_csv(args.csv_ceres, sep='\t', index_col='gene_id')
        values = [df_ceres.ceres[g] if g in df_ceres.index else np.nan for g in df_results.index]
        df_results = df_results.assign(ceres=values)
        print("OK!")

    tsv_output = generate_output_fname(args.output_folder, model.id,
                                       args.experiment, args.solver)
    
    print("Writing experiment results to %s:" % tsv_output, end=" ")
    df_results.index.name = 'gene_id'
    df_results.to_csv(tsv_output, na_rep='NA', sep='\t')
    print("OK!")

main()
