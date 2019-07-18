#!/usr/bin/env python3
import argparse
import os
from os import path

import pandas as pd
from cobra.io import read_sbml_model, write_sbml_model
from cobra.util import solvers
from corda import CORDA
from corda.util import reaction_confidence
from fastcore import Fastcore

from csm4cobra.io import read_json
from csm4cobra.manipulate import set_medium

MEM_METHODS = ['corda', 'fastcore']
SOLVERS = list(solvers.keys())
if 'cglpk' in SOLVERS:
    SOLVERS.remove('cglpk')
    SOLVERS.append('glpk')


def create_parser():
    parser = argparse.ArgumentParser(description='Build a context specific genome-scale metabolic model.')
    parser.add_argument('sbml_fname', action="store", help='SBML file to use a the model reference')
    parser.add_argument('json_confidences', action="store", help='CSV file storing the gene confidences')
    parser.add_argument('model_id', action="store", help='Model id for the builded CSM')

    parser.add_argument('--gene-id', action="store", dest="gene_id", default="gene_symbol",
                        help='Gene identifier name used as column name')

    parser.add_argument('--media', action="store", dest="json_exchanges", default=None,
                        help='JSON file storing the exchange bounds')

    parser.add_argument('--out', action="store", dest="output_folder", default=".",
                        help='Output folder to store the builded CSM')

    parser.add_argument('--method', action='store', dest='method', choices=MEM_METHODS, default='corda',
                        help='Method used to build the CSM: corda / fastcore')

    parser.add_argument('--solver', action="store", dest="solver", choices=SOLVERS,
                        default='glpk', help='LP solver to perform optimizations')

    parser.add_argument('--force-biomass', action='store_true', dest='force_biomass',
                        help='Biomass reactions are set to have the highest confidence')

    parser.add_argument('--use-fbc', action='store', dest='use_fbc', default=True, type=bool,
                        help='Write SBML using FBC package')

    parser.add_argument('--write-confidences', action='store_true', dest='write_confidences',
                        help='Store computed reaction confidences in CSV format')

    curdir = path.dirname(path.realpath(__file__))
    fastcore_params_fname = path.join(curdir, "config/fastcore.json")
    parser.add_argument('--fastcore-params', action="store", dest="fastcore_params_fname",
                        default=fastcore_params_fname, help='JSON file with fastcore parameters')

    return parser


def create_reactions_confidences(model, gene_confidences, column='confidence'):
    # Creating a reactions confidence dict using gene confidendes (df_HPAC.CNS)
    # and the corda function reaction_confidence
    rxn_confidences = dict()
    for r in model.reactions:
        # Creating a confidence dict for each particular reaction
        conf_genes = {g.id: gene_confidences[g.id] for g in r.genes}
        rxn_confidences[r.id] = reaction_confidence(r.gene_reaction_rule, conf_genes)

    return rxn_confidences


def create_reactions_penalties(reactions_confidence, confidence_penalty_map):
    reactions_penalties = {}
    for rxn, conf in reactions_confidence.items():
        reactions_penalties[rxn] = confidence_penalty_map[conf]

    return reactions_penalties


def create_confidences(model, rxns_confidence):
    model_rxns = {r.id for r in model.reactions}
    df = pd.DataFrame(columns=['confidence', 'in_model'], index=rxns_confidence.keys())
    for r, conf in rxns_confidence.items():
        df.confidence[r] = conf
        df.in_model[r] = r in model_rxns

    df.index.name = 'ID'

    return df


def generate_model_id(args, sep="_"):
    if args.force_biomass:
        model_id = sep.join([args.model_id, args.method, args.solver, "biomass"])
    else:
        model_id = sep.join([args.model_id, args.method, args.solver])

    return model_id


def generate_output_fname(args, file_type):
    fname = generate_model_id(args)
    fname = ".".join([fname, file_type])
    fname = os.path.join(args.output_folder, fname)
    return fname


def create_fastcore(model, reactions_confidences, fastcore_params, check_consistency=False):

    core_rxn_cutoff = fastcore_params['core_rxn_cutoff']
    core_reactions = [rxn for rxn, conf in reactions_confidences.items() if conf >= core_rxn_cutoff]

    conf_penalty_map = {int(k): v for k, v in fastcore_params['conf_penalty_map'].items()}
    penalties = create_reactions_penalties(reactions_confidences, confidence_penalty_map=conf_penalty_map)

    for r in model.exchanges:
        penalties[r.id] = 0

    return Fastcore(model, core_reactions, penalties=penalties, check_consistency=check_consistency)



def main():
    parser = create_parser()
    args = parser.parse_args()
    biomass_reaction = "biomass_reaction"

    assert os.path.isfile(args.sbml_fname)
    assert os.path.isfile(args.json_confidences)
    assert os.path.isdir(args.output_folder)

    # Reading Reference Genome-Scale Model
    print("Reading SBML Model from %s:" % args.sbml_fname, end=" ")
    model = read_sbml_model(args.sbml_fname)
    print("OK!")
    # Setting optimization solver
    print("Setting optimization solver %s:" % args.solver, end=" ")
    model.solver = args.solver
    print("OK!")
    # Reading DataFrame including gene confidence
    print("Reading gene confidence from %s:" % args.json_confidences, end=" ")
    gene_confidence = read_json(args.json_confidences)
    # gene_confidence = pd.read_csv(args.csv_fname, delimiter='\t', index_col=args.gene_id)
    # gene_confidence = gene_confidence[args.confidence_column].to_dict()

    print("OK!")
    print("Computing reactions confidences:", end=" ")
    reaction_confidences = create_reactions_confidences(model, gene_confidence)
    print("OK!")

    # If True: biomass reactions are set to have the highest confidence
    if args.force_biomass:
        reaction_confidences[biomass_reaction] = 3
        print("Adding biomass reactions %s to the core set: OK!" % biomass_reaction)

    if args.json_exchanges:
        print("Reading exchange fluxes bounds: %s:" % args.json_exchanges, end=" ")
        media_dict = read_json(args.json_exchanges)
        print("OK!")
        print("Setting exchange fluxes bounds:", end=" ")
        set_medium(model, media_dict, inplace=True)
        print("OK!")

    print("Building context-specific model using \"%s\":" % args.method)
    if args.method == "corda":
        corda = CORDA(model, reaction_confidences)
        corda.build()
        cs_model = corda.cobra_model()

    elif args.method == "fastcore":

        print("Reading fastcore parameters: %s:" % args.fastcore_params_fname, end=" ")
        fastcore_params = read_json(args.fastcore_params_fname)
        print("OK!")
        fc_builder = create_fastcore(model, reaction_confidences, fastcore_params)
        fc_builder.fast_core()
        cs_model = fc_builder.build_context_specific_model()

    print("OK!")

    # Removing isolated metabolites
    isolated_metabolites = [m for m in cs_model.metabolites if len(m.reactions) == 0]
    cs_model.remove_metabolites(isolated_metabolites)

    isolated_genes = [g for g in cs_model.genes if len(g.reactions) == 0]
    for g in isolated_genes:
        cs_model.genes.remove(g)

    model_id = generate_model_id(args)
    cs_model.id = model_id
    cs_model.repair()


    sbml_output = generate_output_fname(args, 'xml')
    print("Writing context-specific model to %s:" % sbml_output, end=" ")
    # write_sbml_model(cs_model, sbml_output, use_fbc_package=args.use_fbc)
    write_sbml_model(cs_model, sbml_output)
    print("OK!")

    if args.write_confidences:
        csv_output = generate_output_fname(args, 'tsv')
        df_rxn_conf = create_confidences(cs_model, reaction_confidences)
        print("Writing reaction confidences to %s:" % csv_output, end=" ")
        df_rxn_conf.to_csv(csv_output, sep='\t')
        print("OK!")


main()