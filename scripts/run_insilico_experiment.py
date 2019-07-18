#!/usr/bin/env python3
import argparse
import os

import pandas as pd

from cobra.io import read_sbml_model
from cobra.util import solvers
from cobra.flux_analysis import pfba
from cobra.flux_analysis import single_gene_deletion
from cobra.flux_analysis import double_gene_deletion

from csm4cobra.io import read_json
from csm4cobra.manipulate import set_medium
from csm4cobra.manipulate import add_atpm_reaction
from csm4cobra.consistency_analysis import UmFinder


EXPERIMENTS = ['test_biomass', 'test_atp_yields', 'test_met_func', 'fba', 'single_gene_ko',
               'double_gene_ko', 'minimal_media', 'consistency_analysis']

SOLVERS = list(solvers.keys())
if 'cglpk' in SOLVERS:
    SOLVERS.remove('cglpk')
    SOLVERS.append('glpk')


def create_parser():
    parser = argparse.ArgumentParser(description='Run an in-silico experiment genome-scale metabolic model.')
    parser.add_argument('sbml_fname', action="store", help='SBML file to use a the model reference')
    
    parser.add_argument('--media', action="store", dest="json_exchanges", default=None,
                        help='JSON file storing the exchange bounds')

    parser.add_argument('--experiment', action="store", dest="experiment", choices=EXPERIMENTS,
                        required=True, help='Choose the in-silico experiment to be performed')

    parser.add_argument('--json-yields', action="store", dest="json_yields",
                        help='Path to json file with ATP yield reference')

    parser.add_argument('--json-met-func', action="store", dest="json_met_func",
                        help='Path to json file with metabolic functions')

    parser.add_argument('--solver', action="store", dest="solver", choices=SOLVERS,
                        default='glpk', help='LP solver to perform optimizations')

    parser.add_argument('--out', action="store", dest="output_folder", default=".",
                        help='Output folder to store the builded CSM')

    parser.add_argument('--processes', action="store", dest="processes", default=1, type=int,
                        help='Number of processesors used for computation')

    return parser


def generate_output_fname(output_folder, model_id, experiment, solver, output_format='tsv'):
    fname = "_".join([model_id, experiment, solver])
    fname = ".".join([fname, output_format])
    fname = os.path.join(output_folder, fname)
    return fname


def test_biomass_production(model, biomass_id='biomass_reaction'):
    # PREPARING MODEL FOR COMPONENTS PRODUCTION TEST
    model = model.copy()
    if 'ATPM' not in model.reactions:
        add_atpm_reaction(model)

    biomass_reaction = model.reactions.get_by_id(biomass_id)
    for r in model.reactions:
        if r.objective_coefficient == 0:
            continue
        r.objective_coefficient = 0

    biomass_components = dict()
    biomass_components['Energy'] = model.reactions.ATPM
    biomass_components['Biomass'] = model.reactions.biomass_reaction
    for m in biomass_reaction.reactants:
        reaction_list = list(m.reactions)
        reaction_list.remove(biomass_reaction)
        biomass_components[m.id] = reaction_list.pop()
        model.add_boundary(m, type='demand', lb=0.0)

    # The individual biomass components:
    # biomass_DNA_c, biomass_RNA_c, biomass_carbohydrate_c,
    # biomass_lipid_c, biomass_protein_c, biomass_other_c
    # Energy (ATPM) & biomass_reaction
    index = sorted(biomass_components.keys())
    df_result = pd.DataFrame(0.0, index=index, columns=['fOpt'])
    df_result.index.name = 'Component'
    for component, demand_rxn in biomass_components.items():
        demand_rxn.objective_coefficient = 1
        solution = model.optimize()
        demand_rxn.objective_coefficient = 0
        df_result.at[component, 'fOpt'] = round(solution.f, 4)

    return df_result


def test_atp_yields(model, atp_yields, zero_cutoff=1e-7, oxygen_bound=20):
    model = model.copy()
    if 'ATPM' not in model.reactions:
        add_atpm_reaction(model)

    for r in model.reactions:
        if r.objective_coefficient == 0:
            continue
        r.objective_coefficient = 0
    
    model.reactions.ATPM.objective_coefficient = 1

    index = sorted(atp_yields.keys())
    columns = ['Predicted', 'Reference', 'Valid']
    df_result = pd.DataFrame(index=index, columns=columns)
    df_result.index.name = 'rxn_id'

    for aminoacid_exchange, amino_acid_atp_yield in atp_yields.items():
        the_medium = {aminoacid_exchange: 1, 'EX_o2_e': oxygen_bound}
        set_medium(model, the_medium, inplace=True)
        solution = pfba(model)

        valid = abs(solution.fluxes.ATPM - amino_acid_atp_yield) < zero_cutoff
        df_result.at[aminoacid_exchange, 'Predicted'] = solution.fluxes.ATPM
        df_result.at[aminoacid_exchange, 'Reference'] = amino_acid_atp_yield
        df_result.at[aminoacid_exchange, 'Valid'] = valid

    return df_result


def test_metabolic_functions(model, met_func, zero_cutoff=1e-7):
    model = model.copy()
    met_func_dict = {}
    for m in met_func:
        if m in model.metabolites:
            try:
                met = model.metabolites.get_by_id(m)
                met_func_dict[m] = model.add_boundary(met, type='demand')
            except:
                print("Demand reactions already present for %s" % m)
                demand_id = "DM_" + m
                met_func_dict[m] = model.reactions.get_by_id(demand_id)

        else:
            met_func_dict[m] = None
    
    for r in model.reactions:
        if r.objective_coefficient == 0:
            continue
        r.objective_coefficient = 0
    
    index = sorted(met_func_dict.keys())
    columns = ['fluxes', 'produced']
    df_result = pd.DataFrame(index=index, columns=columns)
    df_result.index.name = "met_id"
    
    for idx, rxn in met_func_dict.items():
        if rxn:
            rxn.objective_coefficient = 1
            sol = model.optimize()
            df_result.at[idx, 'fluxes'] = sol.f
            df_result.at[idx, 'produced'] = sol.f > zero_cutoff
            rxn.objective_coefficient = 0
        else:
            df_result.at[idx, 'fluxes'] = 0
            df_result.at[idx, 'produced'] = False
    
    return df_result
    

def run_insilico_experiment(model, args, rescale=False):

    df_results = None

    if args.experiment == 'test_atp_yields':
        assert args.json_yields
        yields = read_json(args.json_yields)
        df_results = test_atp_yields(model, yields)
    
    if args.experiment == 'test_met_func':
        print("Reading Metabolic Functions: %s:" % args.json_met_func, end=" ")
        met_func = read_json(args.json_met_func)
        print("OK!")
        df_results = test_metabolic_functions(model, met_func)

    if args.experiment == 'test_biomass':
        df_results = test_biomass_production(model)

    if args.experiment == 'fba':
        solution = pfba(model)
        df_results = solution.to_frame()
        df_results.index.name = 'rxn_id'

    if args.experiment == 'single_gene_ko':
        if rescale:
            for r in model.exchanges:
                r.lower_bound *= 10
                r.upper_bound *= 10
        df_results = single_gene_deletion(model, processes=args.processes)

    if args.experiment == 'double_gene_ko':
        if rescale:
            for r in model.exchanges:
                r.lower_bound *= 10
                r.upper_bound *= 10
        
        genes = model.genes
        df_results = double_gene_deletion(model, gene_list1=genes,
                                          gene_list2=genes, processes=args.processes)

    if args.experiment == 'minimal_media':
        raise NotImplementedError
        # df_results = find_minimal_media(model)

    if args.experiment == 'consistency_analysis':
        gap_finder = UmFinder(model)
        blk_rxns = gap_finder.blocked_reactions
        gap_mets = gap_finder.gap_metabolites
        df_results = gap_finder.unconnected_modules_frame()

    return df_results


def main():

    parser = create_parser()
    args = parser.parse_args()

    assert os.path.isfile(args.sbml_fname)
    assert os.path.isdir(args.output_folder)
    assert os.path.isdir(args.output_folder)

    # Reading SBML genome-scale model
    print("Reading SBML Model from %s:" % args.sbml_fname, end=" ")
    model = read_sbml_model(args.sbml_fname)
    model.solver = args.solver
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
    df_result = run_insilico_experiment(model, args)
    print("OK!")

    tsv_output = generate_output_fname(args.output_folder, model.id,
                                       args.experiment, args.solver)
    
    print("Writing experiment results to %s:" % tsv_output, end=" ")
    df_result.to_csv(tsv_output, sep='\t')
    print("OK!")


main()
