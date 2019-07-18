#!/usr/bin/env python3
import argparse
import os

from cobra.flux_analysis import find_blocked_reactions
from cobra.io import read_sbml_model
from cobra.io import write_sbml_model

from csm4cobra.io import read_json
from csm4cobra.manipulate import set_medium


def create_parser():
    
    parser = argparse.ArgumentParser(description="Trim blocked reactions and gap metabolites from \
                                                  a genome-scale metabolic model")
    
    parser.add_argument('sbml_fname', action="store", help='SBML file to use a the model reference')
    
    parser.add_argument('--out', action="store", dest="sbml_fname_out", required=True,
                         help='SBML file name to for the outputed model')
    
    parser.add_argument('--media', action="store", dest="json_exchanges", default=None,
                         help='JSON file storing the exchange bounds')
    
    parser.add_argument('--open-exchanges', action="store_true", dest="open_exchanges",
                         help="A flag to indicade wheather to relax exchange fluxes bounds. \
                               Ignored if --media is also used")
    
    parser.add_argument('--exchange-prefix', action="store", dest="exchange_prefix", default="EX",
                         help='Prefix for the exchange reaction. Use with open-exhanges')
                         
    parser.add_argument('--flux-bound', action="store", dest="flux_bound", default=1000.,
                         help='Prefix for the exchange reaction. Use with open-exhanges')
                         
    parser.add_argument('--usefbc', action="store_true", help='Write SBML files using FBC package')
    

    
    return parser


def main():

    parser = create_parser()
    args = parser.parse_args()

    assert os.path.isfile(args.sbml_fname)
        
    # Reading SBML genome-scale model
    print("Reading SBML Model from %s:" % args.sbml_fname, end=" ")
    model = read_sbml_model(args.sbml_fname)
    print("OK!")
    
    if args.json_exchanges:
        print("Reading exchange fluxes bounds: %s:" % args.json_exchanges, end=" ")
        media_dict = read_json(args.json_exchanges)
        print("OK!")
        
        print("Setting exchange fluxes bounds")
        set_medium(model, media_dict, inplace=True)
        print("OK!")
    else:
        if args.open_exchanges:
            for r in model.reactions:
                if not r.id.startswith(args.exchange_prefix):
                    continue
                r.lower_bound = -args.flux_bound
                r.upper_bound = args.flux_bound

    print("Finding blocked reactions and gap metabolites:", end=" ")
    blocked = find_blocked_reactions(model)
    blocked = set(blocked)
    gap_metabolites = [m for m in model.metabolites
                        if len(set([r.id for r in m.reactions]) - blocked) == 0]
    print("OK!")
    
    if len(blocked) > 0:
        print("- %i blocked reactions found" % len(blocked))
        print("- %i gap metabolites found" % len(gap_metabolites))
        print("Trimming model", end=" ")
        model.remove_reactions(blocked)
        model.remove_metabolites(gap_metabolites)
        print("OK!")
        print("Writing trimmed model as %s" % args.sbml_fname_out, end=" ")
        write_sbml_model(model, args.sbml_fname_out, use_fbc_package=args.usefbc)
        print("OK!")
    else:
        print("NO blocked reactions found, nothing to do")
    

main()
