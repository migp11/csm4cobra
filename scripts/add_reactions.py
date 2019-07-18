#!/usr/bin/env python3
import argparse

from cobra.io import read_sbml_model
from cobra.io import write_sbml_model

from csm4cobra.manipulate import add_reaction_dict
from csm4cobra.io import read_json


def create_parser():
    parser = argparse.ArgumentParser(description='Add reactions to a genome-scale metabolic model.')
    parser.add_argument('sbml_fname', action="store", help='SBML file to use a the model reference')
    parser.add_argument('json_fname', action="store", help='JSON file with the reactions to be added')

    parser.add_argument('--replace', action="store_true", help='If given, replace reaction if already present in model')
    parser.add_argument('--usefbc', action="store_true", help='Write SBML files using FBC package')
    
    parser.add_argument('--out', action="store", dest="sbml_out_fname", required=True,
                         help='SBML file name to for the outputed model')
    
    return parser


def main():
    
    parser = create_parser()
    args = parser.parse_args()
    
    model = read_sbml_model(args.sbml_fname)
    reactions_dict = read_json(args.json_fname)
    
    for rid,rxn_dict in reactions_dict.items():
        add_reaction_dict(model, rid,rxn_dict, compartment="", replace=args.replace)

    write_sbml_model(model, args.sbml_out_fname, use_fbc_package=args.usefbc)
    

main()
