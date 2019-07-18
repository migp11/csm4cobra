#!/usr/bin/env python3
import sys
import csv
import cobra
import cobra.io


def main():
    sbml_file = sys.argv[1]
    model = cobra.io.read_sbml_model(sbml_file)

    reactions_file = "%s_reactions.csv" % sbml_file
    metabolites_file = "%s_metabolites.csv" % sbml_file

    f = open(reactions_file,'w')
    writer = csv.writer(f, delimiter='\t')
    for r in model.reactions:
        gene = r.gene_reaction_rule
        writer.writerow([r.id,r.name,r.subsystem,r.reaction,gene,r.lower_bound,r.upper_bound])
    f.close()

    f = open(metabolites_file,'w')
    writer = csv.writer(f, delimiter='\t')
    for m in model.metabolites:
        reactions = m.reactions
        writer.writerow([m.id,'_'.join(m.name.split('_')[0:-1]),m.name.split('_')[-1],' '.join([r.id for r in reactions])]) 
    f.close()


main()
