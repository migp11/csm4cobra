#!/usr/bin/env python3
import argparse
import pandas as pd

# from csm4cobra.gene_expression import TISSUES
from csm4cobra.gene_expression import discretize_gene_expression
from csm4cobra.io import read_genes_table


def create_parser():

    parser = argparse.ArgumentParser(description="Convert gene expression (RPMK) to confidences [-1, 0, 1, 2 ,3].")
    parser.add_argument("tsv_sample", action="store", help="TSV file storing the gene expression (RPKM)")

    parser.add_argument("tsv_local_thresholds", action="store", help="TSV file storing the gene expression (RPKM)")

    parser.add_argument("tsv_out", action="store", help="TSV file to save the discertized expression")

    parser.add_argument('--lb_global', action="store", dest="lb_global", required=True,
                        type=float, help='Global lower bound threshold (float > 0)')

    parser.add_argument('--ub_global', action="store", dest="ub_global", required=True,
                        type=float, help='Global upper bound threshold (float > 0)')

    parser.add_argument('--index-name', action="store", dest="idx_name", default="hgnc_id",
                        help='Index name where the gene ids are stored')

    parser.add_argument('--column-name', action="store", dest="col_name", default="rpkm",
                        help='Column name where the RPKM values are stored')

    parser.add_argument('--gene-list', action="store", dest="gene_list", default=None,
                        help='A list of gene\'s ids in plain text (one per line) to discretize')

    # parser.add_argument("tissue", action="store", choices=TISSUES, help="Samples's tissue of origin")

    return parser


def main():
    parser = create_parser()
    args = parser.parse_args()

    df_sample_rpmk = read_genes_table(args.tsv_sample, index=args.idx_name)
    df_local_thresholds = read_genes_table(args.tsv_local_thresholds, index=args.idx_name)

    df_confidences = discretize_gene_expression(df_sample_rpmk, df_local_thresholds,
                                                args.lb_global, args.ub_global,
                                                gene_list=args.gene_list,
                                                col_name=args.col_name)

    df_confidences.to_csv(args.tsv_out, sep='\t')

main()