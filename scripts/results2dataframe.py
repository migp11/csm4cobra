#!/usr/bin/env python3
import os
import glob
import argparse
import pandas as pd


MEM_METHODS = ['corda', 'fastcore']

EXPERIMENTS = ['test_biomass', 'test_atp_yields', 'test_met_func', 'test_fba', 'test_single_gene_ko',
               'test_double_gene_ko', 'test_minimal_media', 'test_consistency', 'binary']


def create_parser():
    parser = argparse.ArgumentParser(description="Join simulation results for different CSM in a single dataframe")
    
    parser.add_argument("--data", action="store", dest="data_folder", default=None,
                        help="Folder where the data in TSV forma is stored")
    
    parser.add_argument("--experiment", action="store", dest="experiment", choices=EXPERIMENTS, required=True,
                        help='Choose the in-silico experiment to be performed')
    
    parser.add_argument("--out", action="store", dest="fname_out", default=".",
                        help="TSV Output file name")
                            
    parser.add_argument("--method", action="store", dest="method", choices=MEM_METHODS, default="corda",
                        help="Method used to build the CSM.")
                            
    parser.add_argument("--idx", action="store", dest="idx_name", default='hgnc_id',
                        help="Colunm name to use as an index")
                            
    parser.add_argument("--col", action="store", dest="col_name", required=True, 
                        help="Colunm name to get the data from")
    
    return parser


def main():
    
    parser = create_parser()
    args = parser.parse_args()
 
    data_dict = {}
    files_globbing = r"*%s*%s*.tsv" % (args.method, args.experiment)
    files_globbing = os.path.join(args.data_folder, files_globbing)
    for fname in glob.glob(files_globbing):
        cell_line_id = os.path.basename(fname).split("_")[0]
        df_aux = pd.read_csv(fname, sep="\t", index_col=args.idx_name)
        data_dict[cell_line_id] = df_aux[args.col_name]

    columns = sorted({i for serie in data_dict.values() for i in serie.index})
    index = sorted(data_dict.keys())

    df_results = pd.DataFrame(index=index, columns=columns)
    for cl_id, serie in data_dict.items():
        for k, v in serie.items():
            df_results.at[cl_id, k] = v
            
    if args.experiment == 'test_single_gene_ko':
        wt_dict = df_results.max(axis=1).to_dict()
        for cl_id, serie in data_dict.items():
            for gene_id in df_results.columns:
                if gene_id in serie:
                    continue
                df_results.at[cl_id, gene_id] = wt_dict[cl_id]

    df_results.index.name = 'cell_line_id'
    df_results.to_csv(args.fname_out, sep="\t")


main()
