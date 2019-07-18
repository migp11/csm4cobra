import json
import pandas as pd


def read_json(json_fname):
    with open(json_fname) as fh:
        result = json.load(fh)
    return result


def read_genes_table(fname, index='hgnc_id'):
    df = pd.read_csv(fname, delimiter='\t')
    df.set_index(index, inplace=True)
    return df
