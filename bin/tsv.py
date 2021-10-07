#!/usr/bin/env python3
"""
Combine FoldX AnalyseComplex output from many complexes
"""
import sys
import argparse
import pandas as pd
from pathlib import Path

def import_complex_dir(path):
    """
    Import tables from an AnalyseComplex output directory
    """
    path = path.rstrip('/')
    interactions = pd.read_csv(f'{path}/interactions.tsv', sep='\t')
    interactions = interactions.rename({'interface_residues': 'number_of_interface_residues'},
                                       axis='columns')
    interface = pd.read_csv(f'{path}/interface_residues.tsv', sep='\t')
    comb = pd.merge(interactions, interface, how='outer', on=['chain', 'position', 'wt', 'mut'])

    comb['complex'] = path.split("/")[-2]
    comb['interface'] = path.split("/")[-1]
    cols = ['complex', 'interface', 'chain', 'position', 'wt', 'mut']
    comb = comb[cols + [c for c in comb.columns if not c in cols]]
    return comb

def main(args):
    """Main"""
    complex_dfs = [import_complex_dir(path) for path in args.dir]
    complexes = pd.concat(complex_dfs)
    sort_cols = ['complex', 'interface', 'chain', 'position', 'wt', 'mut']
    complexes = complexes.sort_values(axis='rows', by=sort_cols).reset_index(drop=True)
    complexes.to_csv(sys.stdout, sep='\t', index=False)

def parse_args():
    """Process input arguments"""
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('dir', metavar='D', nargs='+',
                        help="Directories containing the output of the AnalyseComplex pipeline")

    return parser.parse_args()

if __name__ == "__main__":
    main(parse_args())