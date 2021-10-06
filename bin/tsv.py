#!/usr/bin/env python3
"""
Combine FoldX AnalyseComplex output from many complexes
"""
import sys
import argparse
import pandas as pd
from ruamel.yaml import YAML
from pathlib import Path

def import_complex_dir(path, chains, model):
    """
    Import tables from an AnalyseComplex output directory
    """
    path = path.rstrip('/')
    interactions = pd.read_csv(f'{path}/interactions.tsv', sep='\t')
    interactions = interactions.rename({'interface_residues': 'number_of_interface_residues'},
                                       axis='columns')
    interface = pd.read_csv(f'{path}/interface_residues.tsv', sep='\t')
    comb = pd.merge(interactions, interface, how='outer', on=['chain', 'position', 'wt', 'mut'])

    comb['uniprot'] = [chains[chain]['uniprot'] for chain in comb.chain]
    comb['name'] = [chains[chain]['name'] for chain in comb.chain]
    comb['model'] = model
    comb['int_chain'] = [chain1 if mut == chain2 else chain2 for mut, chain1, chain2 in
                         zip(comb.chain, comb.chain1, comb.chain2)]
    comb['int_uniprot'] = [chains[chain]['uniprot'] for chain in comb.int_chain]
    comb['int_name'] = [chains[chain]['name'] for chain in comb.int_chain]
    cols = ['uniprot', 'name', 'position', 'wt', 'mut',
            'int_uniprot','int_name', 'model', 'chain', 'int_chain']
    comb = comb[cols + [c for c in comb.columns if not c in cols]]
    comb = comb.drop(['chain1', 'chain2'], axis='columns')
    return comb

def main(args):
    """Main"""
    complex_dfs = []
    yaml_loader = YAML(typ='safe')
    for yaml in args.yaml:
        path = Path(yaml)
        yaml = yaml_loader.load(path)
        for interface in yaml['interfaces']:
            complex_dfs.append(import_complex_dir(f'{path.parent}/{interface}', yaml['chains'],
                                                  yaml['model']))

    complexes = pd.concat(complex_dfs)
    sort_cols = ['uniprot', 'name', 'position', 'mut', 'int_uniprot', 'int_name']
    complexes = complexes.sort_values(axis='rows', by=sort_cols).reset_index(drop=True)
    complexes.to_csv(sys.stdout, sep='\t', index=False)

def parse_args():
    """Process input arguments"""
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('yaml', metavar='Y', nargs='+',
                        help="YAML config files indicating location of each interface output")

    return parser.parse_args()

if __name__ == "__main__":
    main(parse_args())