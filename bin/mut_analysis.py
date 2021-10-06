#!/usr/bin/env python3
"""
Parralelised analysis of complex mutants using FoldX AnalyseComplex.
The regular FoldX method of analysing multiple mutant PDBs appears to gradually
use more and more RAM and uses no parralelisation.
"""
import os
import multiprocessing
import subprocess
import itertools
import argparse

def run_analyse_complex(pdb, pdb_dir, interface, output_dir):
    """
    Run FoldX AnalyseComplex command on a single PDB
    """
    command = ['foldx', '--command=AnalyseComplex', f'--pdb={pdb}',
               f'--pdb-dir={pdb_dir}', '--clean-mode=3',
               f'--output-dir={output_dir}',
               f'--analyseComplexChains={interface}']
    subprocess.run(command, capture_output=True)
    print(pdb, 'done', flush=True)

def main(args):
    """
    Run FoldX AnalyseComplex on each input mutant PDB
    """
    if not os.path.isdir(args.output):
        os.mkdir(args.output)

    mutant_pdbs = [i for i in os.listdir(args.pdb) if i.endswith('.pdb')]
    total = len(mutant_pdbs)

    with multiprocessing.Pool(processes=args.processes) as pool:
        print('Opened worker pool with', pool._processes, 'workers', flush=True)
        args = zip(mutant_pdbs,
                   itertools.cycle([args.pdb]),
                   itertools.cycle([args.interface.replace('_', ',')]),
                   itertools.cycle([args.output]))
        pool.starmap(run_analyse_complex, args)

def parse_args():
    """Process arguments"""
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('pdb', metavar='P', help="Directory of Mutant PDB files to process")
    parser.add_argument('interface', metavar='I', help="FoldX interface to process")
    parser.add_argument('output', metavar='O', help="Directory to output results")
    parser.add_argument('--processes', '-p', default=1, type=int,
                        help="Number of processes available")

    return parser.parse_args()

if __name__ == "__main__":
    main(parse_args())
