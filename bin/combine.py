#!/usr/bin/env python3
"""
Combine output files from batch FoldX AnalyseComplex command.
Files must be named as FoldX names them - e.g. Interaction_NAME(_Repair)(_N)_AC.fxout where N is the
PDB number in PDB list or is absent for the WT file. They type of file to process (interfaces,
individual energy, interface residues or summary) is detected from the WT filename.
"""
import sys
import os
import re
import argparse
from dataclasses import dataclass
from pathlib import Path

INTERFACE_RE = re.compile('^Interface_Residues_.*(_Repair)?_[0-9]*_AC.fxout$')
INTERACTION_RE = re.compile('^Interaction_.*(_Repair)?_[0-9]*_AC.fxout$')
INDIV_RE = re.compile('^Indiv_energies_.*(_Repair)?_[0-9]*_AC.fxout$')
SUMMARY_RE = re.compile('^Summary_.*(_Repair)?_[0-9]*_AC.fxout$')

@dataclass
class Mutation:
    chain: str
    position: int
    wt: str
    mut: str

    @staticmethod
    def from_foldx_str(foldx_str):
        """Generate a Mutation from a FoldX formated mutation string"""
        foldx_str = foldx_str.strip(';')
        return Mutation(foldx_str[1], int(foldx_str[2:-1]), foldx_str[0], foldx_str[-1])

def combine_interaction(mutations, foldx_paths, wt_path):
    """Combine Interaction files"""
    with open(wt_path, 'r') as wt_file:
        wt_fields = wt_file.readlines()[-1].strip().split('\t')
    wt_floats = [float(i) for i in wt_fields[3:27]]
    wt_ints = [int(i) for i in wt_fields[27:]]

    float_cols = ['intraclashesgroup1', 'intraclashesgroup2', 'interaction_energy',
                  'backbone_hbond', 'sidechain_hbond', 'van_der_waals', 'electrostatics',
                  'solvation_polar', 'solvation_hydrophobic', 'van_der_waals_clashes',
                  'entropy_sidechain', 'entropy_mainchain', 'sloop_entropy',
                  'mloop_entropy', 'cis_bond', 'torsional_clash', 'backbone_clash',
                  'helix_dipole', 'water_bridge', 'disulfide', 'electrostatic_kon',
                  'partial_covalent_bonds', 'energy_ionisation', 'entropy_complex'
    ]
    int_cols = ['number_of_residues', 'interface_residues', 'interface_residues_clashing',
                'interface_residues_vdw_clashing', 'interface_residues_bb_clashing']
    float_diff = [f'diff_{i}' for i in float_cols]
    int_diff = [f'diff_{i}' for i in int_cols]

    print('chain', 'position', 'wt', 'mut',
          'chain1', 'chain2', *[i for p in zip(float_cols, float_diff) for i in p],
          *[i for p in zip(int_cols, int_diff) for i in p], sep='\t')

    for path in foldx_paths:
        num = int(path.split('_')[-2])
        mutation = mutations[num - 1]
        with open(path, 'r') as foldx_file:
            line = foldx_file.readlines()[-1].strip().split('\t')

        chain1 = line[1]
        chain2 = line[2]
        line_floats = [float(i) for i in line[3:27]]
        line_ints = [int(i) for i in line[27:]]
        diff_floats = [mut - wt for mut, wt in zip(line_floats, wt_floats)]
        diff_ints = [mut - wt for mut, wt in zip(line_ints, wt_ints)]

        print(mutation.chain, mutation.position, mutation.wt,
              mutation.mut, chain1, chain2,
              *[f'{v:.8f}' for p in zip(line_floats, diff_floats) for v in p],
              *[v for p in zip(line_ints, diff_ints) for v in p], sep='\t')

def combine_individual_energies(mutations, foldx_paths):
    """Combine Individual Energy files"""
    print('chain', 'position', 'wt', 'mut', 'group', 'total_energy',
          'backbone_hbond', 'sidechain_hbond', 'van_der_waals', 'electrostatics',
          'solvation_polar', 'solvation_hydrophobic', 'van_der_waals_clashes',
          'entropy_sidechain', 'entropy_mainchain', 'sloop_entropy', 'mloop_entropy',
          'cis_bond', 'torsional_clash', 'backbone_clash', 'helix_dipole',
          'water_bridge', 'disulfide', 'electrostatic_kon', 'partial_covalent_bonds', 'energy_ionisation', 'entropy_complex', sep='\t')
    for path in foldx_paths:
        num = int(path.split('_')[-2])
        mutation = mutations[num - 1]
        with open(path, 'r') as foldx_file:
            lines = foldx_file.readlines()[9:]

        for line in lines:
            print(mutation.chain, mutation.position, mutation.wt, mutation.mut,
                  *line.strip().split('\t')[1:], sep='\t')

def combine_interface_residues(mutations, foldx_paths, wt_path):
    """Combine Interface Residue files"""
    with open(wt_path, 'r') as wt_file:
        wt_residues = set(x[1:] for x in wt_file.readlines()[-1].strip().split('\t'))

    print('chain', 'position', 'wt', 'mut', 'residues_lost',
          'residues_gained', 'interface_residues', sep='\t')
    for path in foldx_paths:
        num = int(path.split('_')[-2])
        mutation = mutations[num - 1]
        with open(path, 'r') as foldx_file:
            interface_residues = foldx_file.readlines()[-1].strip().split('\t')
        interface_set = set(x[1:] for x in interface_residues)
        residues_lost = wt_residues - interface_set
        residues_gained = interface_set - wt_residues
        print(mutation.chain, mutation.position, mutation.wt, mutation.mut,
              ','.join(residues_lost), ','.join(residues_gained),
              ','.join(interface_residues), sep='\t')

def combine_summary(mutations, foldx_paths, wt_path):
    """Combine Summary files"""
    with open(wt_path, 'r') as wt_file:
        wt_energy = [float(i) for i in wt_file.readlines()[-1].strip().split('\t')[3:]]

    print('chain', 'position', 'wt', 'mut', 'group1', 'group2',
          'intraclashesgroup1', 'diff_intraclashesgroup1',
          'intraclashesgroup2', 'diff_intraclashesgroup2',
          'interaction_energy', 'diff_interaction_energy',
          'stabilitygroup1', 'diff_stabilitygroup1',
          'stabilitygroup2', 'diff_stabilitygroup2', sep='\t')
    for path in foldx_paths:
        num = int(path.split('_')[-2])
        mutation = mutations[num - 1]
        with open(path, 'r') as foldx_file:
            summary = foldx_file.readlines()[-1].strip().split('\t')
        chain1 = summary[1]
        chain2 = summary[2]
        summary = [float(i) for i in summary[3:]]
        diff = [mut - wt for mut, wt in zip(summary, wt_energy)]

        print(mutation.chain, mutation.position, mutation.wt, mutation.mut,
              chain1, chain2, *[f'{v:.8f}' for p in zip(summary, diff) for v in p],
              sep='\t')

def read_mutations(path):
    """
    Read mutations from a FoldX individual list file
    """
    with open(path, 'r') as individual_list:
        mutations = [Mutation.from_foldx_str(i.strip()) for i in individual_list]
    return mutations

def detect_filetype(file):
    """
    Detect FoldX filetype from the filenames
    """
    file = Path(file).stem
    filetype = file.split('_')[0].lower()
    if filetype not in ('interface', 'interaction', 'indiv', 'summary'):
        raise ValueError('Unknown type detected')
    return filetype

def main(args):
    """Main"""
    filetype = detect_filetype(args.wt)
    mutations = read_mutations(args.mutations)
    root = args.foldx.rstrip("/")
    files = [f for f in os.listdir(root) if f.endswith('.fxout')]
    files = sorted(files, key=lambda x: int(x.split('_')[-2]))

    if filetype == 'interface':
        files = [f'{root}/{f}' for f in files if INTERFACE_RE.match(f)]
        combine_interface_residues(mutations, files, args.wt)

    elif filetype == 'interaction':
        files = [f'{root}/{f}' for f in files if INTERACTION_RE.match(f)]
        combine_interaction(mutations, files, args.wt)

    elif filetype == 'indiv':
        files = [f'{root}/{f}' for f in files if INDIV_RE.match(f)]
        combine_individual_energies(mutations, files)

    elif filetype == 'summary':
        files = [f'{root}/{f}' for f in files if SUMMARY_RE.match(f)]
        combine_summary(mutations, files, args.wt)

def parse_args():
    """Parse Arguments"""
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('mutations', metavar='M', help="Individual list file")
    parser.add_argument('wt', metavar='W', help="WT FoldX AC File")
    parser.add_argument('foldx', metavar='F', help="FoldX output directory")

    return parser.parse_args()

if __name__ == '__main__':
    main(parse_args())
