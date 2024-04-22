#!/usr/bin/env python
# coding: utf-8

"""
Python code to generate motifs for ASSAM [1]

Code developed inside the Virtuous Projct (https://virtuoush2020.com/)

----------------
Acknowledgements
This project has received funding from the European Union’s Horizon 2020 research and innovation programme under the Marie Skłodowska-Curie grant agreement Action (GA No. 872181)
----------
References
[1] Nadzirin, N., E. J. Gardiner, P. Willett, P. J. Artymiuk, and M. Firdaus-Raih. ‘SPRITE and ASSAM: Web Servers for Side Chain 3D-Motif Searching in Protein Structures’. Nucleic Acids Research 40, no. W1 (1 July 2012): W380–86. https://doi.org/10.1093/nar/gks401.
----------------
Version history:
- Version 1.0 - 13/12/2022
----------------
"""

__version__ = '1.0'
__author__ = 'Marco Cannariato, Lampros Androutsos, Eric Zizzi, Lorenzo Pallante'

import MDAnalysis as mda
import argparse
import os
import sys
import warnings

from motif_search import *

warnings.filterwarnings('ignore')


def parse_options():
    parser = argparse.ArgumentParser(description='''Tool to create motifs for ASSAM analysis''')
    parser.add_argument('--pdb',
                        metavar='complex.pdb',
                        dest='pdb',
                        type=str,
                        help='Protein-ligand complex structure in PDB format.',
                        required=True)
    parser.add_argument('--xtc',
                        metavar='complex.xtc',
                        dest='xtc',
                        type=str,
                        help='Protein-ligand complex trajectory in xtc format.',
                        required=True)
    parser.add_argument('--outdir',
                        metavar='od',
                        dest='outdir',
                        type=str,
                        help='Path to the output directory.',
                        required=True)
    parser.add_argument('--chainProtein',
                        metavar='cP',
                        dest='chainProtein',
                        nargs='+',  # needed if you have more chains in the protein to define the pocket
                        type=str,
                        help='Chain of the protein in the input structure.',
                        required=True)
    parser.add_argument('--chainLigand',
                        metavar='cL',
                        dest='chainLigand',
                        type=str,
                        help='Chain of the ligand in the input structure.',
                        required=True)
    parser.add_argument('--kflag',
                        metavar='k',
                        dest='kflag',
                        type=int,
                        help='Clusters for kmeans: 0 for optimization, other for specific k.',
                        default=0)
    parser.add_argument('--dist',
                        metavar='dist',
                        dest='d',
                        type=float,
                        help='Maximum distance in Angstrom from ligand to identify pocket residues.'
                             ' Default is 10 A.',
                        default=10)

    args, unknown = parser.parse_known_args()
    expected = ["--pdb", "--xtc", "--outdir", "--chainP", "--chainl", "--kflag"
        , "--dist"]

    for cmd in unknown:
        if cmd not in expected:
            print('\n\n{0}\nUnknown command found in your command line: "{1}".\n{0}\n\n'.format('*' * 50, cmd))
            sys.exit(1)

    return args


def check_directories(dirs):
    # check if directories exist and create them if not
    for dir in dirs:
        if not os.path.exists(dir):
            os.makedirs(dir)


def check_files(files):
    for file in files:
        if not os.path.isfile(file):
            print(
                "\n\n{0}\nError! Provided file '{1}' does not exist or is not a file.\n{0}\n\n".format('*' * 50, file))
            sys.exit(3)


def main(args):
    '''Main function to run the script.'''

    # 0- Check if a trajectory is provided, if not, skip clustering and generate motif from the input structure
    if args.xtc == 'None':

        # Froma single pdb input, retrive the file for the pocket residues and the protein-ligand complex (name coherent with the rest of the code) 
        u = mda.Universe(args.pdb)
        sel = "byres (chainID " + " ".join(args.chainProtein) + " and around " + str(
            args.d) + " chainID " + args.chainLigand + ")"
        u.select_atoms(sel + " or (chainID " + args.chainLigand + ")").write("{0}/centroid_0.pdb".format(args.outdir))
        u.select_atoms("chainID " + " ".join(args.chainProtein) + " " + args.chainLigand).write(
            "{0}/protein_centroid_0.pdb".format(args.outdir))

        print("\nCreating motifs from PDB input structure...\n")
    else:
        # 1- Clustering
        print("\nClustering protein-ligand trajectory...\n")

        clustering(pdbFile=args.pdb, trjFile=args.xtc, outPath=args.outdir, chainProtein=args.chainProtein,
                   chainLig=args.chainLigand, kClusters=args.kflag
                   , dist=args.d)

        # 2- Centroids
        print("\nCreating motifs from selected centroids...\n")

    # create directory outdir/motifs
    if not os.path.exists(os.path.join(args.outdir, 'motifs')):
        os.makedirs(os.path.join(args.outdir, 'motifs'))
    # create motifs
    create_motifs(centrPath=args.outdir, outPath=os.path.join(args.outdir, 'motifs'),
                  chainLig=args.chainLigand)
    print("\nDONE!\n\nMotif files ready to be used with Assam and stored in {0}/motifs\n\n".format(args.outdir))


def run_script():
    args = parse_options()
    check_directories([args.outdir])
    if args.xtc != 'None':
        check_files([args.pdb, args.xtc])
    else:
        check_files([args.pdb])
    main(args)


if __name__ == '__main__':
    run_script()
