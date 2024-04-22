#!/usr/bin/env python
# coding: utf-8

"""
Python code to analyse ASSAM output. This code was implemented within the Virtuous Projct (https://virtuoush2020.com/)

----------------
Acknowledgements
This project has received funding from the European Union’s Horizon 2020 research and innovation programme under the Marie Skłodowska-Curie grant agreement Action (GA No. 872181)
----------
References
[1] 
----------------
Version history:
- Version 1.0 - 03/01/2023
----------------
"""

__version__ = '1.0'
__author__ = 'Marco Cannariato, Lampros Androutsos, Eric Zizzi, Lorenzo Pallante'

import argparse
import sys
import warnings

from motif_search import *

warnings.filterwarnings('ignore')


def parse_options():
    parser = argparse.ArgumentParser(description='''Tool to analyse the results file of
    ASSAM webserver. ASSAM results should be named hit_[right/left]_motif_[number].txt''')
    parser.add_argument('--inDir',
                        metavar='inDir',
                        dest='inDir',
                        type=str,
                        help='Input directory containing ASSAM result files.',
                        required=True)
    parser.add_argument('--pdb',
                        metavar='complex.pdb',
                        dest='pdb',
                        type=str,
                        help='Protein-ligand complex structure in PDB format.',
                        required=True)
    parser.add_argument('--clDir',
                        metavar='clDir',
                        dest='clDir',
                        type=str,
                        help='Directory containing Protein-ligand clusters.',
                        required=False)
    parser.add_argument('--outDir',
                        metavar='od',
                        dest='outdir',
                        type=str,
                        help='Path to the output directory.',
                        required=True)
    parser.add_argument('--ligChain',
                        metavar='ln',
                        dest='ligChain',
                        type=str,
                        help='Chain Label of the ligand in the input structure.',
                        default=None)
    parser.add_argument('--protChain',
                        metavar='prt',
                        dest='protChain',
                        type=str,
                        help='Chain Label of the protein in the input structure.',
                        default=None)
    parser.add_argument('--cpu',
                        metavar='cpus',
                        dest='cpus',
                        type=int,
                        help='Number of cpus to be used during sasa and docking.'
                             'Default is 1',
                        default=1)
    parser.add_argument('--verbose',
                        dest='verbose',
                        help='Be verbose.',
                        default=False,
                        action='store_true')
    parser.add_argument('--keep',
                        dest='keep',
                        help='Keep intermediate files (SASA, Docking, PDBs). Default not kept.',
                        default=False,
                        action='store_true')
    parser.add_argument('--dbFolder',
                        metavar='db',
                        dest='db',
                        type=str,
                        help='Path to the folder containing the PDB files of the original database.',
                        default=None)
    parser.add_argument('--sasaThreshold',
                        metavar='sasaThreshold',
                        dest='sasaThreshold',
                        type=float,
                        help='SASA threshold. Default is 0.75.',
                        default=0.75)
    parser.add_argument('--dockThreshold',
                        metavar='dockThreshold',
                        dest='dockThreshold',
                        type=float,
                        help='Docking score threshold. Default is 0.1.',
                        default=0.1)
    parser.add_argument('--maxHits',
                        dest='maxHits',
                        type=int,
                        help='Maximum number of hits to be retrieved. Default is None, meaning all possible hits will be retrieved',
                        default=-1)

    args, unknown = parser.parse_known_args()
    expected = ["--inDir", "--clDir", "--outDir", "--ligChain", "--cpu",
                "--verbose", "--pdb", "--keep", "--protChain", "--dbFolder", "--dockThreshold", "--sasaThreshold"]

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
            sys.exit(1)


def main(args):
    # create folder for the links to PDB files in the original DB
    pdb_link_folder = "{0}/pdb".format(args.outdir)

    # create directory if it does not exist
    if not os.path.exists(pdb_link_folder):
        os.mkdir(pdb_link_folder)

    df_res = pd.DataFrame()
    pdbcount = 0

    # cycle over all blocks of the original DB
    for block in sorted(os.listdir(args.db)):

        # folder of the pdb files of the current block
        pdbFolder = "{0}/{1}/pdb".format(args.db, block)
        # count number of PDB in the blocks
        pdbcount += len(os.listdir(pdbFolder))

        # folder of the ASSAM results of the current block
        ASSAM_resFolder = "{0}/{1}/results".format(args.inDir, block)
        df_res = pd.concat(
            [df_res, analyse_assam_results(resPath=ASSAM_resFolder, pdbPath=pdbFolder, verbose=args.verbose)],
            ignore_index=True)

        # create symbolic link to the PDB files retrieved by ASSSAM
        if not df_res.empty:  # check if the dataframe is not empty (no ASSAM results)
            for pdb in df_res['PDBid']:
                # check if the link already exists (one PDB can be retrieved multiple times for each motif)
                if not os.path.exists("{0}/{1}.pdb".format(pdb_link_folder, pdb)):
                    os.symlink("{0}/{1}.pdb".format(pdbFolder, pdb), "{0}/{1}.pdb".format(pdb_link_folder, pdb))

    # if no results from ASSAM, exit
    if df_res.empty:
        print("\n\n{0}\nError! No results from ASSAM.\n{0}\n\n".format('*' * 50))
        sys.exit(1)

    # order results by ascending RMSD
    df_res = df_res.sort_values(by=['RMSD'])

    # Save ASSAM results
    df_res.to_csv(args.outdir + os.sep + 'Results_ASSAM.txt', sep='\t', index=False)

    # drop duplicates -> preserve only one hit per PDB (lowest RMSD value)
    df_res_uniq = df_res.drop_duplicates(subset=['PDBid'])

    # if maxHits is provided, keep only the first maxHits
    if args.maxHits != -1:
        print('Keeping only the first {0} hits from ASSAM'.format(args.maxHits))
        df_res_max = df_res_uniq.iloc[:args.maxHits]
        df_res_max.to_csv(args.outdir + os.sep + 'Results_ASSAM_best{0}Hits.txt'.format(args.maxHits), sep='\t',
                          index=False)
    else:
        df_res_max = df_res_uniq

    #############################
    # Filter 1 - calculating SASA
    sasaPath = "{0}/SASA".format(args.outdir)

    # create directory if it does not exist
    if not os.path.exists(sasaPath):
        os.mkdir(sasaPath)

    print("Computing SASA of hit sites and filtering (SASA threshold = {0})...".format(args.sasaThreshold))
    df_sasa = sasa_filter_joblib(pdbPath=pdb_link_folder, outPath=sasaPath, clPath=args.clDir,
                                 df=df_res_max, sasaThreshold=args.sasaThreshold, cpus=args.cpus)

    df_sasa.to_csv("{0}/Results_sasa.txt".format(args.outdir), sep='\t', index=False)

    ################################
    # Filter 2 - calculating Docking
    dockPath = "{0}/docking".format(args.outdir)

    # create directory if it does not exist
    if not os.path.exists(dockPath):
        os.mkdir(dockPath)

    print(
        "\nPerforming docking on hit sites and filtering (Docking Score threshold = {0})...".format(args.dockThreshold))
    df_dock = docking_filter_plants_joblib(df=df_sasa, pdbPath=pdb_link_folder, outPath=dockPath, clPath=args.clDir,
                                           chainLig=args.ligChain, chainProtein=args.protChain,
                                           dockingThreshold=args.dockThreshold, cpus=args.cpus)

    # order results by ascending docking score
    df_dock = df_dock.sort_values(by=['Docking'])

    # save results
    df_dock.to_csv("{0}/Results.txt".format(args.outdir), sep='\t', index=False)

    # print statistics
    print("\n\n*{0:5d} PDB files in the original database".format(pdbcount))
    print("*{0:5d} hits sites from ASSAM analysis.".format(len(df_res)))
    print("*{0:5d} unique PDBs".format(len(df_res_uniq)))
    if args.maxHits != -1:
        print("*{0:5d} hits sites kept from ASSAM analysis due to MaxHits request.".format(len(df_res_max)))
    print("*{0:5d} sites passed the SASA filter (SASA Threshold = {1}).".format(len(df_sasa), args.sasaThreshold))
    print(
        "*{0:5d} sites passed the docking filter (Docking Threshold = {1}).\n".format(len(df_dock), args.dockThreshold))

    # Removing intermediate files if keep is not set
    if not args.keep:
        # shutil.rmtree(pdbPath)
        shutil.rmtree(sasaPath)
        shutil.rmtree(dockPath)

    print("\nDONE!\n\nCheck {0}/Results.txt for all hits!\n\n".format(args.outdir))


def run_script():
    args = parse_options()
    check_directories([args.outdir, args.inDir, args.clDir])
    check_files([args.pdb])
    main(args)


if __name__ == '__main__':
    run_script()
# python ../code/analyse_assam.py --inDir ../assam_linux/results/ --pdb ../Testing_files/solute.pdb --outDir out --ligName LIG --cpu 8 --clDir prova/
