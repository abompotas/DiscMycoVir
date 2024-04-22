#!/usr/bin/env python
# coding: utf-8

"""
Python code to run the master script of the CBP tool
This code was implemented within the Virtuous Projct (https://virtuoush2020.com/)

----------------
Acknowledgements
This project has received funding from the European Union’s Horizon 2020 research and innovation programme under the Marie Skłodowska-Curie grant agreement Action (GA No. 872181)
----------
References
[1] 
----------------
Version history:
- Version 1.0 - 08/02/2023
----------------
"""

__version__ = '1.0'
__author__ = 'Lorenzo Pallante, Marco Cannariato, Lampros Androutsos, Eric Zizzi'

import argparse
import multiprocessing
import numpy as np
import os
import pandas as pd
import shutil
import subprocess
import sys
# time the code
import time
from art import *

from motif_search import _plot_plip_interactions, _plot_filtering, retrieve_info, _write_js_file, _write_vmd_file, \
    functional_annotation_charts


def code_name():
    # print the name of the code and center it
    # tprint("CBP", font="tarty1", chr_ignore=True)
    tprint("Virtuous   Pocketome", font="tarty1", chr_ignore=True)
    tprint(
        "Code to compare the binding pocket of a target protein-ligand complex with the entire human proteome and identify the most similar binding clefts".upper(),
        font="awesome", chr_ignore=True)

    # print authors and version
    print("\n")
    print("Version: " + __version__)
    print("Authors: " + __author__)
    print("\n")


# parse arguments
def parse_options():
    # variables needed for script - generate_assam_motifs.py
    # DEFAULT VALUES
    kflag = 0  # # number of clusters to be generated; 0 = optimization with silhouette score, 1 = user defined number of clusters
    dist = 10  # Maximum distance from ligand to identify pocket residues (Default 10A)
    cpus = int(multiprocessing.cpu_count() / 2)  # variable for script - run_assam.py
    sasaThreshold = 0.75  # Threshold for the SASA filter step (Default 0.75)
    dockThreshold = 0.1  # Threshold for the docking filter step (Default 0.1)
    outfolder = os.getcwd() + 'VirtuousPocketome'  # Folder where the output will be saved

    parser = argparse.ArgumentParser(description=code_name())
    parser.add_argument('--outfolder', type=str, default=outfolder, help='Folder where the output will be saved')
    parser.add_argument('--complexPDB', type=str, default=None, help='PDB file of the complex')
    parser.add_argument('--complexXTC', type=str, default=None, help='XTC file of the complex')
    parser.add_argument('--chainProtein', type=str, default=None, help='Chain of the protein')
    parser.add_argument('--chainLigand', type=str, default=None, help='Chain of the ligand')
    parser.add_argument('--db2bcompared', type=str, default=None,
                        help='Folder containing the files of the proteins to be compared')
    parser.add_argument('--cpus', type=int, default=cpus, help='Number of cpus to be used')
    parser.add_argument('--kflag', type=int, default=kflag,
                        help='Number of clusters to be generated; 0 = optimization with silhouette score, 1 = user defined number of clusters')
    parser.add_argument('--dist', type=int, default=dist,
                        help='Maximum distance from ligand to identify pocket residues (Default 10A)')
    parser.add_argument('--sasaThreshold', type=float, default=sasaThreshold,
                        help='Threshold for the SASA filter step (Default 0.75)')
    parser.add_argument('--dockThreshold', type=float, default=dockThreshold,
                        help='Threshold for the docking filter step (Default 0.1)')
    parser.add_argument('--extensive', type=bool, default=False,
                        help='If True, the code will run the extensive analysis parsing all the possible ASSAM hits (Default False)')
    # add a config file: if the file is passed, other  arguments should be ignored
    parser.add_argument('--config', type=str, default=None, help='Config file')
    return parser.parse_args()


def run_cbp(args):
    start_time = time.time()

    # define required variables to None
    required_args = {}
    required_args['complex_pdb'] = False
    required_args['chain_protein'] = False
    required_args['chain_ligand'] = False
    required_args['db_2b_compared'] = False

    # define optional arguments
    optional_args = {}
    optional_args['complex_xtc'] = False
    optional_args['cpus'] = False
    optional_args['kflag'] = False
    optional_args['dist'] = False
    optional_args['sasaThreshold'] = False
    optional_args['dockThreshold'] = False
    optional_args['outfolder'] = False
    optional_args['extensive'] = False

    config_flag = False
    if 'config' in args:
        if args.config is not None:
            config_flag = True

    # define the variables
    if config_flag:
        print('Reading config file "{0}"...'.format(args.config))

        # read the config file
        with open(args.config, 'r') as f:
            for line in f.readlines():
                if line.startswith(('#', ' ', '\n')):
                    continue
                elif line.startswith('complex_pdb'):
                    complex_pdb = line.split('=')[1].strip()
                    required_args['complex_pdb'] = True
                elif line.startswith('complex_xtc'):
                    complex_xtc = line.split('=')[1].strip()
                    optional_args['complex_xtc'] = True
                elif line.startswith('chain_protein'):
                    chain_protein = line.split('=')[1].strip()
                    required_args['chain_protein'] = True
                elif line.startswith('chain_ligand'):
                    chain_ligand = line.split('=')[1].strip()
                    required_args['chain_ligand'] = True
                elif line.startswith('outfolder'):
                    outfolder = line.split('=')[1].strip()
                    optional_args['outfolder'] = True
                elif line.startswith('db_2b_compared'):
                    db_2b_compared = line.split('=')[1].strip()
                    required_args['db_2b_compared'] = True
                elif line.startswith('cpus'):
                    cpus = int(line.split('=')[1].strip())
                    optional_args['cpus'] = True
                elif line.startswith('kflag'):
                    kflag = int(line.split('=')[1].strip())
                    optional_args['kflag'] = True
                elif line.startswith('dist'):
                    dist = int(line.split('=')[1].strip())
                    optional_args['dist'] = True
                elif line.startswith('sasaThreshold'):
                    sasaThreshold = float(line.split('=')[1].strip())
                    optional_args['sasaThreshold'] = True
                elif line.startswith('dockThreshold'):
                    dockThreshold = float(line.split('=')[1].strip())
                    optional_args['dockThreshold'] = True
                elif line.startswith('extensive'):
                    extensive = line.split('=')[1].strip()
                    # string to boolean
                    if extensive == 'True' or extensive == 'true':
                        extensive = True
                    elif extensive == 'False' or extensive == 'false':
                        extensive = False
                    else:
                        print('Unknown parameter: {}'.format(line))
                        sys.exit(1)
                    optional_args['extensive'] = True
                else:
                    print('Unknown parameter: {}'.format(line))
                    sys.exit(1)
    else:
        print('Reading passed arguments ...')

        outfolder = args['outfolder']
        complex_pdb = args['complexPDB']
        required_args['complex_pdb'] = True
        chain_protein = args['chainProtein']
        required_args['chain_protein'] = True
        chain_ligand = args['chainLigand']
        required_args['chain_ligand'] = True
        db_2b_compared = args['db2bcompared']
        required_args['db_2b_compared'] = True
        if 'complexXTC' in args:
            complex_xtc = args['complexXTC']
            optional_args['complex_xtc'] = True
        if 'cpus' in args:
            cpus = args['cpus']
            optional_args['cpus'] = True
        if 'kflag' in args:
            kflag = args['kflag']
            optional_args['kflag'] = True
        if 'dist' in args:
            dist = args['dist']
            optional_args['dist'] = True
        if 'sasaThreshold' in args:
            sasaThreshold = args['sasaThreshold']
            optional_args['sasaThreshold'] = True
        if 'dockThreshold' in args:
            dockThreshold = args['dockThreshold']
            optional_args['dockThreshold'] = True
        if 'extensive' in args:
            extensive = args['extensive']
            optional_args['extensive'] = True

    if optional_args['complex_xtc'] == False:
        complex_xtc = "None"

    # check if all arguments are passed
    for key, value in required_args.items():
        if not value:
            print('Error: missing argument {0}'.format(key))
            sys.exit(1)

    # check if the output folder exists
    if not os.path.exists(outfolder):
        os.makedirs(outfolder)

    # define needed directories
    motifs_folder = outfolder + '01-gen_assam_motifs' + os.sep + 'motifs' + os.sep
    cluster_folder = outfolder + '01-gen_assam_motifs' + os.sep

    # Define the steps to be run
    step1 = True  # Motifs generation
    step2 = True  # ASSAM Search
    step3 = True  # Multi-step filtering
    step4 = True  # Analysis

    # 1 - generate_assam_motifs.py
    # run the script to generate the motifs for ASSAM
    print('\nSTEP 1 - MOTIFS SEARCH')
    print('Generating motifs for ASSAM...')

    if step1:

        res = subprocess.run(['python', 'generate_assam_motifs.py',
                              '--pdb', complex_pdb,
                              '--xtc', complex_xtc,
                              '--chainProtein', chain_protein,
                              '--chainLigand', chain_ligand,
                              '--outdir', outfolder + '01-gen_assam_motifs',
                              '--kflag', str(kflag),
                              '--dist', str(dist)])

        if res.returncode != 0:
            print('Error in motif generation')
            return

    # print elapsed time
    print('Elapsed time: {0:.2f} seconds'.format(time.time() - start_time))

    # 2 - Run ASSAM
    #     The user can provide only the pdb folder containing proteins to be compared or both the precalculated vek and master files
    print('\nStep 2 - SIMILARIY SEARCH')
    print('Running ASSAM similarity search algorithm ...')

    if step2:
        res = subprocess.run(['python', 'run_assam.py',
                              '--dbFolder', db_2b_compared,
                              '--motifsFolder', motifs_folder,
                              '--outdir', outfolder + '02-run_assam',
                              '--cpu', str(cpus),
                              '--verbose'])

        if res.returncode != 0:
            print('Error in ASSAM run')
            return

    # print elapsed time
    print('Elapsed time: {0:.2f} seconds'.format(time.time() - start_time))

    # 3 - Analyse ASSAM output
    #     ex. python /home/lorenzo/Documenti/GitHub/CBP/code/analyse_assam.py --inDir 02-run_assam/results/ --pdb ../Testing_files/solute.pdb --clDir 01-gen_assam/ --outDir 03-analyse_assam --ligChain C --cpu 8 --verbose --keep --pdbFolder ../database/20230201/pdb_clean/
    print('\nStep 3 - MULTI-STEP FILTERING')
    print('Analysing ASSAM output and filtering hits according to SASA and Docking ...')

    if step3:

        # if the extensive run is not requested, set maxHits to 1000
        if not extensive:
            maxHits = 100
        else:
            maxHits = -1  # all hits will be considered

        res = subprocess.run(['python', 'analyse_assam.py',
                              '--inDir', outfolder + '02-run_assam',
                              '--pdb', complex_pdb,
                              '--clDir', cluster_folder,
                              '--outDir', outfolder + '03-analyse_assam',
                              '--ligChain', chain_ligand,
                              '--protChain', chain_protein,
                              '--cpu', str(cpus),
                              '--verbose',
                              '--keep',
                              '--dbFolder', db_2b_compared,
                              '--sasaThreshold', str(sasaThreshold),
                              '--dockThreshold', str(dockThreshold),
                              '--maxHits', str(maxHits)])

        if res.returncode != 0:
            print('Error in ASSAM analysis')
            return

    # 4 - ANALYSE RESULTS
    # Create plots, visualisation states and genomics analysis

    print('\nStep 4 - ANALYSIS')
    print('Analysing results ...')

    if step4:

        # create the folder for the analysis
        analysis_folder = outfolder + '04-analysis' + os.sep
        if not os.path.exists(analysis_folder):
            os.makedirs(analysis_folder)
        # create folder for the visualisation
        vis_folder = analysis_folder + 'visualisation' + os.sep
        if not os.path.exists(vis_folder):
            os.makedirs(vis_folder)

        # copy Results.txt and all_interactions.csv to the analysis folder
        shutil.copy(outfolder + '03-analyse_assam' + os.sep + 'Results.txt', analysis_folder + 'Results.txt')
        shutil.copy(outfolder + '01-gen_assam_motifs' + os.sep + 'all_interactions.csv',
                    vis_folder + 'all_interactions.csv')

        # A) PLOT THE INTERACTIONS OF THE MOTIFS
        _plot_plip_interactions(vis_folder + 'all_interactions.csv', analysis_folder + 'interactions.png')

        # B) PLOT THE MULTI-STEP FILTERING PROCESS
        # count original pdbs
        pdbcount = 0
        # cycle over all blocks of the original DB
        for block in sorted(os.listdir(db_2b_compared)):
            # count number of PDB in the blocks
            pdbcount += len(os.listdir("{0}/{1}/pdb".format(db_2b_compared, block)))

        assam_results_file = outfolder + '03-analyse_assam' + os.sep + 'Results_ASSAM.txt'
        sasa_results_file = outfolder + '03-analyse_assam' + os.sep + 'Results_sasa.txt'
        dock_results_file = outfolder + '03-analyse_assam' + os.sep + 'Results.txt'

        files = [assam_results_file, sasa_results_file, dock_results_file]

        # count number of hits in each file
        hits = []
        for file in files:
            with open(file) as f:
                count = sum(1 for _ in f)
            hits.append(count)

        # plot the filtering
        _plot_filtering(tot_pdb=pdbcount, assam_pdb=hits[0], sasa_pdb=hits[1], dock_pdb=hits[2],
                        outputfile=analysis_folder + 'filtering.png')

        # C) GENERATE .vmd and .js files for the visualisation     
        # copy centroid_*.pdb files from the 01-gen_assam_motifs folder to the analysis folder
        for file in os.listdir(outfolder + '01-gen_assam_motifs'):
            if file.startswith('protein_centroid_') and file.endswith('.pdb'):
                shutil.copy(outfolder + '01-gen_assam_motifs' + os.sep + file, vis_folder + file)

        # generate the .js file
        _write_js_file(vis_folder, chain_protein, chain_ligand)
        # generate the .vmd file
        _write_vmd_file(vis_folder, chain_protein, chain_ligand)

        # C) RETRIEVE INFO FROM HITS
        df_list = pd.read_csv(analysis_folder + 'Results.txt', sep='\t')
        ids, infos = retrieve_info(df_list, keypath=os.getcwd())
        np.savetxt(analysis_folder + 'UniProt_ids.txt', ids, fmt='%s')

        # Open the file for reading
        with open(analysis_folder + 'UniProt_ids.txt', 'r') as file:
            # Read the IDs from the file and store them in a list
            ids = file.read().splitlines()

        # Join the IDs using a comma separator to get them in the desired format
        ids_formatted = ','.join(ids)

        # create folder for the enrichment analysis
        enrich_folder = analysis_folder + 'enrichment_analysis' + os.sep
        if not os.path.exists(enrich_folder):
            os.makedirs(enrich_folder)

        functional_annotation_charts(ids_formatted, analysis_folder + 'UniProt_ids.txt', outputfolder=enrich_folder)

        # remove file all_interactions.csv
        os.remove(vis_folder + 'all_interactions.csv')


# run main
if __name__ == '__main__':
    start_time = time.time()

    # clear screen
    os.system('clear')

    # parse arguments
    args = parse_options()

    # run the CBP tool
    run_cbp(args)

    # print elapsed time
    print('\nElapsed time: %.2f s' % (time.time() - start_time))
