# -*- coding: utf-8 -*-

"""
Python code to search for motifs in a protein-ligand complex

Code developed inside the Virtuous Projct (https://virtuoush2020.com/)

----------------
Acknowledgements
This project has received funding from the European Union’s Horizon 2020 research and innovation programme under the Marie Skłodowska-Curie grant agreement Action (GA No. 872181)
----------
References
[1]
----------------
Version history:
- Version 1.0 - 13/12/2022
----------------
"""

__version__ = '1.0'
__author__ = 'Marco Cannariato, Lampros Androutsos, Eric Zizzi, Lorenzo Pallante'

import MDAnalysis as mda
import MDAnalysis.analysis.encore as encore
import concurrent.futures
import fileinput
import glob
from gromacs_py import gmx
import matplotlib.pyplot as plt
import multiprocessing
import numpy as np
import os
import pandas as pd
import re
import requests
import seaborn as sns
import shutil
import subprocess
import subprocess
import sys
from joblib import Parallel, delayed
from matplotlib.colors import ListedColormap
from openmm.app import PDBFile
from pdbfixer import PDBFixer
from plip.structure.preparation import PDBComplex
from sklearn import metrics
from ssl import SSL_ERROR_EOF
from suds import *
from suds.client import Client
from tqdm import tqdm

try:
    from openbabel import pybel
except ImportError:
    print('\n{0}\nOpenbabel is not present in the current environment\n{0}\n'.format('#' * 80))
    sys.exit(1)

# suppress openbabel warnings
pybel.ob.obErrorLog.StopLogging()
pybel.ob.obErrorLog.SetOutputLevel(0)

# suppress warnings
import warnings

warnings.filterwarnings('ignore')

# suppress simtk warnings
import logging

logging.disable(logging.WARNING)


def define_paths():
    """
    Function to define paths to programs
    """
    global MGLPATH, prepLig, prepRec, altLoc, MGL_pythonsh, pdbqt2pdb

    # search for MGLTools path
    for l in os.environ['PATH'].split(':'):
        r = re.search('MGL', l)
        if r:
            flag = True
            splitted = list()
            for ll in l.split('/'):
                if re.search('MGL', ll):
                    splitted.append(ll)
                    flag = False
                elif flag:
                    splitted.append(ll)
                else:
                    break
            MGLPATH = '/'.join(splitted)

    # define paths to MGLTools scripts             
    prepLig = "{0}/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py".format(MGLPATH)
    prepRec = "{0}/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py".format(MGLPATH)
    altLoc = "{0}/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_pdb_split_alt_confs.py".format(MGLPATH)
    MGL_pythonsh = "{0}/bin/pythonsh".format(MGLPATH)
    pdbqt2pdb = "{0}/MGLToolsPckgs/AutoDockTools/Utilities24/pdbqt_to_pdb.py".format(MGLPATH)


def CheckGMXExit(result, logpath="gromacs.log"):
    """
    This function checks the ouput code of a subprocess call to a gromacs
    utility and warns and exits if it is non-zero
    """
    if result.returncode != 0:
        print("FATAL: Error encountered in {0} Check {1} for details!".format(result.args[1], logpath))
        sys.exit(1)


def clustering(pdbFile, trjFile, outPath, chainProtein, chainLig, kClusters=0, dist=10):
    '''
    Function to cluster the trajectory of Protein-Ligand complex using MDAnalysis Kmeans algorithm.

    ---------------------------
    Input:
        pdbFile   ->   Input structure containing protein-ligand complex
        trjFile   ->   Equilibrium trajectory of Protein-Ligand complex
        outPath   ->   Path of output folder
        chainProtein  ->  Chain of the protein
        chainLig   ->   Chain of the ligand
        kClusters ->   Number of clusters/centroids to form (Default: 0-kmeans optimization)
        dist      ->   Maximum distance from ligand to identify pocket residues (Default 10A)

    Output:
        Centroids of clusters with at least 5% of frames, numbered according to gromacs clustering output
    '''

    # Clustering of input trajectory
    outStr = "{0}/centroids.pdb".format(outPath)

    # create mda universe
    u = mda.Universe(pdbFile, trjFile)

    # select protein and ligand using chain IDs and distance
    chainProtein = [i.upper() for i in chainProtein]
    chainLig = chainLig.upper()
    sel = "byres (chainID " + " ".join(chainProtein) + " and around " + str(dist) + " chainID " + chainLig + ")"

    # Check if the selection contains atoms
    if len(u.select_atoms(sel).atoms) == 0:
        print(
            "No atoms selected: check the chain IDs of protein and ligand and their relative distance in Angstrom!\n\n")
        sys.exit(1)

    # Calculating the RMSD matrix
    # rmsfer = RMSD(sel, verbose=True).run()     ------ this computes the rmsd over time, not the matrix

    # Build a universe containing only the selection to be used in clustering
    u_cluster = mda.core.universe.Merge(u.select_atoms(sel).atoms)
    indexes_pocket = "index " + " ".join([str(i) for i in u.select_atoms(sel).atoms.indices])
    coordinates = []
    for t in u.trajectory:
        coordinates.append(u.select_atoms(indexes_pocket).positions)
    coordinates = np.array(coordinates)
    u_cluster.load_new(coordinates, order='fac')
    u_cluster.filename = "{0}/top_for_cluster.pdb".format(outPath)
    u_cluster.atoms.write(u_cluster.filename)

    # Flatten coordinate matrix into n_frame x n_coordinates. Will be used to compute silhouette
    coordinates = np.reshape(coordinates, (coordinates.shape[0], -1))

    if kClusters == 0:
        silhouette = {}
        print("Tuning the number of clusters...")
        for k in tqdm(range(2, 12)):
            # clusters = KMeans(n_clusters=k, random_state=0).fit(rmsfer.results.rmsd)
            clusters = encore.cluster(u_cluster, method=encore.KMeans(n_clusters=k, verbose=0))
            labels = np.zeros(u_cluster.trajectory.n_frames)
            for i in clusters.get_ids():
                labels[clusters.clusters[i].elements] = i
            silhouette[k] = (metrics.silhouette_score(coordinates, labels))

            # print("Silhouette Coefficient: %0.3f" % metrics.silhouette_score(rmsfer.results.rmsd, clusters.labels_))

        # run clustering using the optimal number of clusters and save the centroids
        best_k = max(silhouette, key=silhouette.get)
        print("The optimal number of cluster is {0}".format(int(best_k)))
        clusters = encore.cluster(u_cluster, method=encore.KMeans(n_clusters=best_k, verbose=0))

    else:
        print("Clustering using {0} clusters as requested by user".format(kClusters))
        clusters = encore.cluster(u_cluster, method=encore.KMeans(n_clusters=kClusters, verbose=0))

    i = 0
    for cluster in clusters.get_centroids():
        u.trajectory[cluster]
        u.select_atoms("(" + indexes_pocket + ") or (chainID " + chainLig + ")").write(
            "{0}/centroid_{1}.pdb".format(outPath, i))
        u.select_atoms("chainID " + " ".join(chainProtein) + " " + chainLig).write(
            "{0}/protein_centroid_{1}.pdb".format(outPath, i))
        i += 1


def create_motifs(centrPath, outPath, chainLig):
    '''
    Function to create the motifs, ensemble of residues interacting with ligand, for each cluster.

    ---------------------------
    Input:
        centrPath   ->   Path to the folder containing the centroids
        outPath     ->   Path of output folder
        chainLig    ->   Chain of the ligand

    Output:
        Motifs of provided clusters, numbered according to gromacs clustering output
    '''
    # Dictionaty of interaction labels
    interaction_labels = {'hydroph_interaction': 'HI',
                          'hbond': 'HB', 'pistack': 'PS',
                          'pication': 'PC', 'saltbridge': 'SB',
                          'halogenbond': 'HA', 'waterbridge': 'WB',
                          'metal_complex': 'MC'}
    all_interactions = list()

    # Get the list of centroids in acending order
    centroids = sorted(glob.glob("{0}/centroid_*".format(centrPath)))

    # Analyze each centroid
    k = 0
    for centroid in centroids:
        u = mda.Universe(centroid)
        sel = u.select_atoms('chainID {0}'.format(chainLig))

        # format string for PLIP
        ss = "{0}:{1}:{2}".format(np.unique(sel.resnames)[0], chainLig, np.unique(sel.resnums)[0])

        # Load structure into PLIP
        mol = PDBComplex()
        mol.load_pdb(centroid)
        # Analyze interactions
        mol.analyze()
        # Extract the interactions of the ligand
        interactions = mol.interaction_sets[ss]
        idres = interactions.interacting_res

        # Split the residue ID in chain and number
        chain_res = [res[-1] for res in idres]
        num_res = [int(res[:-1]) for res in idres]

        # order according to the num_res
        chain_res = [x for _, x in sorted(zip(num_res, chain_res))]
        num_res.sort()

        # Create the motif
        site_resgroups = u.atoms.split('residue')
        motif = mda.AtomGroup([], u)
        # cycle over the interacting residues
        for num, chain in zip(num_res, chain_res):
            for site in site_resgroups:
                if (num == np.unique(site.resids)[0]) and (chain == np.unique(site.chainIDs)[0]):
                    motif = motif + site

        # write the motif
        motif.write("{0}/motif_{1}.pdb".format(outPath, k))

        # gather the interaction types
        for ii in interactions.all_itypes:
            res = ii.reschain + str(ii.resnr)
            tp = interaction_labels[type(ii).__name__]
            all_interactions.append(str(k) + '-' + res + '-' + tp)

            # increment the counter of the centroids
        k += 1

    try:
        # remove duplicates from all_interactions
        all_interactions = list(set(all_interactions))

        # summarize plip results in a dataframe
        motif_ids = np.array([x.split('-')[0] for x in all_interactions])
        res = np.array([x.split('-')[1] for x in all_interactions])
        tp = np.array([x.split('-')[2] for x in all_interactions])

        data = pd.DataFrame({'motif': motif_ids,
                             'res': res, 'tp': tp}).pivot_table(index='res', columns='motif', values='tp',
                                                                aggfunc=lambda x: ','.join(x)).reset_index().fillna('')
        data['num'] = [int(x[1:]) for x in data['res']]
        data = data.sort_values(by=['num']).drop(columns=['num'])
        # set res column as index
        data = data.set_index('res')
        data.to_csv("{0}/all_interactions.csv".format(centrPath))
    except:
        # let's not crash if there is any problem in summarizing the plip results
        pass


def analyse_assam_results(resPath, pdbPath, verbose=True):
    '''
    Function to analyse ASSAM results and save the results in a dataframe

    ---------------------------
    Input: 
        resPath     ->   Path to the folder containing the ASSAM results
        pdbPath     ->   Path to the folder containing the PDB files
        verbose     ->   Print the results of the analysis

    Output:
        Dataframe containing the results of the analysis
    '''

    # Get the list of files
    resFolders = sorted(glob.glob("{0}/*".format(resPath)))

    # Get the list of PDB files
    files = list()
    for folder in resFolders:
        for file in os.listdir(folder):
            if file.endswith(".LPL") or file.endswith(".LPS"):
                files.append(os.path.join(folder, file))

    # Create the dataframe  
    df_results = pd.DataFrame()

    # Analyse each file
    for file in files:

        df = pd.DataFrame()
        motifName = file.split(os.sep)[-2]  # save the name of the motif
        if os.path.splitext(file)[-1] == '.LPL':
            superimp = 'left'
        elif os.path.splitext(file)[-1] == '.LPS':
            superimp = 'right'
        ids, rmsds, matches = list(), list(), list()
        flag = False
        for line in fileinput.input(file):
            # print(line)
            if 'KERBm' in line:
                ids.append([line.split()[i] for i in [3, 6]])
            elif 'RMSD' in line:
                rmsds.append(line.split()[2])
            elif 'TRANSFORM' in line:
                flag = True
                foo = list()
            elif flag:
                if ">" in line:
                    splitted = line.split()
                    foo.append(" ".join(splitted[1:10]))
                else:
                    matches.append(foo)
                    flag = False

        if ids == []:
            continue
        if flag and len(matches) == len(ids) - 1:
            matches.append(foo)
        elif not flag and len(matches) != len(ids):
            print("\nSomething went wrong reading {0}!\n".format(file))
            sys.exit(1)

        # check last hit
        if int(ids[-1][-1]) != len(matches[-1]):
            ids[-1][-1] = str(len(matches[-1]))

        query = list()
        hit = list()
        for match in matches:
            temp_q = list()
            temp_m = list()
            for ind in match:
                ind = ind.replace('>', '')
                ind = ind.replace('<', '')
                q, m = ind.split('matches')
                qs = q.split()
                if not (qs[0].isalpha() or qs[0].isdigit()):
                    qs = [qs[0][0], qs[0][1:], qs[1]]
                elif qs[1].isalpha():
                    qs = ['None', qs[0], qs[1]]
                qs[1] = qs[1].replace('+', '')
                qs[1] = qs[1].replace('-', '')
                q = ":".join(qs)
                ms = m.split()
                if not (ms[0].isalpha() or ms[0].isdigit()):
                    ms = [ms[0][0], ms[0][1:], ms[1]]
                elif ms[1].isalpha():
                    ms = ['None', ms[0], ms[1]]
                ms[1] = ms[1].replace('+', '')
                ms[1] = ms[1].replace('-', '')
                m = ":".join(ms)
                temp_q.append(q)
                temp_m.append(m)
            query.append(','.join(temp_q))
            hit.append(','.join(temp_m))

        for n, i in enumerate(ids):
            pdbid = i[0]
            if query[n].split(':')[0] == 'None':
                for line in fileinput.input("{0}/{1}.pdb".format(pdbPath, pdbid)):
                    if 'ATOM' in line:
                        if line[21].isalpha():
                            query[n] = query[n].replace('None', line[21])
                        break

        df['motif'] = [motifName for i in range(len(ids))]
        df['supimp'] = [superimp for i in range(len(ids))]
        df['PDBid'] = [ids[i][0] for i in range(len(ids))]
        df['matchSize'] = [ids[i][-1] for i in range(len(ids))]
        df['query'] = query
        df['hit'] = hit
        df['RMSD'] = rmsds
        df_sorted = df.sort_values(by=['RMSD'])
        df_results = pd.concat([df_results, df_sorted], ignore_index=True)

    # check if assam results is not empty
    if not df_results.empty:
        # Check hits with None as chain and correct them if necessary
        df_s = df_results.sort_values(by=['RMSD'])
        df_final = complete_df(df_s, pdbPath)
    else:
        df_final = df_results

    return df_final


def pdb_parse(pdb_file):
    '''
    Function to parse a pdb file and extract the following information:
    - molid
    - chain
    - organism
    - organism_scientific
    - doi

    Parameters
    ----------
    pdb_file : str
        Path to the pdb file.

    Returns
    -------
    df : pandas dataframe
        Dataframe containing the information extracted from the pdb file.
    '''
    # initialize the pdb dictionary and useful variables
    pdb_dict = {}
    molid = []
    organism = []
    chain = []
    header_chain = list()
    header_organism = list()
    organism_scientific = list()
    molecule = list()

    # read the pdb file
    with open(pdb_file, 'r') as f:

        lines = f.readlines()

        for i, line in enumerate(lines):

            # get molid, chain and organism
            if line.startswith('COMPND'):
                if 'MOL_ID:' in line:
                    header_chain.append(line)
                    molid.append(line.replace(';', ':').split(':')[1].strip())
                elif line.split()[2] == 'CHAIN:':
                    header_chain.append(line)
                    # chain.append(line.replace(';', ':').split(':')[1].strip())
                    temp_line = line.replace(';', ':').split(':')[1].strip()

                    # check if line ends with ; otherwise other chains are stores in the next line
                    k = 0

                    # read the next line
                    while lines[i + k].strip()[-1] != ';':
                        k += 1
                        line_next = lines[i + k]
                        header_chain.append(line_next)
                        temp_line = temp_line + ' ' + line_next.replace(';', ' ').split()[2].strip()

                    chain.append(temp_line)

                elif line.split()[2] == 'MOLECULE:':
                    header_chain.append(line)
                    molecule.append(line.replace(';', ':').split(':')[1].strip())

            # get the organism
            elif line.startswith('SOURCE'):
                if 'MOL_ID:' in line:
                    header_organism.append(line)
                if line.split()[2] == 'ORGANISM_TAXID:':
                    header_organism.append(line)
                if line.split()[2] == 'ORGANISM_SCIENTIFIC:':
                    header_organism.append(line)

            # retrieve doi
            elif line.startswith('JRNL'):
                if line.split()[1] == 'DOI':
                    doi = line.split()[-1].strip()

            # retrieve pdb id
            elif line.startswith('HEADER'):
                pdb_id = line.split()[-1].strip()

    # check the correspondance between molid/chain and molid/organism
    for i in np.argwhere(['MOL_ID' in x for x in header_organism]).ravel():
        try:
            if 'ORGANISM_TAXID' in header_organism[i + 2]:
                organism.append(header_organism[i + 2].replace(';', ':').split(':')[1].strip())
            else:
                # if the organism is not specified, we assume it is human
                organism.append('9606')
        except IndexError:
            organism.append('9606')

        # check the correspondance between molid/chain and molid/organism_scientific
        try:
            if 'ORGANISM_SCIENTIFIC' in header_organism[i + 1]:
                organism_scientific.append(header_organism[i + 1].replace(';', ':').split(':')[1].strip())
            else:
                # if the organism is not specified, we assume it is human
                organism_scientific.append('HOMO SAPIENS')
        except IndexError:
            organism_scientific.append('HOMO SAPIENS')

    # set default values if no information is found
    if 'doi' not in locals():
        doi = 'N/A'

    if 'pdb_id' not in locals():
        pdb_id = 'N/A' * len(molid)

    pdb_dict['pdb_id'] = pdb_id
    pdb_dict['molid'] = molid
    pdb_dict['molecule'] = molecule
    pdb_dict['chain'] = chain
    pdb_dict['organism'] = organism
    pdb_dict['organism_scientific'] = organism_scientific
    pdb_dict['doi'] = doi * len(molid)

    df = pd.DataFrame(pdb_dict)
    return df


def complete_df(df_old, pdbPath):
    '''
    Function to complete the table containing information about ASSAM hits (v2 using function pdb_parse)
    '''

    # create copy of the original dataframe
    df = df_old.copy()

    parsed_doi = []
    parsed_organism_scientific = []
    parsed_organism = []
    parsed_molecule = []

    for i, ids in enumerate(df.PDBid):
        ids = str(ids)

        # find file in the pdbPath with the same name as the pdb id
        file = "{0}/{1}.pdb".format(pdbPath, ids)
        # parse the pdb file for information
        df_parsed = pdb_parse(file)

        # get the chains that are in the hit list in the original datadrame
        hit_chains = list(set([h.split(':')[0] for h in df.iloc[i].hit.split(',')]))

        # select only information related to the hit chains
        df_parsed_hit = df_parsed[df_parsed.chain.str.contains('|'.join(hit_chains))]
        parsed_organism_scientific.append(', '.join(list(set(df_parsed_hit.organism_scientific.values))))
        parsed_organism.append(', '.join(list(set(df_parsed_hit.organism.values))))
        parsed_molecule.append(', '.join(list(set(df_parsed_hit.molecule.values))))
        parsed_doi.append(', '.join(list(set(df_parsed_hit.doi.values))))

    df['class'] = parsed_molecule
    df['organism'] = parsed_organism_scientific
    df['taxid'] = parsed_organism
    df['doi'] = parsed_doi
    return (df)


def _run_sasa(pdbCode, match, pdbPath, outPath):
    '''
    Function to be able to run the SASA filter with Parallel

    '''
    # for pdbCode,match in tqdm(zip(df['PDBid'].values,df['hit'].values), total=len(df['PDBid'].values)):
    pdbFile = "{0}/{1}.pdb".format(pdbPath, pdbCode)
    u = mda.Universe(pdbFile)
    outName = "index_{}".format(pdbCode)
    outFile = "{0}/index_{1}.ndx".format(outPath, pdbCode)
    outLog = "{0}/gromacs_{1}.log".format(outPath, pdbCode)
    outArea = "{0}/area_{1}.xvg".format(outPath, pdbCode)

    # create GROMACS selection
    chains = [i.split(':')[0] for i in match.split(',')]
    resid = [i.split(':')[1] for i in match.split(',')]
    if 'None' in chains:
        gmx_sel = ["keep 1"]  # keep protein
        gmx_sel.extend('r ' + ' '.join(map(str, resid)))
    else:
        gmx_sel = ["keep 1", "del 0"]  # remove all groups
        gmx_sel.extend(["chain {0} & r {1}".format(c, r) for c, r in zip(chains, resid)])
        gmx_sel.extend([" | ".join([str(i) for i in range(len(chains))])])
        gmx_sel.extend(["name {0} Pocket".format(len(chains))])
        gmx_sel.extend(['chain ' + ' '.join(map(str, chains))])
        gmx_sel.extend(["name {0} Protein".format(len(chains) + 1)])

    gmx_sel.append('q')

    # Before computing sasa, keep only coordinates to avoid problems with
    # headers of huge PDB files
    coordFile = "{0}/coord_{1}.pdb".format(outPath, pdbCode)
    glog = open(outLog, 'a')
    with open(coordFile, 'w') as f:
        process = subprocess.run(['pdb_keepcoord', pdbFile],
                                 stdout=f,
                                 stderr=glog)
    glog.close()
    CheckGMXExit(process, outLog)

    # create index file (redirect stdout and stderror to log)
    with open(outLog, 'w') as glog:
        try:
            md_sys = gmx.GmxSys(name='gromacs_{}.log'.format(pdbCode), coor_file=coordFile)
            md_sys.add_ndx('\n'.join(gmx_sel) + '\n', ndx_name=outName, folder_out=outPath)
        except:
            print("FATAL: Error encountered in {0} Check {1} for details!".format('gmx make_ndx', glog.name))
            sys.exit(1)

    # now compute SASA 
    with open(outLog, 'a') as glog:
        process = subprocess.run(['gmx', 'sasa',
                                  '-f', coordFile,
                                  '-s', coordFile,
                                  '-o', outArea,
                                  '-surface', 'Protein',
                                  '-output', 'Pocket',
                                  '-n', outFile],
                                 stdout=glog,
                                 stderr=glog)
        CheckGMXExit(process, outLog)
    os.remove(outLog)
    os.remove(outFile)
    os.remove(coordFile)


def sasa_filter_joblib(pdbPath, outPath, clPath, df, sasaThreshold=0.2, cpus=int(multiprocessing.cpu_count() / 2)):
    '''
    Function to implement the SASA filter to select the best hit in Parallel

    :param pdbPath: Path to the folder containing the PDB files
    :param outPath: Path to the folder where the output files will be saved
    :param clPath: Path to the folder containing the cluster files
    :param df: Pandas dataframe containing the informations about the hits
    :param sasaThreshold: % of the original SASA in the query centroids to preserve hits (Default: 0.2)
    :return: Pandas dataframe containing the informations about the hits after the SASA filter

    '''
    # First, compute SASA
    areas = list()

    # run the sasa filter in parallel
    Parallel(n_jobs=cpus)(delayed(_run_sasa)(pdbCode, match, pdbPath, outPath) for pdbCode, match in
                          tqdm(zip(df['PDBid'].values, df['hit'].values), total=len(df['PDBid'].values)))

    # cycle over the out areas files 
    outAreas = sorted(glob.glob("{0}/area_*.xvg".format(outPath)))

    for outArea in outAreas:

        a = np.loadtxt(outArea, comments=['#', '@'], usecols=(2))
        if a.shape:
            areas.append(a.mean())
        else:
            areas.append(a)

        # remove gromacs output file
        os.remove(outArea)

    # create a new dataframe with the sasa values
    df_new = df.copy()
    df_new['SASA'] = areas

    # Second, average SASA in receptor centroids
    # create dictonary with the sasa values for each cluster
    cl_dict = dict()

    # cycle over the clusters
    centroids = sorted(glob.glob("{0}/protein_centroid_*.pdb".format(clPath)))
    for centroid in centroids:
        dd = dict()
        # number of the centroid and the relative motif
        num = os.path.basename(centroid).split('.')[0].split('_')[-1]
        # create protein index for centroid 
        u = mda.Universe(centroid)
        outFile = "{0}/index.ndx".format(outPath)
        outLog = "{0}/gromacs.log".format(outPath)
        outArea = "{0}/area_cl_{1}.xvg".format(outPath, num)
        w = mda.selections.gromacs.SelectionWriter(outFile, mode='w')
        prot = u.select_atoms("protein")
        w.write(prot.residues.atoms, name="Protein", mode='w', number=0)

        # retrieve all the hit res for the current motif
        hitReslist = df_new[df_new['motif'] == 'motif_{0}'.format(num)]['query'].unique()

        # if the list is empty, no hits for the current motif -> skip to the following motif
        if len(hitReslist) == 0:
            continue

        for i, hitRes in enumerate(hitReslist):
            sel_string = _create_selection(hitRes)
            sel = u.select_atoms(sel_string)
            w.write(sel.residues.atoms, name=hitRes, mode='w', number=i + 1)
        w.close()

        with open(outLog, 'a') as glog:
            process = subprocess.run(['gmx', 'sasa',
                                      '-f', centroid,
                                      '-s', centroid,
                                      '-o', outArea,
                                      '-surface', 'Protein',
                                      '-output', '"' + '"; "'.join(hitReslist) + '"',
                                      '-n', outFile],
                                     stdout=glog,
                                     stderr=glog)
            CheckGMXExit(process, outLog)
        os.remove(outLog)
        os.remove(outFile)

        # save the sasa value in the dictionary for each motif hit
        for i, hitRes in enumerate(hitReslist):
            dd[hitRes] = np.loadtxt(outArea, comments=['#', '@'], usecols=(i + 2))
        cl_dict['motif_{0}'.format(num)] = dd
    # Third, filter according to SASA
    # Preserve only hits with SASA > sasaThreshold * SASA in the centroid
    # the comparison is between the SASA of the hit and the SASA of the relative motif
    # default value of sasaThreshold is 0.2
    df_filtered = df_new[[r['SASA'] > cl_dict[r['motif']][r['query']] * sasaThreshold for _, r in df_new.iterrows()]]

    return df_filtered


def sasa_filter(pdbPath, outPath, clPath, df, sasaThreshold=0.2):
    '''
    Function to implement the SASA filter to select the best hit

    :param pdbPath: Path to the folder containing the PDB files
    :param outPath: Path to the folder where the output files will be saved
    :param clPath: Path to the folder containing the cluster files
    :param df: Pandas dataframe containing the informations about the hits
    :param sasaThreshold: % of the original SASA in the query centroids to preserve hits (Default: 0.2)
    :return: Pandas dataframe containing the informations about the hits after the SASA filter

    '''
    # First, compute SASA
    areas = list()
    df_new = df.copy()
    for pdbCode, match in tqdm(zip(df['PDBid'].values, df['hit'].values), total=len(df['PDBid'].values)):
        pdbFile = "{0}/{1}.pdb".format(pdbPath, pdbCode)
        u = mda.Universe(pdbFile)
        outName = "index_{}".format(pdbCode)
        outFile = "{0}/index.ndx".format(outPath)
        outLog = "{0}/gromacs.log".format(outPath)
        outArea = "{0}/area_{1}.xvg".format(outPath, pdbCode)

        # create GROMACS selection
        chains = [i.split(':')[0] for i in match.split(',')]
        resid = [i.split(':')[1] for i in match.split(',')]
        if 'None' in chains:
            gmx_sel = ["keep 1"]  # keep protein
            gmx_sel.extend('r ' + ' '.join(map(str, resid)))
        else:
            gmx_sel = ["keep 1", "del 0"]  # remove all groups
            gmx_sel.extend(["chain {0} & r {1}".format(c, r) for c, r in zip(chains, resid)])
            gmx_sel.extend([" | ".join([str(i) for i in range(len(chains))])])
            gmx_sel.extend(["name {0} Pocket".format(len(chains))])
            gmx_sel.extend(['chain ' + ' '.join(map(str, chains))])
            gmx_sel.extend(["name {0} Protein".format(len(chains) + 1)])

        gmx_sel.append('q')

        # Before computing sasa, keep only coordinates to avoid problems with
        # headers of huge PDB files
        coordFile = "{0}/coord_{1}.pdb".format(outPath, pdbCode)
        glog = open(outLog, 'a')
        with open(coordFile, 'w') as f:
            process = subprocess.run(['pdb_keepcoord', pdbFile],
                                     stdout=f,
                                     stderr=glog)
        glog.close()
        CheckGMXExit(process, outLog)

        # create index file (redirect stdout and stderror to log)
        with open(outLog, 'w') as glog:
            try:
                md_sys = gmx.GmxSys(name='gromacs_{}.log'.format(pdbCode), coor_file=coordFile)
                md_sys.add_ndx('\n'.join(gmx_sel), ndx_name=outName, folder_out=outPath)
            except:
                print("FATAL: Error encountered in {0} Check {1} for details!".format('gmx make_ndx', glog.name))
                sys.exit(1)

        # now compute SASA 
        with open(outLog, 'a') as glog:
            process = subprocess.run(['gmx', 'sasa',
                                      '-f', coordFile,
                                      '-s', coordFile,
                                      '-o', outArea,
                                      '-surface', 'Protein',
                                      '-output', 'Pocket',
                                      '-n', outFile],
                                     stdout=glog,
                                     stderr=glog)
            CheckGMXExit(process, outLog)
        os.remove(outLog)
        os.remove(outFile)
        os.remove(coordFile)
        a = np.loadtxt(outArea, comments=['#', '@'], usecols=(2))
        if a.shape:
            areas.append(a.mean())
        else:
            areas.append(a)

        # remove gromacs output file
        os.remove(outArea)

    df_new['SASA'] = areas

    # Second, average SASA in receptor centroids
    # create dictonary with the sasa values for each cluster
    cl_dict = dict()

    # cycle over the clusters
    centroids = sorted(glob.glob("{0}/protein_centroid_*.pdb".format(clPath)))
    for centroid in centroids:
        dd = dict()
        # number of the centroid and the relative motif
        num = os.path.basename(centroid).split('.')[0].split('_')[-1]
        # create protein index for centroid 
        u = mda.Universe(centroid)
        outFile = "{0}/index.ndx".format(outPath)
        outLog = "{0}/gromacs.log".format(outPath)
        outArea = "{0}/area_cl_{1}.xvg".format(outPath, num)
        w = mda.selections.gromacs.SelectionWriter(outFile, mode='w')
        prot = u.select_atoms("protein")
        w.write(prot.residues.atoms, name="Protein", mode='w', number=0)

        # retrieve all the hit res for the current motif
        hitReslist = df_new[df_new['motif'] == 'motif_{0}'.format(num)]['query'].unique()

        # if the list is empty, no hits for the current motif -> skip to the following motif
        if len(hitReslist) == 0:
            continue

        for i, hitRes in enumerate(hitReslist):
            sel_string = _create_selection(hitRes)
            sel = u.select_atoms(sel_string)
            w.write(sel.residues.atoms, name=hitRes, mode='w', number=i + 1)
        w.close()

        with open(outLog, 'a') as glog:
            process = subprocess.run(['gmx', 'sasa',
                                      '-f', centroid,
                                      '-s', centroid,
                                      '-o', outArea,
                                      '-surface', 'Protein',
                                      '-output', '"' + '"; "'.join(hitReslist) + '"',
                                      '-n', outFile],
                                     stdout=glog,
                                     stderr=glog)
            CheckGMXExit(process, outLog)
        os.remove(outLog)
        os.remove(outFile)

        # save the sasa value in the dictionary for each motif hit
        for i, hitRes in enumerate(hitReslist):
            dd[hitRes] = np.loadtxt(outArea, comments=['#', '@'], usecols=(i + 2))
        cl_dict['motif_{0}'.format(num)] = dd
    # Third, filter according to SASA
    # Preserve only hits with SASA > sasaThreshold * SASA in the centroid
    # the comparison is between the SASA of the hit and the SASA of the relative motif
    # default value of sasaThreshold is 0.2
    df_filtered = df_new[[r['SASA'] > cl_dict[r['motif']][r['query']] * sasaThreshold for _, r in df_new.iterrows()]]

    return df_filtered


def _prepare_pdb(pdbFile, outPath, hitRes, outfile):
    '''
    Function to prepare the PDB file for the analysis

    :param pdbFile: Path to the PDB file
    :param outPath: Path to the folder where the output files will be saved
    :param hitRes: List of the residues in the hits
    :param outfile: Name of the output file
    :return: None

    '''

    localPath = os.getcwd()
    outLog = os.path.abspath("{0}/pdbTools.log".format(outPath))

    # extract one frame
    with open(outLog, 'a') as glog:

        # change directory
        os.chdir(outPath)

        # keep only protein atoms
        f = open('nohetatm.pdb', 'w')
        process = subprocess.run(['pdb_delhetatm', pdbFile],
                                 stdout=f,
                                 stderr=glog)
        f.close()
        CheckGMXExit(process, outLog)

        # keep only protein atoms
        f = open('onlyprotein.pdb', 'w')
        process = subprocess.run(['pdb_delresname',
                                  '-DA,DC,DG,DT,DI,A,C,G,U,I',
                                  'nohetatm.pdb'],
                                 stdout=f,
                                 stderr=glog)
        f.close()
        CheckGMXExit(process, outLog)
        os.remove('nohetatm.pdb')

        # renumber only atoms
        with open("renumbered.pdb", 'w') as f:
            process = subprocess.run(['pdb_reatom', 'onlyprotein.pdb'],
                                     stdout=f,
                                     stderr=glog)
        CheckGMXExit(process, outLog)
        os.remove("onlyprotein.pdb")

        # Removes all non-coordinate records from the file.
        # Keeps only MODEL, ENDMDL, END, ATOM, HETATM, CONECT
        with open("coordinates.pdb", 'w') as f:
            process = subprocess.run(['pdb_keepcoord', 'renumbered.pdb'],
                                     stdout=f,
                                     stderr=glog)
        CheckGMXExit(process, outLog)
        os.remove("renumbered.pdb")

        # Keep only chains of interest
        chains = [hit.split(':')[0] for hit in hitRes.split(',')]
        if 'None' in chains:
            os.rename("coordinates.pdb", "receptor.pdb")
        else:
            selChain = ','.join(set(chains))
            with open("receptor.pdb", 'w') as f:
                process = subprocess.run(['pdb_selchain',
                                          '-{0}'.format(selChain),
                                          'coordinates.pdb'],
                                         stdout=f,
                                         stderr=glog)
                CheckGMXExit(process, outLog)
                os.remove("coordinates.pdb")
    os.chdir(localPath)
    os.rename("{0}/receptor.pdb".format(outPath), outfile)  # output of the function
    os.remove(outLog)


def _fix_pdb(pdbfile, outFile, verbose=False):
    '''
    Function to fix a PDB file
    
    :param pdbfile: Path to the PDB file
    :param outFile: Path to the output file
    :param verbose: Verbose output
    :return: absolute path to the fixed PDB file
    
    '''

    # retrieve absolute path of the pdb file
    pdbpath = os.path.abspath(pdbfile)

    # 0 - Check input file:
    if pdbfile is None:
        print("FATAL: No PDB file provided!")
        sys.exit(1)
    else:
        if verbose:
            print("INFO: PDB file is %s" % pdbpath)
        # Check if file exists:
        if not os.path.isfile(pdbpath):
            print("FATAL: PDB file %s does not exist!" % pdbpath)
            sys.exit(1)

    # Load PDB file:
    if verbose:
        print("INFO: Loading PDB file...")

    # 1 - Remove all water molecules:
    if verbose:
        print("INFO: Removing all water molecules...")
        # Read PDB file:
    with open(pdbpath, 'r') as f:
        lines = f.readlines()
    # Remove all lines containing HOH:
    lines = [l for l in lines if 'HOH' not in l]
    # Write new PDB file:
    with open(outFile, 'w') as f:
        f.writelines(lines)

    # 2 - Run Fixer:
    fixer = PDBFixer(filename=pdbpath)
    if verbose:
        print("INFO: Fixing PDB file...")
        print("INFO: Checking missing residues...")
    fixer.findMissingResidues()
    fixer.findNonstandardResidues()
    if verbose:
        print("INFO: Checking and replacing non-standard residues...")
    fixer.replaceNonstandardResidues()
    fixer.findMissingAtoms()
    if verbose:
        print("INFO: Adding missing atoms...")
    fixer.addMissingAtoms()
    if verbose:
        print("INFO: Adding missing hydrogens at pH = 7.4 ...")
    fixer.addMissingHydrogens(7.4)
    if verbose:
        print("INFO: PDB file fixed!")
    # Save PDB file:
    outFile_path = os.path.abspath(outFile)
    if verbose:
        print("INFO: Saving preprocessed PDB file to {0}...".format(outFile_path))
    PDBFile.writeFile(fixer.topology, fixer.positions, open(outFile, 'w+'), keepIds=True)

    if verbose:
        print("INFO: PDB file saved!")
        print("INFO: Done!")


def FindExec(name):
    '''
    Function to find the path to an executable within the directiories of $PATH, considering also partial matches.
    Warns and exits in case the executable is not found.
    :param name: Name of the executable (e.g. "SPORES")
    '''
    find_flag = False
    for l in np.unique(os.environ['PATH'].split(':')):
        if os.path.isdir(l):
            for g in glob.glob(f'{l}/*'):
                if re.search(name, os.path.basename(g)):
                    this_path = g
                    find_flag = True
                    break
        elif re.search(name, l):
            this_path = l
            find_flag = True
            break

    if not find_flag:
        print(f'Path to {name} executable not found.\n\n\tAborting.')
        sys.exit(1)

    return this_path


def _createConf(outDir, conFile, Protein='protein.mol2', Ligand='ligand.mol2', center=0, radius=10):
    """
    Function to create the configuration file to run docking with PLANTS
    :param outDir: name of the output directory. If False, the code will create a directory in the same path as 'Protein' and with the name dock_'Ligands'.
    
    :param conFile: output filename of configuration file
    :param Protein: name of the prepared protein mol2 file
    :param Ligand: name of the prepared ligand to dock
    :param center: center of the binding site. 
    :param radius: radius of binding site
    :return: output directory where the docking results will be inserted.
    """
    with open(conFile, 'w+') as f:
        f.write(f'''
# scoring function and search settings
search_speed speed1

# input
protein_file {Protein}
ligand_file {Ligand}

# write single mol2 files (e.g. for RMSD calculation)
write_multi_mol2 0
write_protein_conformations 0

# binding site definition
bindingsite_center {center[0]} {center[1]} {center[2]}
bindingsite_radius {radius}
# cluster algorithm
cluster_structures 1
cluster_rmsd 2.0

outside_binding_site_penalty 100.0
scoring_function chemplp
aco_ants 20
aco_evap 0.25
output_dir {outDir}''')
    return outDir


def docking_filter_plants(df, pdbPath, outPath, clPath, chainLig, chainProtein, dockingThreshold=0.1,
                          cpus=int(multiprocessing.cpu_count() / 2)):
    '''
    Function to filter protein hits using the docking affinity score

    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame containing the protein hits
    pdbPath : str
        Path to the pdb file of the proteins
    outPath : str
        Path to the output directory
    clPath : str
        Path to the directory containing the centroids of the original protein-ligand complex
    chainLig : str
        Chain of the ligand (this must be unique)
    chainProtein : str
        Chain of the protein (this must be unique)
    dockingThreshold : float, optional
        Threshold for the docking affinity score. The default is 0.1 (10% of the original docking score to keep the hits during filtering)
    cpus : int, optional
        Number of cpus to use. The default is half of the available cpus.

    Returns
    -------
    pandas.DataFrame
        DataFrame containing the filtered protein hits

    '''

    #  ====== DOCKING ON CENTROIDS OF THE ORIGINAL PROTEIN-LIGAND COMPLEX ======

    # retrieve centroid of the original protein-ligand complex
    centroids = glob.glob("{0}/protein_centroid_*.pdb".format(clPath))

    # list to store the docking scores of the original protein-ligand complex
    original_scores = []

    # cycle over the centroids
    for centroid in centroids:

        u = mda.Universe(centroid)

        # ==== LIGAND PREPARATION ====
        ligSel = u.select_atoms('chainID {0}'.format(chainLig))
        ligPDBFile = "{0}/ligand.pdb".format(outPath)
        ligSel.write(ligPDBFile)

        # Prepare ligand for docking with PLANTS (conversion from pdb to mol2)
        ligMol2 = "{0}/lig_to_prepare.mol2".format(outPath)
        ligFile = "{0}/to_dock.mol2".format(outPath)
        Mol2File = pybel.Outputfile("mol2", ligMol2)
        mol = next(pybel.readfile("pdb", ligPDBFile))
        # change the name of molecule in the mol2 file
        mol.OBMol.SetTitle('LIG')
        Mol2File.write(mol)

        # clear pdb file
        os.remove(ligPDBFile)

        # Preparation of the ligand with SPORES maintaining the protonation from PDB
        spores_exec = FindExec("SPORES")
        outLog = "{0}/prepLig.log".format(outPath)
        with open(outLog, 'a+') as glog:
            process = subprocess.run([spores_exec,
                                      '--mode', 'settypes',
                                      ligMol2,
                                      ligFile],
                                     stdout=glog,
                                     stderr=glog)
            CheckGMXExit(process, outLog)

        os.remove(outLog)
        os.remove(ligMol2)

        # ==== PROTEIN PREPARATION ====

        # Extract the ligand from the original protein-ligand complex pdb file
        ligSel = u.select_atoms('chainID {0}'.format(chainLig))

        # # Extract the protein from the original protein-ligand complex pdb file
        protSel = u.select_atoms('chainID {0}'.format(chainProtein))
        protPDBFile = "{0}/receptor.pdb".format(outPath)
        protSel.write(protPDBFile)

        # convert original receptor to mol2 using SPORES
        protMol2File = "{0}/receptor.mol2".format(outPath)
        curdir = os.getcwd()
        os.chdir(outPath)
        with open(outLog, 'a+') as glog:
            process = subprocess.run([spores_exec,
                                      '--mode', 'splitpdb',
                                      protPDBFile],
                                     stdout=glog,
                                     stderr=glog)
            CheckGMXExit(process, outLog)
        os.remove(outLog)

        # keep only protein 
        for g in glob.glob(os.getcwd() + '/*.mol2'):
            if re.search('protein', g):
                os.rename(g, protMol2File)
            elif re.search('ligand', g):
                os.remove(g)
            elif re.search('water', g):
                os.remove(g)
            elif re.search('metal', g):
                os.remove(g)
        os.chdir(curdir)

        # ==== RESCORING BETWEEN LIGAND AND ORIGINAL PROTEIN ====

        # define site for docking on original receptor
        center = ligSel.atoms.center_of_geometry()
        radius = ligSel.atoms.radius_of_gyration() + 5

        # check if radius is nan
        if np.isnan(radius):
            sys.exit(
                "Radius for docking is nan: please check that the ligand is correctly defined in the pdb file (element in the last column and chainID)")

        # create configuration file
        conFile = "{0}/plantsconfig".format(outPath)
        dockPath = _createConf(outPath + "/run", conFile, Protein=protMol2File, Ligand=ligFile,
                               center=center, radius=radius)

        # add line to config file to only rescore the pose and not to perform docking
        with open(conFile, 'a+') as g:
            g.write('\nrescore_mode no_simplex')

        # Perform rescoring with PLANTS
        plants_exec = FindExec("PLANTS")
        outLog = "{0}/dock.log".format(outPath)
        with open(outLog, 'w+') as glog:
            process = subprocess.run([plants_exec,
                                      '--mode', 'rescore',
                                      conFile],
                                     stdout=glog,
                                     stderr=glog)
            CheckGMXExit(process, outLog)
        os.remove(outLog)

        # Get docking score and clean other data
        original_scores.append(pd.read_csv("{0}/bestranking.csv".format(dockPath))['TOTAL_SCORE'].values[0])

        # clean up
        os.remove(protMol2File)
        os.remove(protPDBFile)
        shutil.rmtree(dockPath)
        os.remove(conFile)

    # save original scores to a txt file
    with open("{0}/original_scores.txt".format(outPath), 'w+') as f:
        for score in original_scores:
            f.write("{0}\n".format(score))

    # compute mean score of the centroids of the original protein-ligand complex
    # this value will be used to filter the protein hits
    original_score = np.mean(original_scores)

    # ================================================================
    # ==== PROTEIN HITS ====

    # Prepare hit receptor for docking
    df_new = df.copy()

    with concurrent.futures.ThreadPoolExecutor(max_workers=cpus) as executor:

        futures = []

        for pdbCode, hitRes in tqdm(zip(df['PDBid'].values, df['hit'].values), total=len(df)):
            pdbFile = os.path.abspath("{0}/{1}.pdb".format(pdbPath, pdbCode))
            recFile = os.path.abspath("{0}/{1}.mol2".format(outPath, pdbCode))

            # prepare receptor pdb
            inFile = outPath + "/{0}_receptor.pdb".format(pdbCode)
            _prepare_pdb(pdbFile, outPath, hitRes, outfile=inFile)
            # fix pdb (missing atoms, protonation, etc.)
            inFile_fixed = outPath + "/{0}_receptor_fixed.pdb".format(pdbCode)
            _fix_pdb(inFile, outFile=inFile_fixed)
            outLog = "{0}/prepRec.log".format(outPath)

            # convert pdb to mol2 using SPORES
            curdir = os.getcwd()
            os.chdir(outPath)

            with open(outLog, 'a+') as glog:
                process = subprocess.run([spores_exec,
                                          '--mode', 'settypes',
                                          inFile_fixed,
                                          recFile],
                                         stdout=glog,
                                         stderr=glog)
            CheckGMXExit(process, outLog)
            os.remove(outLog)

            # remove bad receptor from SPORES if present
            if os.path.exists('{0}/{1}_receptor_fixed_bad.mol2'.format(outPath, pdbCode)):
                os.remove('{0}/{1}_receptor_fixed_bad.mol2'.format(outPath, pdbCode))

            os.chdir(curdir)
            # Definition of the site using mda syntax
            site_sel = _create_selection(hitRes, format='mda')

            # Define the site using mda
            u = mda.Universe(inFile)
            site = u.select_atoms(site_sel)
            center = site.atoms.center_of_geometry()
            radius = site.atoms.radius_of_gyration() + 5

            # create configuration file
            conFile = "{0}/plantsconfig_{1}".format(outPath, pdbCode)
            dockPath = _createConf(outPath + "/run_{0}".format(pdbCode), conFile, Protein=recFile, Ligand=ligFile,
                                   center=center, radius=radius)

            # Perform docking with PLANTS
            plants_exec = FindExec("PLANTS")
            outLog = "{0}/dock_{1}.log".format(outPath, pdbCode)
            # with open(outLog, 'w+') as glog:
            #     process = subprocess.run([plants_exec,
            #                                 '--mode','screen',
            #                                 conFile],
            #                                 stdout=glog,
            #                                 stderr=glog)
            #     CheckGMXExit(process,outLog)
            # os.remove(outLog)
            future = executor.submit(subprocess.run, "{0} --mode screen {1}".format(plants_exec, conFile), shell=True,
                                     stdout=open(outLog, 'w+'), stderr=open(outLog, 'w+'))
            futures.append(future)

        # wait for all jobs to be completed
        concurrent.futures.wait(futures)

        # check all the instances have been completed with no errors
        for pdb, future in zip(df['PDBid'].values, futures):
            if future.exception() is not None:
                raise future.exception()
            elif future.result().returncode != 0:
                raise Exception("PLANTS failed for {0}".format(pdb))
            elif not future.done():
                raise Exception("PLANTS not completed for {0}".format(pdb))

    # ================================================================

    # retrieve docking scores from all hits
    scores = list()

    for pdbCode, hitRes in zip(df['PDBid'].values, df['hit'].values):

        dockPath = "{0}/run_{1}".format(outPath, pdbCode)
        recFile = os.path.abspath("{0}/{1}.mol2".format(outPath, pdbCode))
        conFile = "{0}/plantsconfig_{1}".format(outPath, pdbCode)

        # Get docking score and clean other data
        scores.append(pd.read_csv("{0}/bestranking.csv".format(dockPath))['TOTAL_SCORE'].values[0])

        # convert the ligand from mol2 to pdb
        ligMol2 = glob.glob("{0}/LIG_*.mol2".format(dockPath))[0]
        dOut = "{0}/docked.pdb".format(outPath)
        Mol2File = pybel.Outputfile("pdb", dOut)
        mol = next(pybel.readfile("mol2", ligMol2))
        Mol2File.write(mol)

        # concatenate receptor and ligand pdb filesm to create a complex pdb file
        outfile = "{0}/{1}_complex.pdb".format(outPath, pdbCode)
        inFile = outPath + "/{0}_receptor.pdb".format(pdbCode)
        inFile_fixed = outPath + "/{0}_receptor_fixed.pdb".format(pdbCode)
        process = subprocess.run(['cat', inFile_fixed, dOut], stdout=open(outfile, 'w+'))

        # removing lines not starting with ATOM
        with open(outfile, 'r') as f:
            lines = f.readlines()
        with open(outfile, 'w') as f:
            for line in lines:
                if line.startswith('ATOM'):
                    f.write(line)

        outLog = "{0}/dock_{1}.log".format(outPath, pdbCode)

        # clean up: only preserving the complex pdb file
        os.remove(inFile)
        shutil.rmtree(dockPath)
        os.remove(recFile)
        os.remove(dOut)
        os.remove(conFile)
        os.remove(outLog)
        os.remove(inFile_fixed)

    os.remove(ligFile)

    # Filter the dataframe
    df_new['Docking'] = scores
    t = original_score + dockingThreshold * np.abs(
        original_score)  # get only hits with a docking score lower than the threshold
    df_filtered = df_new[df_new['Docking'].values.astype(float) < t]

    return df_filtered


def _run_plants_docking(pdbCode, hitRes, ligFile, pdbPath, outPath):
    # create folder for the output
    outPath = os.path.abspath(outPath + "/{0}".format(pdbCode))
    os.makedirs(outPath, exist_ok=True)

    # define paths
    pdbFile = os.path.abspath("{0}/{1}.pdb".format(pdbPath, pdbCode))
    recFile = os.path.abspath("{0}/{1}.mol2".format(outPath, pdbCode))
    spores_exec = FindExec("SPORES")
    plants_exec = FindExec("PLANTS")

    # prepare receptor pdb
    inFile = outPath + "/{0}_receptor.pdb".format(pdbCode)
    _prepare_pdb(pdbFile, outPath, hitRes, outfile=inFile)
    # fix pdb (missing atoms, protonation, etc.)
    inFile_fixed = outPath + "/{0}_receptor_fixed.pdb".format(pdbCode)
    _fix_pdb(inFile, outFile=inFile_fixed)
    outLog = "{0}/prepRec.log".format(outPath)

    # convert pdb to mol2 using SPORES
    curdir = os.getcwd()
    os.chdir(outPath)

    with open(outLog, 'a+') as glog:
        process = subprocess.run([spores_exec,
                                  '--mode', 'settypes',
                                  inFile_fixed,
                                  recFile],
                                 stdout=glog,
                                 stderr=glog)
    CheckGMXExit(process, outLog)
    os.remove(outLog)

    # remove bad receptor from SPORES if present
    if os.path.exists('{0}/{1}_receptor_fixed_bad.mol2'.format(outPath, pdbCode)):
        os.remove('{0}/{1}_receptor_fixed_bad.mol2'.format(outPath, pdbCode))

    os.chdir(curdir)
    # Definition of the site using mda syntax
    site_sel = _create_selection(hitRes, format='mda')

    # Define the site using mda
    u = mda.Universe(inFile)
    site = u.select_atoms(site_sel)
    center = site.atoms.center_of_geometry()
    radius = site.atoms.radius_of_gyration() + 5

    # create configuration file
    conFile = "{0}/plantsconfig_{1}".format(outPath, pdbCode)
    dockPath = _createConf(outPath + "/run_{0}".format(pdbCode), conFile, Protein=recFile, Ligand=ligFile,
                           center=center, radius=radius)

    # Perform docking with PLANTS
    plants_exec = FindExec("PLANTS")
    outLog = "{0}/dock_{1}.log".format(outPath, pdbCode)
    with open(outLog, 'w+') as glog:
        process = subprocess.run([plants_exec,
                                  '--mode', 'screen',
                                  conFile],
                                 stdout=glog,
                                 stderr=glog)
        CheckGMXExit(process, outLog)


def docking_filter_plants_joblib(df, pdbPath, outPath, clPath, chainLig, chainProtein, dockingThreshold=0.1,
                                 cpus=int(multiprocessing.cpu_count() / 2)):
    '''
    Function to filter protein hits using the docking affinity score

    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame containing the protein hits
    pdbPath : str
        Path to the pdb file of the proteins
    outPath : str
        Path to the output directory
    clPath : str
        Path to the directory containing the centroids of the original protein-ligand complex
    chainLig : str
        Chain of the ligand (this must be unique)
    chainProtein : str
        Chain of the protein (this must be unique)
    dockingThreshold : float, optional
        Threshold for the docking affinity score. The default is 0.1 (10% of the original docking score to keep the hits during filtering)
    cpus : int, optional
        Number of cpus to use. The default is half of the available cpus.

    Returns
    -------
    pandas.DataFrame
        DataFrame containing the filtered protein hits

    '''

    #  ====== DOCKING ON CENTROIDS OF THE ORIGINAL PROTEIN-LIGAND COMPLEX ======

    # retrieve centroid of the original protein-ligand complex
    centroids = glob.glob("{0}/protein_centroid_*.pdb".format(clPath))

    # list to store the docking scores of the original protein-ligand complex
    original_scores = []

    # cycle over the centroids
    for centroid in centroids:

        u = mda.Universe(centroid)

        # ==== LIGAND PREPARATION ====
        ligSel = u.select_atoms('chainID {0}'.format(chainLig))
        ligPDBFile = "{0}/ligand.pdb".format(outPath)
        ligSel.write(ligPDBFile)

        # Prepare ligand for docking with PLANTS (conversion from pdb to mol2)
        ligMol2 = "{0}/lig_to_prepare.mol2".format(outPath)
        ligFile = "{0}/to_dock.mol2".format(outPath)
        Mol2File = pybel.Outputfile("mol2", ligMol2)
        mol = next(pybel.readfile("pdb", ligPDBFile))
        # change the name of molecule in the mol2 file
        mol.OBMol.SetTitle('LIG')
        Mol2File.write(mol)

        # clear pdb file
        os.remove(ligPDBFile)

        # Preparation of the ligand with SPORES maintaining the protonation from PDB
        spores_exec = FindExec("SPORES")
        outLog = "{0}/prepLig.log".format(outPath)
        with open(outLog, 'a+') as glog:
            process = subprocess.run([spores_exec,
                                      '--mode', 'settypes',
                                      ligMol2,
                                      ligFile],
                                     stdout=glog,
                                     stderr=glog)
            CheckGMXExit(process, outLog)

        os.remove(outLog)
        os.remove(ligMol2)

        # ==== PROTEIN PREPARATION ====

        # Extract the ligand from the original protein-ligand complex pdb file
        ligSel = u.select_atoms('chainID {0}'.format(chainLig))

        # # Extract the protein from the original protein-ligand complex pdb file
        protSel = u.select_atoms('chainID {0}'.format(chainProtein))
        protPDBFile = "{0}/receptor.pdb".format(outPath)
        protSel.write(protPDBFile)

        # convert original receptor to mol2 using SPORES
        protMol2File = "{0}/receptor.mol2".format(outPath)
        curdir = os.getcwd()
        os.chdir(outPath)
        with open(outLog, 'a+') as glog:
            process = subprocess.run([spores_exec,
                                      '--mode', 'splitpdb',
                                      protPDBFile],
                                     stdout=glog,
                                     stderr=glog)
            CheckGMXExit(process, outLog)
        os.remove(outLog)

        # keep only protein 
        for g in glob.glob(os.getcwd() + '/*.mol2'):
            if re.search('protein', g):
                os.rename(g, protMol2File)
            elif re.search('ligand', g):
                os.remove(g)
            elif re.search('water', g):
                os.remove(g)
            elif re.search('metal', g):
                os.remove(g)
        os.chdir(curdir)

        # ==== RESCORING BETWEEN LIGAND AND ORIGINAL PROTEIN ====

        # define site for docking on original receptor
        center = ligSel.atoms.center_of_geometry()
        radius = ligSel.atoms.radius_of_gyration() + 5

        # check if radius is nan
        if np.isnan(radius):
            sys.exit(
                "Radius for docking is nan: please check that the ligand is correctly defined in the pdb file (element in the last column and chainID)")

        # create configuration file
        conFile = "{0}/plantsconfig".format(outPath)
        dockPath = _createConf(outPath + "/run", conFile, Protein=protMol2File, Ligand=ligFile,
                               center=center, radius=radius)

        # add line to config file to only rescore the pose and not to perform docking
        with open(conFile, 'a+') as g:
            g.write('\nrescore_mode no_simplex')

        # Perform rescoring with PLANTS
        plants_exec = FindExec("PLANTS")
        outLog = "{0}/dock.log".format(outPath)
        with open(outLog, 'w+') as glog:
            process = subprocess.run([plants_exec,
                                      '--mode', 'rescore',
                                      conFile],
                                     stdout=glog,
                                     stderr=glog)
            CheckGMXExit(process, outLog)
        os.remove(outLog)

        # Get docking score and clean other data
        original_scores.append(pd.read_csv("{0}/bestranking.csv".format(dockPath))['TOTAL_SCORE'].values[0])

        # clean up
        os.remove(protMol2File)
        os.remove(protPDBFile)
        shutil.rmtree(dockPath)
        os.remove(conFile)

    # save original scores to a txt file
    with open("{0}/original_scores.txt".format(outPath), 'w+') as f:
        for score in original_scores:
            f.write("{0}\n".format(score))

    # compute mean score of the centroids of the original protein-ligand complex
    # this value will be used to filter the protein hits
    original_score = np.mean(original_scores)

    # ================================================================
    # ==== PROTEIN HITS ====

    # Run docking on all the hits in Parallel
    Parallel(n_jobs=cpus)(
        delayed(_run_plants_docking)(pdbCode, hitRes, ligFile, pdbPath, outPath) for pdbCode, hitRes in
        tqdm(zip(df['PDBid'].values, df['hit'].values), total=len(df['PDBid'].values)))

    # ================================================================

    # retrieve docking scores from all hits
    scores = list()

    for pdbCode in df['PDBid'].values:

        resPath = "{0}/{1}".format(outPath, pdbCode)
        dockPath = "{0}/run_{1}".format(resPath, pdbCode)
        conFile = "{0}/plantsconfig_{1}".format(resPath, pdbCode)

        # Get docking score and clean other data
        scores.append(pd.read_csv("{0}/bestranking.csv".format(dockPath))['TOTAL_SCORE'].values[0])

        # convert the ligand from mol2 to pdb
        ligMol2 = glob.glob("{0}/LIG_*.mol2".format(dockPath))[0]
        dOut = "{0}/docked.pdb".format(resPath)
        Mol2File = pybel.Outputfile("pdb", dOut)
        mol = next(pybel.readfile("mol2", ligMol2))
        Mol2File.write(mol)

        # concatenate receptor and ligand pdb filesm to create a complex pdb file
        outfile = "{0}/{1}_complex.pdb".format(outPath, pdbCode)
        inFile_fixed = resPath + "/{0}_receptor_fixed.pdb".format(pdbCode)
        process = subprocess.run(['cat', inFile_fixed, dOut], stdout=open(outfile, 'w+'))

        # removing lines not starting with ATOM
        with open(outfile, 'r') as f:
            lines = f.readlines()
        with open(outfile, 'w') as f:
            for line in lines:
                if line.startswith('ATOM'):
                    f.write(line)

        outLog = "{0}/dock_{1}.log".format(resPath, pdbCode)

        # clean up: only preserving the complex pdb file
        shutil.rmtree(resPath)

    os.remove(ligFile)

    # Filter the dataframe
    df_new = df.copy()
    df_new['Docking'] = scores
    t = original_score + dockingThreshold * np.abs(
        original_score)  # get only hits with a docking score lower than the threshold
    df_filtered = df_new[df_new['Docking'].values.astype(float) < t]

    return df_filtered


def docking_filter(df, pdbPath, outPath, complFile, chainLig, cpus=None):
    '''
    Function to filter protein hits using the docking affinity score (VINA)

    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame containing the protein hits
    pdbPath : str
        Path to the pdb file of the protein
    outPath : str
        Path to the output directory
    complFile : str
        Path to the file containing the protein-ligand coordinates
    chainLig : str
        Chain of the ligand (this must be unique)
    cpus : int, optional
        Number of cpus to use, by default half of the available cpus

    Returns
    -------
    pandas.DataFrame
        DataFrame containing the filtered protein hits

    '''
    # if cpus not defined , use half of the available cpus
    if not cpus:
        cpus = int(multiprocessing.cpu_count() / 2)

    # Define MGL tools paths
    define_paths()

    # Extract the ligand from the original protein-ligand complex pdb file
    u = mda.Universe(complFile)
    ligSel = u.select_atoms('chainID {0}'.format(chainLig))
    ligPDBFile = "{0}/ligand.pdb".format(outPath)
    ligSel.write(ligPDBFile)

    # Prepare ligand for docking
    outLog = "{0}/prepLig.log".format(outPath)
    ligFile = "{0}/ligand.pdbqt".format(outPath)
    ligPDBFile = "{0}/ligand.pdb".format(outPath)

    with open(outLog, 'a+') as glog:
        process = subprocess.run([MGL_pythonsh, prepLig,
                                  '-l', ligPDBFile,
                                  '-A', "bonds_hydrogens",
                                  '-o', ligFile],
                                 stdout=glog,
                                 stderr=glog)
        CheckGMXExit(process, outLog)
    os.remove(outLog)

    # Prepare receptor for docking
    scores = list()
    df_new = df.copy()
    for pdbCode, hitRes in tqdm(zip(df['PDBid'].values, df['hit'].values), total=len(df)):
        pdbFile = os.path.abspath("{0}/{1}.pdb".format(pdbPath, pdbCode))
        recFile = os.path.abspath("{0}/{1}.pdbqt".format(outPath, pdbCode))

        # prepare receptor pdb
        inFile = "{0}/{1}_receptor.pdb".format(outPath, pdbCode)
        _prepare_pdb(pdbFile, outPath, hitRes, outfile=inFile)
        inFile_fixed = False
        outLog = "{0}/prepRec.log".format(outPath)

        # convert pdb to pdbqt using autodock tools          
        with open(outLog, 'a+') as glog:
            process = subprocess.run([MGL_pythonsh, prepRec,
                                      '-r', inFile,
                                      '-A', "bonds_hydrogens",
                                      '-o', recFile],
                                     stdout=glog,
                                     stderr=glog)

        # if MGL Tools fails, try to fix the pdb file with the fixer
        if process.returncode != 0:
            with open(outLog, 'a+') as glog:
                glog.write(
                    "\nWARNING: MGL Tools failed to prepare the receptor, trying to fix the PDB file {0}.pdb...\n".format(
                        pdbCode))

            inFile_fixed = _fix_pdb(inFile, "{0}/receptor_fixed.pdb".format(outPath, pdbCode))

            with open(outLog, 'a+') as glog:
                process = subprocess.run([MGL_pythonsh, prepRec,
                                          '-r', inFile_fixed,
                                          '-o', recFile],
                                         stdout=glog,
                                         stderr=glog)
            CheckGMXExit(process, outLog)

        # Definition of the site using mda syntax
        site_sel = _create_selection(hitRes, format='mda')

        # Define the site using mda
        u = mda.Universe(pdbFile)
        site = u.select_atoms(site_sel)

        # retrieve maximum and minimum x, y, z coordinates
        max_x = site.atoms.positions[:, 0].max()
        min_x = site.atoms.positions[:, 0].min()
        max_y = site.atoms.positions[:, 1].max()
        min_y = site.atoms.positions[:, 1].min()
        max_z = site.atoms.positions[:, 2].max()
        min_z = site.atoms.positions[:, 2].min()

        # Define the center and the side (plus a margin to include all residues) of the search box (triclinic box)
        center = np.array([min_x + (max_x - min_x) / 2, min_y + (max_y - min_y) / 2, min_z + (max_z - min_z) / 2])
        margin = 5
        side = np.array([max_x - min_x, max_y - min_y, max_z - min_z]) + margin

        # Perform docking with vina
        vOut = os.path.abspath("{0}/{1}_vina.pdbqt".format(outPath, pdbCode))
        vLog = os.path.abspath("{0}/{1}_vina.log".format(outPath, pdbCode))
        with open(outLog, 'a+') as glog:
            process = subprocess.run(['vina',
                                      '--receptor', recFile,
                                      '--ligand', ligFile,
                                      '--out', vOut,
                                      '--log', vLog,
                                      '--exhaustiveness', '64',
                                      '--center_x', str(center[0]),
                                      '--center_y', str(center[1]),
                                      '--center_z', str(center[2]),
                                      '--size_x', str(side[0]),
                                      '--size_y', str(side[1]),
                                      '--size_z', str(side[2]),
                                      '--cpu', str(cpus),
                                      '--num_modes', '1'],
                                     stdout=glog,
                                     stderr=glog)
            CheckGMXExit(process, outLog)

        # Retrieve the score from the log file
        line = open(vLog, 'r').readlines()[-2]
        scores.append(float(line.split()[1]))

        # convert the ligand pdbqt file to pdb
        process = subprocess.run([MGL_pythonsh, pdbqt2pdb, '-f', vOut, '-o', vOut.replace('.pdbqt', '.pdb')],
                                 stdout=open(outLog, 'a+'), stderr=open(outLog, 'a+'))
        # concatenate receptor and ligand pdb filesm to create a complex pdb file
        if inFile_fixed:
            process = subprocess.run(['cat', inFile_fixed, vOut.replace('.pdbqt', '.pdb')],
                                     stdout=open('{0}/{1}_complex.pdb'.format(outPath, pdbCode), 'w+'))
        else:
            process = subprocess.run(['cat', inFile, vOut.replace('.pdbqt', '.pdb')],
                                     stdout=open('{0}/{1}_complex.pdb'.format(outPath, pdbCode), 'w+'))

        # clean up
        os.remove(inFile)
        if inFile_fixed:
            os.remove(inFile_fixed)
        os.remove(recFile)
        os.remove(vOut)
        os.remove(vLog)
        os.remove(vOut.replace('.pdbqt', '.pdb'))

    os.remove(ligFile)
    os.remove(ligPDBFile)

    # Filter the dataframe
    df_new['Docking'] = scores
    df_filtered = df_new[df_new['Docking'].values.astype(float) < 0]
    return df_filtered


def _create_selection(hitRes, format='mda'):
    '''
    Convert from hitRes format to selection format for MDAnalysis or Vina
    
    with format='mda' -> chain1:num1:res1,chain2:num2:res2, ... to chainID chain1 and resid num1 or chianID chain2 and resid num2 ...
    with format='vina' -> chain1:num1:res1,chain2:num2:res2, ... to chain1:res1num1,res2num2 ; chain2:...
    
    ;param hitRes: string with the residues to select
    ;param format: format of the selection (Default: 'mda', Allowed: 'mda', 'vina')
    ;return: string with the selection
    '''
    hits = hitRes.split(',')
    selection = {}
    chains, resnums, resnames = list(), list(), list()
    ifchain = True  # if there is no chain, the selection will be without chain

    # Initialization of the dictionary and retriving elements
    for hit in hits:
        chain, resnum, resname = hit.split(':')
        chains.append(chain)
        resnums.append(resnum)
        resnames.append(resname)
        if chain == 'None':
            ifchain = False  # if there is no chain, the selection will be without chain
        if chain in selection.keys():
            continue
        else:
            selection[chain] = list()

    # Format for VINA
    if format == 'vina':

        # Ordering of elements and complete dictionary
        zipl = zip(resnums, resnames, chains)
        sortl = sorted(zipl)
        list1, list2, list3 = [list(tuple) for tuple in zip(*sortl)]
        for resnum, resname, chain in zip(list1, list2, list3):
            selection[chain].append("{0}{1}".format(resname, resnum))

        # create selection string
        if ifchain:
            string = list()
            for chain in selection.keys():
                prov = chain + ":" + ",".join(selection[chain])
                string.append(prov)
            string = " ".join(string)
        else:
            for chain in selection.keys():
                string = ",".join(selection[chain])

    # Format for MDAnalysis
    elif format == 'mda':

        # Ordering of elements and complete dictionary
        # not considering the resname since not necessary for mda selection
        zipl = zip(resnums, chains)
        sortl = sorted(zipl)
        list1, list2 = [list(tuple) for tuple in zip(*sortl)]
        for resnum, chain in zip(list1, list2):
            selection[chain].append("{0}".format(resnum))

        # create selection string
        if ifchain:
            string = list()
            for chain in selection.keys():
                # chainID 'X' and resnum res 1 5 70
                string.append("(chainID {0} and resid ".format(chain) + " ".join(selection[chain]) + ")")
            string = " or ".join(string)
        else:
            for chain in selection.keys():
                string = ("resid " + " ".join(selection[chain]))
    else:
        print("Error! format not allowed!\n")
        sys.exit(55)

    return string


def _get_recursively(search_dict, field, target):
    """
    Takes a dict with nested lists and dicts,
    and searches all dicts for a key of the field
    provided.
    """
    fields_found = []

    for key, value in search_dict.items():

        if (key == field) & (value == target):
            fields_found.append(value)

        elif isinstance(value, dict):
            results = _get_recursively(value, field, target)
            for result in results:
                fields_found.append(result)

        elif isinstance(value, list):
            for item in value:
                if isinstance(item, dict):
                    more_results = _get_recursively(item, field, target)
                    for another_result in more_results:
                        fields_found.append(another_result)

    return fields_found


def retrieve_info(df_list, keypath=False):
    '''
    Function to retrieve information about the proteins identified by the procedure. In detail, this function will get information about:
    - Molecular processes 
    - Biological processes
    - Subcellular locations

    :param df_list: list of dataframes containing the output of the research
    :return: dictionary containing the available information about identified proteins.
    '''
    # initialization of the function
    SERVER_URL = "https://www.ebi.ac.uk/pdbe/api"
    UNIPROT = "/mappings/uniprot"

    # get keywords to search
    if not keypath:
        keypath = '.'
    bioProcess = np.loadtxt(f'{keypath}/info/BioProcess.txt', dtype=str, delimiter=',')
    subLocation = np.loadtxt(f'{keypath}/info/SubLocation.txt', dtype=str, delimiter=',')
    molProcess = np.loadtxt(f'{keypath}/info/MolProcess.txt', dtype=str, delimiter=',')

    # concatenate the dataframes if mutiple are given
    if isinstance(df_list, list) & len(df_list) > 1:
        df = pd.concat(df_list, ignore_index=True)
    elif isinstance(df_list, list) & len(df_list) == 1:
        df = df_list[0].copy()
    else:
        df = df_list.copy()

    # initialize the list of UniProt IDs
    ids = list()
    print(f'Number of PDBs to be searched: {df.shape[0]}')
    print('Searching corresponding UniProt IDs...')
    for _, r in tqdm(df.iterrows(), total=df.shape[0]):
        # get pdb id
        pdbCode = r['PDBid']
        # get chains of interest in the pdb 
        chain2search = [i.split(':')[0] for i in r['hit'].split(',')]

        # get the possible UniProt ID matching the pdb
        full_url = "%s/%s/%s?pretty=false" % (SERVER_URL, UNIPROT, pdbCode.lower())
        json_results = requests.get(full_url).json()

        # select the UniProt ID containing the chains of interest
        for l in list(json_results[pdbCode.lower()]['UniProt'].keys()):
            if 'HUMAN' not in json_results[pdbCode.lower()]['UniProt'][l]['identifier']:
                continue
            for chain in chain2search:
                flag = len(_get_recursively(json_results[pdbCode.lower()]['UniProt'][l], 'chain_id', chain)) > 0
                if flag:
                    ids.append(l)

    # delete duplicate IDs
    ids = np.unique(np.array(ids))
    print(f'Number of unique UniProt IDs: {len(ids)}')
    # retrieve information about the ids
    biop = []
    molp = []
    sloc = []
    for id in ids:
        # get the results from UniProt
        url = f'https://rest.uniprot.org/uniprotkb/stream?query=(accession:{id})&format=json'
        res = requests.get(url).json()['results'][0]

        try:
            for i in res['keywords']:
                try:
                    if (i['category'] == 'Biological process') & (i['name'] in bioProcess):
                        biop.append(i['name'])
                    if (i['category'] == 'Molecular function') & (i['name'] in molProcess):
                        molp.append(i['name'])
                except KeyError as e:
                    print(str(e))
        except KeyError as e:
            print(str(e))
        try:
            for i in res['comments']:
                if i['commentType'] == 'SUBCELLULAR LOCATION':
                    if 'subcellularLocations' in i.keys():
                        for l in i['subcellularLocations']:
                            if l['location']['value'] in subLocation:
                                sloc.append(l['location']['value'])
        except KeyError as e:
            print(str(e))

    # summarize results
    infos = {"BiologicalProcess": biop, "MolecularProcess": molp,
             "SubcellularLocation": sloc}
    return ids, infos


def get_input_info(file):
    """
    Function to get information about the input pdb file. It recognizes the number of chains and:
        - if there are more than two chains, it will display the chains detected, specifying which are ligand and which are protein and
            let the user choose the protein and ligand chains.
        - if there are two chains, it will search for a ligand and a protein chain. If two protein chains are found, it will consider 
            the largest one as the protein and the other as the ligand (peptide).
    :param file: file containing the pdb of interest
    :return: protein chain, ligand chain
    """
    # check if the file is a pdb file
    if not file.endswith('.pdb'):
        print('The input file is not a pdb file!')
        sys.exit(1)

    # load the pdb file
    u = mda.Universe(file)

    # get the total number of chains
    n_chains = len(u.segments)
    # if chains are 2, find ligand and protein chains
    if n_chains < 2:
        print('The number of chains is {}, which is less than 2!'.format(n_chains))
        sys.exit(1)
    elif n_chains == 2:
        # check if a ligand is found
        if len(u.select_atoms('not protein')) > 0:
            ligChain = u.select_atoms('not protein').chainIDs[0]
            protChain = u.select_atoms('protein').chainIDs[0]
            print('Ligand chain: ', ligChain)
            print('Protein chain: ', protChain)
        else:
            # count the number of residues in each chain
            n_residues = [len(i.residues) for i in u.segments]
            chains = [i.atoms.chainIDs[0] for i in u.segments]

            # find the chain with the most residues
            protChain = chains[np.argmax(n_residues)]
            # find the chain with the least residues
            ligChain = chains[np.argmin(n_residues)]
            print('Considering the smallest protein as the ligand chain...')
            print('Ligand chain: ', ligChain)
            print('Protein chain: ', protChain)
    else:
        print('The number of chains is {}, not 2!'.format(n_chains))

        # display the chains and if they are protein or ligand
        for i in u.segments:
            tp = 'protein' if len(i.residues) > 1 else 'ligand'
            print('Chain {} - {}'.format(i.atoms.chainIDs[0], tp))

        # let the user choose the chain with the ligand
        ligChain = input('Please enter the chain with the ligand: ')

        # ask the user if the user wants select also the protein chain
        protChain = input('Do you want to select the protein chain? (y/n): ')
        if protChain == 'y':
            protChain = input('Please enter the chain with the protein: ')
        else:
            print('The protein chain will be selected automatically...')
            # find the chain(s) with the protein, which is (are) the ones with a distance from the ligand smaller than 5A
            protChain = np.unique(u.select_atoms('protein and around 10.0 (chainID {})'.format(ligChain)).chainIDs)
            print('Protein chain: ', protChain)

    return protChain, ligChain


def _plot_plip_interactions(interactions_file, output_file):
    color_dict = {"HI": "forestgreen", "SB": "mediumorchid", "HB": "gold", "PI": "dodgerblue", "": "white",
                  "PC": "tomato"}

    data = pd.read_csv(f'{interactions_file}', index_col=0).fillna('')

    # check if all residues are of the same chain, if so remove the chain information from the residue name
    if len(np.unique([i[0] for i in data.index])) == 1:
        data.index = [i[1:] for i in data.index]

    # here we want to invert index and columns
    data = data.T

    # plot the interactions
    unique_chars, matrix = np.unique(data.values, return_inverse=True)
    fig, ax = plt.subplots(figsize=(8 / 2.54, 4 / 2.54))
    sns.heatmap(matrix.reshape(data.shape), annot=data.values, annot_kws={'fontsize': 8, 'weight': 'bold'}, fmt='',
                linecolor='white', lw=1, clip_on=False,
                cmap=ListedColormap([color_dict[char] for char in unique_chars]), ax=ax,
                xticklabels=data.columns, yticklabels=data.index, cbar=False)

    ax.set_ylabel('Motif')
    ax.set_xlabel('Residue')
    ax.spines['top'].set_visible(True);
    ax.spines['bottom'].set_visible(True)
    ax.spines['right'].set_visible(True);
    ax.spines['left'].set_visible(True)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45)
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0)

    fig.tight_layout()
    fig.savefig(output_file, dpi=300)


def _plot_filtering(tot_pdb, assam_pdb, sasa_pdb, dock_pdb, outputfile):
    # set the palette that will be used (avoid color blindness)
    palette = ["cornflowerblue", "tomato", "mediumorchid"]

    # plot the number of hits and pdb during filtering
    fig, ax = plt.subplots(figsize=(8 / 2.54, 5 / 2.54))
    ax.bar(['Total', 'ASSAM', 'SASA', 'Docking'], [tot_pdb, assam_pdb, sasa_pdb, dock_pdb], color=palette[0], width=0.6)
    ax.set_ylabel('Number of Structures')
    ax.set_yscale('log')

    # insert the number of hits as text
    ax.text(0, tot_pdb + tot_pdb / 4, str(tot_pdb) + "\nhits", color='black', fontweight='bold', ha='center')
    ax.text(1, assam_pdb + assam_pdb / 4, str(assam_pdb) + "\nhits", color='black', fontweight='bold', ha='center')
    ax.text(2, sasa_pdb + sasa_pdb / 4, str(sasa_pdb) + "\nhits", color='black', fontweight='bold', ha='center')
    ax.text(3, dock_pdb + dock_pdb / 4, str(dock_pdb) + "\nhits", color='black', fontweight='bold', ha='center')
    fig.tight_layout()
    sns.despine()

    # save the figure
    fig.savefig(outputfile, dpi=300)


def _write_vmd_file(centrPath, chainProtein, chainLig):
    '''
    Write a vmd script to visualize the protein-ligand complex highlighting the motif residues.
    '''

    # Dictionaty of interaction colors as VMD color ID
    interaction_colors = {'HI': '7',  # hydrophic - green
                          'HB': '3',  # hydrogen bond - orange
                          'PS': '1',  # pi stacking - red
                          'PC': '9',  # oi cation - pink
                          'SB': '11',  # saltbridge - purple
                          'HA': '4',  # halogen bond - yellow
                          'WB': '10',  # water bridge - cyan
                          'MC': '14'  # metal complex - ochra
                          }

    data = pd.read_csv("{0}/all_interactions.csv".format(centrPath)).fillna(' ')

    # cycle over the centroids (columns) and write the vmd script to visualize the protein-ligand complex
    for c in data.columns[1:]:
        pdbfile = "{0}/protein_centroid_{1}.pdb".format(centrPath, c)
        vmdFile = "{0}/centroid_{1}.vmd".format(centrPath, c)

        # initialize the vmd script
        with open(vmdFile, 'w') as f:
            f.write("display projection Orthographic\ndisplay depthcue off\n")
            f.write("display shadows off\ncolor Display {Background} white\n")
            f.write("axes location off\n")
            f.write("mol new {} type pdb first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all\n".format(pdbfile))
            f.write("mol delrep 0 top\nmol representation NewCartoon\n")

            # display the protein chain in silver
            sel = "chain " + " ".join(list(chainProtein))
            f.write("mol color ColorID 6\nmol selection { " + sel + " }\nmol material AOChalky\n")
            f.write("mol addrep top\nmol selupdate 0 top 0\nmol colupdate 0 top 0\n")
            f.write("mol scaleminmax top 0 0 0\nmol smoothrep top 0 0\n")
            f.write("mol drawframes top 0 {now}\nmol showrep top 0 1\n\n\n")

            # display the ligand chain in blue2 color of vmd (23)
            sel = "chain " + chainLig + " and not hydrogen"
            f.write("mol representation Licorice 0.5 12.0 12.0\n")
            f.write("mol color ColorID 23\nmol selection { " + sel + " }\nmol material AOChalky\n")
            f.write("mol addrep top\nmol selupdate 0 top 0\nmol colupdate 0 top 0\n")
            f.write("mol scaleminmax top 0 0 0\nmol smoothrep top 0 0\n")
            f.write("mol drawframes top 0 {now}\nmol showrep top 0 1\n\n\n")

            # display the motif residues coloured by interaction type
            for i in data.index:
                if data.loc[i, c] != ' ':
                    resChain = data.loc[i, 'res'][0]
                    resNum = data.loc[i, 'res'][1:]

                    sel = "chain " + resChain + " and resid " + resNum + " and not hydrogen"
                    f.write("mol representation Licorice 0.5 12.0 12.0\n")
                    f.write("mol color ColorID " + interaction_colors[data.loc[i, c]] + "\n")
                    f.write("mol selection { " + sel + " }\nmol material AOChalky\n")
                    f.write("mol addrep top\nmol selupdate 0 top 0\nmol colupdate 0 top 0\n")
                    f.write("mol scaleminmax top 0 0 0\nmol smoothrep top 0 0\n")
                    f.write("mol drawframes top 0 {now}\nmol showrep top 0 1\n\n\n")


def _write_js_file(centrPath, chainProtein, chainLig):
    '''
    Write a js script to visualize the protein-ligand complex highlighting the motif residues using NGLView.
    '''

    # Import the data with the interaction types
    data = pd.read_csv("{0}/all_interactions.csv".format(centrPath)).fillna(' ')

    # create dictionary of interaction colors for NGLView
    interaction_colors = {'HI': '#03ff00',  # hydrophic - green
                          'HB': '#ff9103',  # hydrogen bond - orange
                          'PS': '#ff0000',  # pi stacking - red
                          'PC': '#ffbab9',  # oi cation - pink
                          'SB': '#b100b1',  # saltbridge - purple
                          'HA': '#ffff00',  # halogen bond - yellow
                          'WB': '#42c8c8',  # water bridge - cyan
                          'MC': '#9e5f01'  # metal complex - ochra
                          }

    interaction_types = interaction_colors.keys()

    # cycle over the centroids (columns) and write the vmd script to visualize the protein-ligand complex
    for ii, c in enumerate(data.columns[1:]):

        pdbfile = "{0}/protein_centroid_{1}.pdb".format(centrPath, c)
        jsFile = "{0}/centroid_{1}.js".format(centrPath, c)

        interactions_strings = []

        for i, inter in enumerate(interaction_types):
            resNum = [str(x)[1:] for x in data.res[data.iloc[:, ii + 1] == inter].values]
            resChain = [str(x)[0] for x in data.res[data.iloc[:, ii + 1] == inter].values]
            # string corresponding to the different types of inrteractions
            interactions_strings.append(' or '.join(['({0} and :{1})'.format(i, j) for i, j in zip(resNum, resChain)]))
            # create a string with all the interactions for the present centroid
            interactions_strings_tot = ' or '.join([x for x in interactions_strings if x != ''])

        # Define the JavaScript code as a multiline string using triple quotes
        js_code = f'''
        // Define the color scheme for the interaction types
        document.addEventListener("DOMContentLoaded", function () {{
            var motifs_color = NGL.ColormakerRegistry.addSelectionScheme(
                [
                    ["#03ff00", "{interactions_strings[0]}"],
                    ["#ff9103", "{interactions_strings[1]}"],
                    ["#ff0000", "{interactions_strings[2]}"],
                    ["#ffbab9", "{interactions_strings[3]}"],
                    ["#b100b1", "{interactions_strings[4]}"],
                    ["#ffff00", "{interactions_strings[5]}"],
                    ["#42c8c8", "{interactions_strings[6]}"],
                    ["#9e5f01", "{interactions_strings[7]}"],
                ],
                "Motifs Color Scheme"
            );

            // Create NGL Stage object
            var stage = new NGL.Stage("viewport", {{ backgroundColor: "#ffffff" }});

            // Handle window resizing
            window.addEventListener("resize", function (event) {{
                stage.handleResize();
            }}, false);

            // Load the PDB file and create representations
            stage.loadFile("{pdbfile}").then(function (o) {{
                // Receptor as cartoon representation
                o.addRepresentation("cartoon", {{
                    sele: ":{chainProtein}",
                    color: "#c1c1c1",
                    material: "flat",
                }});

                // Ligand as licorice representation
                o.addRepresentation("licorice", {{
                    sele: "{chainLig}",
                    color: "#4193d5",
                    radius: 0.4,
                }});

                // Motif residues as licorice representation
                o.addRepresentation("licorice", {{
                    sele: "{interactions_strings_tot}",
                    color: motifs_color,
                    radius: 0.3,
                }});

                // Label motif residues with residue IDs
                o.addRepresentation("label", {{
                    sele: "{interactions_strings_tot} and .CA",
                    color: motifs_color,
                    labelType: "res",
                    xOffset: 1,
                    yOffset: 1,
                    zOffset: 1,
                    radius: 4,
                }});

                // Automatically adjust camera position
                stage.autoView();
            }});

            // Disable mouse scrolling
            window.addEventListener("scroll", function (event) {{
                stage.mouseControls.remove("scroll");
            }}, false);
        }});
        '''

        # Write the JavaScript code to a file
        with open(jsFile, 'w') as file:
            file.write(js_code)


def functional_annotation_charts(uniprotIds, uniprotIds_file, outputfolder):
    url = 'https://david.ncifcrf.gov/webservice/services/DAVIDWebService?wsdl'
    client = Client(url)
    client.wsdl.services[0].setlocation(
        'https://david.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap11Endpoint/')

    # authenticate user email
    client.service.authenticate('l.androutsos@insybio.com')

    # add a list

    inputIds = '{}'.format(uniprotIds)
    idType = 'UNIPROT_ACCESSION'
    listName = 'make_up'
    listType = 0
    client.service.addList(inputIds, idType, listName, listType)

    # setCategories
    categories = ['GOTERM_BP_DIRECT', 'GOTERM_CC_DIRECT', 'GOTERM_MF_DIRECT', 'pathway']
    for category in categories:
        if category == 'pathway':
            categoryString = str(client.service.setCategories('KEGG_PATHWAY,REACTOME_PATHWAY'))
        else:
            categoryString = str(client.service.setCategories(category))

        # getChartReport
        thd = 0.1
        ct = 2
        chartReport = client.service.getChartReport(thd, ct)

        # parse and print chartReport
        if category == 'pathway':
            resF = outputfolder + '{}.chartReport.txt'.format('KEGG_PATHWAY_REACTOME_PATHWAY')
        else:
            resF = outputfolder + '{}.chartReport.txt'.format(categoryString)

        with open(resF, 'w') as fOut:
            fOut.write(
                'Category\tTerm\tCount\t%\tPvalue\tGenes\tList Total\tPop Hits\tPop Total\tGeneRatio\tBonferroni\tp.adjust\tFDR\n')
            for simpleChartRecord in chartReport:
                categoryName = simpleChartRecord.categoryName
                termName = simpleChartRecord.termName
                listHits = simpleChartRecord.listHits
                percent = simpleChartRecord.percent
                ease = simpleChartRecord.ease
                Genes = simpleChartRecord.geneIds
                listTotals = simpleChartRecord.listTotals
                popHits = simpleChartRecord.popHits
                popTotals = simpleChartRecord.popTotals
                foldEnrichment = simpleChartRecord.foldEnrichment
                bonferroni = simpleChartRecord.bonferroni
                benjamini = simpleChartRecord.benjamini
                FDR = simpleChartRecord.afdr
                rowList = [categoryName, termName, str(listHits), str(percent), str(ease), Genes, str(listTotals),
                           str(popHits), str(popTotals), str(foldEnrichment), str(bonferroni), str(benjamini), str(FDR)]
                fOut.write('\t'.join(rowList) + '\n')
        # print('write file:', resF, 'finished!')
        reports_file = outputfolder + 'report_file.txt'
        github_code = os.getcwd() + os.sep

        try:
            with open(reports_file, 'a') as runtime_report:
                subprocess.run(['Rscript', github_code + 'enrichment_david_plots.R', resF, category.replace("'", ""),
                                uniprotIds_file, outputfolder], check=True, stdout=runtime_report,
                               stderr=runtime_report)
        except Exception:
            print("Exception occured.")
