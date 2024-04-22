#!/usr/bin/env python
# coding: utf-8

# code to clean the DB of PDB preserving only one model per pdb

import os
from tqdm import tqdm
import warnings
import argparse
import numpy as np
import subprocess
import glob
import re
import fileinput
import shutil

# suppress warnings
warnings.filterwarnings('ignore')

def convert_to_list(string):
    if ', ' in string:
        return string.split(', ')
    else:
        return [string]


# parsing options
parser = argparse.ArgumentParser(description='''Tool to clean the DB of PDB preserving only one model per pdb''')
parser.add_argument('--inDir',
                    metavar='inDir',
                    dest='inDir',
                    type=str,
                    help='Input directory containing PDB files.',
                    required=True)
parser.add_argument('--outDir',
                    metavar='od',
                    dest='outdir',
                    type=str,
                    help='Path to the output directory.',
                    required=True)
parser.add_argument('--errDir',
                    metavar='ef',
                    dest='errDir',
                    type=str,
                    help='Path to the error directory.',
                    required=True)
args = parser.parse_args()


# define the directories
start_directory = args.inDir # directory where the pdbs are
end_directory = args.outdir # directory where the pdbs will be saved
error_folder = args.errDir # directory where the error files will be saved

# create directory if not exists
if not os.path.exists(end_directory):
    os.makedirs(end_directory)

# create list for pdb files with error in cleaning
error_list = list()
altloc_error_list = list()
error_file = os.path.join(error_folder, 'error_clean.txt')
altloc_error_file = os.path.join(error_folder, 'altloc_error.txt')

# cycle over all files in the original directory containing the pdbs
for filename in tqdm(os.listdir(start_directory)):

    # define the path of the pdb
    pdb_file = os.path.join(start_directory, filename)
    out_file = os.path.join(end_directory, filename)

    ##########################################
    # STEP 1 - REMOVE NON HUMAN PARTS from PDB
    
    # initialize the pdb dictionary and useful variables
    pdb_dict = {}
    molid = []
    organism = []
    chain = []
    header_chain  = list()
    header_organism = list()

    # read the pdb file
    with open(pdb_file, 'r') as f:
        for line in f.readlines():
            
            # get molid, chain and organism
            if line.startswith('COMPND'):
                if 'MOL_ID' in line:
                    header_chain.append(line)
                    molid.append(line.replace(';', ':').split(':')[1].strip())
                if 'CHAIN:' in line:
                    header_chain.append(line)   
                    chain.append(line.replace(';', ':').split(':')[1].strip())

            # get the organism
            if line.startswith('SOURCE'):
                if 'MOL_ID' in line:
                    header_organism.append(line)
                if 'ORGANISM_TAXID' in line:
                    header_organism.append(line)

    # check the correspondance between molid/chain and molid/organism
    for i in np.argwhere(['MOL_ID' in x for x in header_organism]).ravel():
        try: 
            if 'ORGANISM_TAXID' in header_organism[i+1]:
                organism.append(header_organism[i+1].replace(';', ':').split(':')[1].strip())
            else:
                # if the organism is not specified, we assume it is human
                organism.append('9606')
        except IndexError:
            organism.append('9606')

    # modify the pdb file if it contains chains not belonging to human proteome
    if np.any(organism != '9606'): 
        # create dictionary
        pdb_dict['molid'] = molid
        pdb_dict['chain'] = chain
        pdb_dict['organism'] = organism

        # remove chains not belonging to human proteome
        not_human_chains = []
        for i in range(len(pdb_dict['organism'])):
            if pdb_dict['organism'][i] != '9606':
                not_human_chains.append(pdb_dict['chain'][i])

        # concatenate chains
        not_human_chains = ', '.join(not_human_chains)
        # convert to list
        not_human_chains = convert_to_list(not_human_chains)
        
        # remove chains not belonging to human proteome
        with open(pdb_file, 'r') as f:
            with open(out_file, 'w') as f_out:
                for line in f.readlines():
                    # line starting with ATOM or HETATM
                    if line.startswith(tuple(['ATOM', 'HETATM', 'TER', 'ANISOU'])):
                        if line[21] not in not_human_chains:
                            f_out.write(line)      
                            res = line[21]              
                    else:                
                        f_out.write(line)
    
    ################################################################################
    # STEP 2 - REMOVE PDB FILES with no ATOM or HETATM left from the previous step
    # -> mainly protein chimeras!
    f = open('pdbinfo.txt','w')
    # Count number of atoms
    process = subprocess.run(['pdb_wc', '-a', out_file,],
                                stdout=f)
    f.close()

    line = open('pdbinfo.txt', 'r').readlines()[0]
    
    # if no atoms are zero, remove the pdb
    if line.split()[2] == '0':
        os.remove(out_file)
        os.remove('pdbinfo.txt')
        error_list.append(filename)
        continue

    os.remove('pdbinfo.txt')


    #########################################
    # STEP3 - REMOVE MULTIPLE MODELS from PDB
    # read the pdb
    with open(out_file, "r") as f:
        lines = f.readlines()

    frame = []
    model_counter = 0

    for line in lines:
        if line.startswith("MODEL"):

            model_counter += 1
            if model_counter > 1:
                break
            
        frame.append(line)

    with open(out_file, "w") as f:
        f.writelines(frame)

    #################################################
    # STEP4 - REMOVE ANISOU AND HETATM lines from PDB

    # read the pdb
    lines = open(out_file,'r').readlines()
    f = open(out_file,'w')
    for line in lines:
        if not line.startswith('ANISOU') and not line.startswith('HETATM'):
            f.write(line)
    f.close()


    #################################
    # STEP5 - REMOVE ALTLOC from PDB

    noaltLoc_file = os.path.join(end_directory, 'noaltLoc.pdb')

    # select atoms with highest occupacy if necessary
    f = open('pdbinfo.txt','w')
    process = subprocess.run(['pdb_wc', '-o', out_file], stdout=f)
    f.close()

    # check if alternative locations are present
    flag_altloc = open('pdbinfo.txt', 'r').readlines()[0].split()[2] # True if altloc are present
    os.remove('pdbinfo.txt')
    
    if flag_altloc == 'True':
        #  remove alternative locations if present with pdb_tools
        f = open(noaltLoc_file,'w')
        process = subprocess.run(['pdb_selaltloc', out_file], stdout=f, stderr=subprocess.PIPE)
        f.close()

        # check exit status
        if process.returncode != 0:
            
            altloc_error_list.append(filename)
            # discard the pdb
            os.remove(out_file)
            os.remove(noaltLoc_file)
            continue
    else:
        # if no alternative locations are present, copy the pdb
        shutil.copy(out_file, noaltLoc_file)
    
    

    #################
    # WRITE OUT file
    shutil.copy(noaltLoc_file, out_file)
    os.remove(noaltLoc_file)


# print the list of pdb files with error in cleaning
open(error_file, 'w').write('\n'.join(error_list))
open(altloc_error_file, 'w').write('\n'.join(altloc_error_list))