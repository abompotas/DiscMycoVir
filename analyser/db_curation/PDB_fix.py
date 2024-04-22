#!/usr/bin/env python
# coding: utf-8

# code to fix the DB of PDB adding missing atoms and residues, adding hydrigens, etc..

import argparse
import os
from tqdm import tqdm
import sys
from pdbfixer import PDBFixer
from openmm.app import PDBFile

def fix_pdb(pdbfile, outFile, verbose=False):

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
            print("INFO: PDB file is %s" %pdbpath)
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
    with open(pdbpath,'r') as f:
        lines = f.readlines()
    # Remove all lines containing HOH:
    lines = [l for l in lines if 'HOH' not in l]
    # Write new PDB file:
    with open(outFile,'w') as f:
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
    PDBFile.writeFile(fixer.topology, fixer.positions, open(outFile, 'w+'))

    if verbose:
        print("INFO: PDB file saved!")
        print("INFO: Done!")
    
    return outFile_path

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
error_file = os.path.join(error_folder, 'error_fix.txt')

# cycle over all files in the original directory containing the pdbs
for filename in tqdm(os.listdir(start_directory)):  

    pdbfile = os.path.join(start_directory, filename)
    outFile = os.path.join(end_directory, filename)

    try:
        fix_pdb (pdbfile, outFile, verbose=False)
    except:
        error_list.append(filename)
        
    # copy pdb header from original pdb file
    with open(pdbfile, 'r') as f:
        header_lines = f.readlines()
        header_lines = [l for l in header_lines if l.startswith(('HEADER', 'TITLE', 'COMPND', 'SOURCE', 'JRNL'))]

    # prepend the header to the fixed pdb file
    with open(outFile, 'r') as f:
        lines = f.readlines()
        # keep only lines starting with ATOM
        lines = [l for l in lines if l.startswith('ATOM')]

    with open(outFile, 'w+') as f:
        f.writelines(header_lines)
        f.writelines(lines)


# save the list of pdb files with error in cleaning
with open(error_file, 'w') as f:
    for item in error_list:
        f.write("%s" % item)
