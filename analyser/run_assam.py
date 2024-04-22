#!/usr/bin/env python
# coding: utf-8

"""
Python code to run ASSAM [1] on calculated motifs (using the generate_assam_motifs.py code).

Code developed inside the Virtuous Projct (https://virtuoush2020.com/)

----------------
Acknowledgements
This project has received funding from the European Union’s Horizon 2020 research and innovation programme under the Marie Skłodowska-Curie grant agreement Action (GA No. 872181)
----------
References
[1] Nadzirin, N., E. J. Gardiner, P. Willett, P. J. Artymiuk, and M. Firdaus-Raih. ‘SPRITE and ASSAM: Web Servers for Side Chain 3D-Motif Searching in Protein Structures’. Nucleic Acids Research 40, no. W1 (1 July 2012): W380–86. https://doi.org/10.1093/nar/gks401.
----------------
Version history:
- Version 1.0 - 26/01/2023
----------------
"""

__version__ = '1.0'
__author__ = 'Marco Cannariato, Lampros Androutsos, Eric Zizzi, Lorenzo Pallante'

import argparse
import os
import subprocess
import sys
import warnings

from motif_search import *

warnings.filterwarnings('ignore')


def parse_options():
    parser = argparse.ArgumentParser(description='''Tool to run ASSAM)''')
    parser.add_argument('--dbFolder',
                        metavar='db',
                        dest='db',
                        type=str,
                        help='Folder containing the pdbs to be screened for similar motifs')
    parser.add_argument('--motifsFolder',
                        metavar='motifs',
                        dest='motifs',
                        type=str,
                        help='Folder containing the motifs',
                        required=True)
    parser.add_argument('--outdir',
                        metavar='od',
                        dest='outdir',
                        type=str,
                        help='Path to the output directory.',
                        required=True)
    parser.add_argument('--version',
                        action='version',
                        version='%(prog)s ' + __version__)
    parser.add_argument('--author',
                        action='version',
                        version='%(prog)s ' + __author__)
    parser.add_argument('--verbose',
                        action='store_true',
                        help='Prints verbose output.')
    parser.add_argument('--cpu',
                        metavar='cpus',
                        dest='cpus',
                        type=int,
                        help='Number of cpus to be used for parallel version of Assam search.'
                             'Default is 1.',
                        default=1)

    args, unknown = parser.parse_known_args()
    expected = ["dbFolder", "motifsFolder", "outdir"]

    for cmd in unknown:
        if cmd not in expected:
            print('\n\n{0}\nUnknown command found in your command line: "{1}".\n{0}\n\n'.format('*' * 50, cmd))
            sys.exit(1)

    return args


def _run_assam(block, args):
    # define folder of the present script, which contains also the assam.sh script
    script_dir = os.path.dirname(os.path.realpath(__file__))

    out_dir = args.outdir + os.sep + block

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    pdb_dir = args.db + os.sep + block + os.sep + 'pdb'
    vek_dir = args.db + os.sep + block + os.sep + 'vek'
    mas_dir = args.db + os.sep + block

    # Use provided vek files
    # print('Running ASSAM with vek and MASTER.CURRENT files from ' + args.db + block + '...')
    with open(out_dir + os.sep + 'assam_run_' + block + '.log', 'w') as f_obj:

        if args.verbose:
            res = subprocess.call(
                ['bash', script_dir + os.sep + 'lib' + os.sep + 'assam_linux_v2' + os.sep + 'assam.sh', '-p', pdb_dir,
                 '-k', vek_dir, '-m', mas_dir, '-t', args.motifs, '-o', out_dir, '-v'], stdout=f_obj, stderr=f_obj)
        else:
            res = subprocess.call(
                ['bash', script_dir + os.sep + 'lib' + os.sep + 'assam_linux_v2' + os.sep + 'assam.sh', '-p', pdb_dir,
                 '-k', vek_dir, '-m', mas_dir, '-t', args.motifs, '-o', out_dir], stdout=f_obj, stderr=f_obj)

        if res != 0:
            print('Something went wrong running assam.sh. Check file {0}.'.format(
                out_dir + os.sep + 'assam_run_' + block + '.log'))
            sys.exit(1)


def main(args):
    '''Main function to run the script.'''
    # define folder of the present script which contains also the assam.sh script
    script_dir = os.path.dirname(os.path.realpath(__file__))

    # check if the output folder exists, if not create it
    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)

    #  Control if the proper db folder structure
    # db must contain a seried of folders (blocks), each one containing a pdb folder
    # if a vek folder and a MASTER.CURRENT file are provided use them, otherwise generate them
    vek_to_create = False

    if args.db:
        if not os.path.exists(args.db):
            print('The provided database folder does not exist.')
            sys.exit(1)
        else:
            # check if the db folder contains the proper structure
            for block in sorted(os.listdir(args.db)):
                if not os.path.isdir(args.db + os.sep + block):
                    print(
                        'The provided database folder contains a file. Please provide a folder containing a series of folders (blocks).')
                    sys.exit(1)
                else:
                    if not os.path.exists(args.db + os.sep + block + os.sep + 'pdb'):
                        print('The provided database folder contains a block folder without a pdb folder.')
                        sys.exit(1)
                    elif (not os.path.exists(args.db + os.sep + block + os.sep + 'vek')) & (
                            not os.path.exists(args.db + os.sep + block + os.sep + 'MASTER.CURRENT')):
                        vek_to_create = True
                    elif (os.path.exists(args.db + os.sep + block + os.sep + 'vek')) & (
                            not os.path.exists(args.db + os.sep + block + os.sep + 'MASTER.CURRENT')):
                        print(args.db + os.sep + block + os.sep + 'MASTER.CURRENT')
                        print(
                            'The provided database folder contains a block folder with a vek folder but without a MASTER.CURRENT file.')
                        sys.exit(1)
                    elif (not os.path.exists(args.db + os.sep + block + os.sep + 'vek')) & (
                            os.path.exists(args.db + os.sep + block + os.sep + 'MASTER.CURRENT')):
                        print(
                            'The provided database folder contains a block folder without a vek folder but with a MASTER.CURRENT file.')
                        sys.exit(1)

    #  if the user provided only a pdb folder, you need also to create vek and master files before running assam
    if vek_to_create:
        #  Run the script to generate the vek files
        print('Generating vek files from pdb files...')
        # invoking database.sh script
        # bash /home/lorenzo/Documenti/GitHub/CBP/code/assam_linux_v2/database.sh -p pdb2compare/ -o vek2compare/ -m mast2compare/ 
        for block in sorted(os.listdir(args.db)):
            # create the vek folder in the same location of the output folder
            out_dir = args.outdir + os.sep + block
            vek_dir = out_dir + os.sep + 'vek'
            pdb_dir = args.db + os.sep + block + os.sep + 'pdb'
            if not os.path.exists(vek_dir):
                os.makedirs(vek_dir)

            with open(out_dir + os.sep + 'assam_database.log', 'w') as f_obj:
                if args.verbose:
                    res = subprocess.call(
                        ['bash', script_dir + os.sep + 'lib' + os.sep + 'assam_linux_v2' + os.sep + 'database.sh',
                         '-p', pdb_dir, '-o', vek_dir, '-m', out_dir, '-v'], stdout=f_obj, stderr=f_obj)
                else:
                    res = subprocess.call(
                        ['bash', script_dir + os.sep + 'lib' + os.sep + 'assam_linux_v2' + os.sep + 'database.sh',
                         '-p', pdb_dir, '-o', vek_dir, '-m', out_dir], stdout=f_obj, stderr=f_obj)

                if res != 0:
                    print('Error in generating vek files. Check file {0} for more details.'.format(
                        out_dir + os.sep + 'assam_database.log'))
                    sys.exit(1)

            print('Done!')

            # Run ASSAM
            # invoking assam.sh script
            print('Running ASSAM...')
            with open(out_dir + os.sep + 'assam_run.log', 'w') as f_obj:
                if args.verbose:
                    res = subprocess.call(
                        ['bash', script_dir + os.sep + 'lib' + os.sep + 'assam_linux_v2' + os.sep + 'assam.sh', '-p',
                         pdb_dir, '-k', vek_dir, '-m', out_dir, '-t', args.motifs, '-o', out_dir, '-v'], stdout=f_obj,
                        stderr=f_obj)
                else:
                    res = subprocess.call(
                        ['bash', script_dir + os.sep + 'lib' + os.sep + 'assam_linux_v2' + os.sep + 'assam.sh', '-p',
                         pdb_dir, '-k', vek_dir, '-m', out_dir, '-t', args.motifs, '-o', out_dir], stdout=f_obj,
                        stderr=f_obj)

                if res != 0:
                    print('Something went wrong running assam.sh. Check file {0} file.',
                          out_dir + os.sep + 'assam_run.log')
                    sys.exit(1)

            print('Done!')

    else:
        # Old Serial version
        # for block in sorted(os.listdir(args.db)):
        #     out_dir = args.outdir + os.sep + block
        #     if not os.path.exists(out_dir):
        #         os.makedirs(out_dir)
        #     pdb_dir = args.db + os.sep + block + os.sep + 'pdb'
        #     vek_dir = args.db + os.sep + block + os.sep + 'vek'
        #     mas_dir = args.db + os.sep + block 
        #     # Use provided vek files
        #     print('Running ASSAM with vek and MASTER.CURRENT files from ' + args.db + block + '...')
        #     with open(out_dir + os.sep + 'assam_run.log','w') as f_obj:

        #         if args.verbose:
        #             res = subprocess.call(['bash', script_dir + os.sep + 'assam_linux_v2' + os.sep + 'assam.sh', '-p', pdb_dir, '-k', vek_dir, '-m', mas_dir, '-t', args.motifs, '-o', out_dir, '-v'],stdout=f_obj, stderr=f_obj)
        #         else:
        #             res = subprocess.call(['bash', script_dir + os.sep + 'assam_linux_v2' + os.sep + 'assam.sh', '-p', pdb_dir, '-k', vek_dir, '-m', mas_dir, '-t', args.motifs, '-o', out_dir],stdout=f_obj, stderr=f_obj)

        #         if res != 0:
        #             print('Something went wrong running assam.sh. Check file {0}.'.format(out_dir + os.sep + 'assam_run.log'))
        #             sys.exit(1)         

        # parallelize the assam run
        Parallel(n_jobs=args.cpus)(delayed(_run_assam)(block, args) for block in sorted(os.listdir(args.db)))
        print('Done!')


def run_script():
    args = parse_options()
    main(args)


if __name__ == '__main__':
    run_script()
