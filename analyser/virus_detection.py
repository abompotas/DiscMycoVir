import argparse
import multiprocessing
import os
import shutil
import subprocess
import sys
# time the code
import time


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

    parser = argparse.ArgumentParser()
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


def run_pipeline(args):
    print(args.keys())
    if args['single_paired'] == "pair" and not (args["reverse_file"] and args["forward_file"]):
        print("ERROR: MISSING Revers and forward files")
    #step 1: TrimmomaticSE (paired-single)
    #step 2: fastq -> fasta
    #step 3: Trinity (paired-single)
    #step 4: BWA-MEM
    #step 5: SAMTOOLS
    #step 6: BLAST


# run main
if __name__ == '__main__':
    start_time = time.time()

    # clear screen
    os.system('clear')

    # parse arguments
    args = parse_options()

    # run the CBP tool
    run_pipeline(args)

    # print elapsed time
    print('\nElapsed time: %.2f s' % (time.time() - start_time))
