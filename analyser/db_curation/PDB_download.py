
#!/usr/bin/env python
# coding: utf-8

from Bio.PDB.PDBList import PDBList
import pandas as pd
from tqdm import tqdm
import urllib.request
import time

start_time = time.time()

# set the output folder
folder = '/home/lorenzo/00-work/CBP/database/pdb/'

# set the error file
error_file = '/home/lorenzo/00-work/CBP/database/error_download.txt'

# read csv file with pdb ids
x = '/home/lorenzo/Documenti/GitHub/CBP/Files/PDB_ids_20230201.csv'
df = pd.read_csv(x, header=None)
listt = df[0].values.tolist()

# download pdb files
for id in tqdm(listt):
    url = "https://files.rcsb.org/download/" + id.lower() + ".pdb"
    filename = folder + id.lower() + ".pdb"

    # download pdb from url
    try: 
        urllib.request.urlretrieve(url, filename)

    except:
        # save pdb code in a txt file
        with open(error_file, 'a') as f:
            f.write(id + '\n')

# get elapsed time
elapsed_time = time.time() - start_time