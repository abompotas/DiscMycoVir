# %%
from Bio import Entrez
import numpy as np
import datetime
from tqdm import tqdm
import concurrent.futures
import time

start_time = time.time()

date = datetime.datetime.now().strftime("%Y%m%d")
out_folder = '/home/lorenzo/Documenti/GitHub/CBP/Files/'
out_file = out_folder + 'PDB_ids_' + date + '.csv'

Entrez.email = "marco.cannariato@polito.it"
record = Entrez.read(Entrez.einfo())
handle = Entrez.esearch(db='structure',term='Human[Organism]',retmax=100000)
record = Entrez.read(handle)
ids = record['IdList']
handle.close()

pdbs = list()

def fetch_pdb(id):
    handle = Entrez.esummary(db='structure',id=id)
    record = Entrez.read(handle)
    handle.close()
    return record[0]['PdbAcc']

# serial version
for i in tqdm(ids):
    pdbs.append(fetch_pdb(i))

# save ids in a csv file
np.savetxt(out_file,np.array(pdbs),fmt="%s")

# get elapsed time
elapsed_time = time.time() - start_time
print('Elapsed time: ' + "%.2f" %elapsed_time + ' seconds')

