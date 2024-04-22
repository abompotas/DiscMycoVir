
# Create Database with human proteome PDB
python PDB_fetch.py
python PDB_download.py 
# cleaning DB
python /home/lorenzo/Documenti/GitHub/CBP/code/PDB_clean.py --inDir /home/lorenzo/00-work/CBP/database/20230309/pdb/ --outDir /home/lorenzo/00-work/CBP/database/20230309/pdb_clean/ --errDir /home/lorenzo/00-work/CBP/database/20230309/
# fix missing atoms and residues, add protonation, etc
python /home/lorenzo/Documenti/GitHub/CBP/code/PDB_fix.py --inDir /home/lorenzo/00-work/CBP/database/20230309/pdb_clean/ --outDir /home/lorenzo/00-work/CBP/database/20230309/pdb_fix/ --errDir /home/lorenzo/00-work/CBP/database/20230309/
# --> not Done! Too much time to be completed --> fix pdb only before docking

# Convert to vek and master files for ASSAM (optional)
bash /home/lorenzo/Documenti/GitHub/CBP/code/assam_linux_v2/database.sh -p pdb_clean/ -o vek_clean/ -v -m mas_clean/

# Divide into blocks
bash /home/lorenzo/Documenti/GitHub/CBP/code/PDB_divideblocks.sh

# run the pipeline (master.py)
1 - generate assam motifs
2 - run assam 
3 - analyse assam results 
