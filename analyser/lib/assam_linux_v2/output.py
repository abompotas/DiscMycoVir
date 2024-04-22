#!/usr/bin/env python
# -*- coding: utf-8 -*-
#A program to parse ASSAM output file (.LPL, .LPS)

# 23/01/2023 - Modified by LP to define pdb and bin directories


import sys, os, re, argparse, csv
import numpy
from string import ascii_letters

trans = {'ALA':'A','CYS':'C','CYH':'C','CSS':'C','ASP':'D','GLU':'E','PHE':'F','GLY':'G','HIS':'H','ILE':'I','LYS':'K','LEU':'L','MET':'M','ASN':'N','PRO':'P','GLN':'Q','ARG':'R','SER':'S','THR':'T','VAL':'V','TRP':'W','TYR':'Y','UNK':'X','MSE':'M'}

def generate_results(filedir,side,pdbdir,bindir):
    if side=="R":
        filename = "AUTO.LPS"
        transformed_path = "transformr"
        csvoutput = "sumr.csv"
    if side=="L":
        filename = "AUTO.LPL"
        transformed_path = "transforml"
        csvoutput = "suml.csv"
    patternstry = re.findall(r"KERBm.*",open(filename,"r").read())
    if patternstry !=[]:
        print((filedir, side, "generate results.. "))
        patterns = open(filename,"r").read().split("KERBm")
        if not os.path.exists(transformed_path): os.makedirs(transformed_path)
        csv_writer = csv.writer(open(csvoutput,"w"),delimiter=",")
        for bil,pat in enumerate(patterns[1:]):
            pdbid = re.findall(r"FOUND in [0-9][a-zA-Z0-9][a-zA-Z0-9][a-zA-Z0-9]",pat)[0].replace("FOUND in ","").strip()
            macromolecule = re.findall(r"FOUND in.*size =.*",pat)[0].split("^")[-1].split("#")[0].lstrip().rstrip().replace(";","")
            rmsd = re.findall(r"RMSD =.*",pat)[0].replace("RMSD =     ","")
            matches1 = [line.split("matches")[-1].replace(">","").replace("<","").lstrip().rstrip() for line in re.findall(r"<.*matches.*>",pat)]
            matches = [line.split()[-1]+" "+line[:5] for line in matches1]
            hetatm = [line.split("(")[1].split("from")[0].lstrip().rstrip()+" ( "+line.split("(")[1].split("from")[1].lstrip() if len(line.split("("))>1 else "None" for line in re.findall(r"<.*matches.*>.*",pat)]
            #hetresidues1 = [het.split("(")[1].replace(")","").rstrip().lstrip() for het in hetatm if het !="None"]
            hetresidues1 = [het.split("(")[1].replace(")","").rstrip().lstrip() for het in hetatm if het !="None" if "site" not in het if re.findall(r"^\d+\.\d+",het) !=[] if float(re.findall(r"^\d+\.\d+",het)[0])<5.0 ]
            hetresidues = sorted(list(set([het.replace(het.split()[0],"").lstrip().rstrip(ascii_letters) for het in hetresidues1])))
            matrix = re.findall(r"TRANSFORM.*",pat)[0].replace("TRANSFORM","").split()
            resno = str(len(matches))
            transformed_id = side+"pat-"+str(bil+1)
            results = [filedir]+[transformed_id]+[pdbid]+[macromolecule]+[resno]+["; ".join(matches)]+["; ".join(hetatm)]+[rmsd]+[",".join(matrix)]
            csv_writer.writerow(results)
            hitpdb = transformed_path+"/"+side+"patt-"+str(bil+1)+".pdb"
            pdbfile = pdbdir+pdbid+".pdb"
            hitlines = [ line for line in open(pdbfile,"r").readlines() if line[:4]=="ATOM" if line[17:26] in matches ]+[ line for line in open(pdbfile,"r").readlines() if line[:6]=="HETATM" if line[17:26] in hetresidues ]
            with open(hitpdb,"w") as output: output.write("".join(hitlines))
            with open("{}hebpat.sh".format(side),"a") as matmultout:
                if os.path.isfile(hitpdb):
                    matrix_rot = [" "*(8-len(s))+s for s in matrix[:9]]
                    matrix_trans = [" "*(4)+" "*(8-len(s))+s for s in matrix[9:]]
                    transformedlist = " "+" ".join(matrix_rot)+"".join(matrix_trans)
                    matmultinput = ["#hit no {}".format(bil+1),"", bindir + "/matmult1r <<EOF","hebe.pdb","","{}/{}hebe-{}.pdb".format(transformed_path,side,str(bil+1)),transformedlist,"Y","EOF","","# end"]
                    matmultout.write("\n".join(matmultinput))
    else:
        print((filedir, side, "No hits found.."))

#in ucsf chimera: color blue #0.1; color red #0.2; color yellow ligand

# Create the parser
my_parser = argparse.ArgumentParser(description='Get first model from NMR structure')

# Add the arguments
my_parser.add_argument('filedir',type=str,help='Pattern name')
my_parser.add_argument('side',type=str,help='R/L superposition')
my_parser.add_argument('pdbdir',type=str,help='PDB directory')
my_parser.add_argument('bindir',type=str,help='ASSAM Bin directory')


# Execute the parse_args() method
args = my_parser.parse_args()
print ("Get results from the output file")

generate_results(args.filedir,args.side, args.pdbdir, args.bindir)

