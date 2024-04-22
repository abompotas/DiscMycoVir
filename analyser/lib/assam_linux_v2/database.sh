#!/bin/bash

#convert pdb to vek
# version edited by LP
# v.2.0 - 2023-01-13: change folders input and output

version="2.0"

script_name="ASSAM Conversion - PDB to VEC - v $version"

usage="
$script_name

USAGE: database.sh -p <PDB directory input> -o <VEK directory output> -m <MASTER.CURRENT output directory> -v <verbose>

"

# define the script directory
script_path=$(dirname $(readlink -f $0))

# define default path for MASTER.CURRENT
master_path=$(pwd)

# define default verbose option
verbose=false

# parsing input arguments using getopts
while getopts ':hp:o:m:v' option; do
  case "$option" in
    h) echo "$usage"
       exit
       ;;
    p) pdb_folder=$(readlink -m $OPTARG); 
       ;;
    o) vek_folder=$(readlink -m $OPTARG);
      ;;
    m) master_path=$(readlink -m $OPTARG);
        ;;
    v) verbose=true
        ;;
    :) printf "\nFATAL: missing argument for -%s\n" "$OPTARG" >&2
       echo -e "use -h for help\n" >&2
       exit 1
       ;;
   \?) printf "\nFATAL: illegal option: -%s\n" "$OPTARG" >&2
       echo -e "use -h for help\n" >&2
       exit 1
       ;;
  esac
done
shift $((OPTIND - 1))

# Check if the input directory exists
if [ ! -d $pdb_folder ]; then
    echo "ERROR: The input directory does not exist"
    echo "$usage"
    exit 1
fi

# Check if the output directory exists and create it if not
if [ ! -d $vek_folder ]; then
    mkdir $vek_folder
fi 

# Check if the MASTER.CURRENT directory exists and create it if not
if [ ! -d $master_path ]; then
    mkdir $master_path
fi

# Create vector file directories, delete directories if exist
[[ -d VMC_D ]] && rm -r VMC_D
mkdir -p VMC_D

# Get the list of PDB files
ls -1 $pdb_folder > pdb.list

#set environment
if command -v export &> /dev/null
then
    export ASS_EXE="$script_path/bin"
    export MAS_LIST="$master_path/MASTER.CURRENT"
    export VEC_D="$vek_folder"
    export VMC_D="$(pwd)/VMC_D"
    export PIK_D="$(pwd)/PIK_D"
    export PDB_D="$pdb_folder"
else
    setenv ASS_EXE="$script_path/bin"
    setenv MAS_LIST="$master_path/MASTER.CURRENT"
    setenv VEC_D="$vek_folder"
    setenv VMC_D="$(pwd)/VMC_D"
    setenv PIK_D="$(pwd)/PIK_D"
    setenv PDB_D="$pdb_folder"
fi

# Run pdbaamc where it will compile the function for individual PDB file
$ASS_EXE/pdbaamc   <<EOF
pdb.list
EOF

printf -- "\nConverting PDB to VEC\n"
# if verbose is true print the output of the script otherwise redirect it to /dev/null
if [ "$verbose" == true ]; then
    source pdbaamc.com 
else 
    source pdbaamc.com >/dev/null 2>&1
fi

# Get the list of VEK files
ls -1 $VEC_D > vek.list

# Run mastermod to compile the list of vector files that will be used by ASSAM later
DATE=$(date +"%m-%d-%Y")
$ASS_EXE/mastermod  <<EOF
vek.list
MASTER.CURRENT
$DATE
EOF

# move the master file of the master path is different from the current path
[[ $master_path != $(pwd) ]] && mv MASTER.CURRENT $master_path/MASTER.CURRENT

#remove unused directory and files
rm -f chet.sh mas.bak newvek_lino6.sftp newvek_mac.sftp newvek_macPair.sftp newvek.ftp pdb.list pdbaamc.com pdbaamc.log vek.list MASTER.DNARNA
rm -rf VMC_D
