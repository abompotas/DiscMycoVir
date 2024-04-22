#assam.sh : Amino acid 3D pattern searching in protein structures
#Directories:
#    - patterns: stores user-defined patterns (3-12 residues)
#    - pdb: stores PDB-formatted protein structures
#    - VEC_D: stores matrices of protein structures for comparison

# version edited by LP
# v.2.0 - 2023-01-13: change folders input and output
version="2.0"

script_name="ASSAM Motifs Search - v $version"

# exit on error
set -e

usage="
$script_name

USAGE: assam.sh -p <PDB directory input> -k <VEK directory input> -t <PATTERNS folder> -o <resutlts output folder> -m <MASTER.CURRENT> -v <verbose>

"

# define the script directory
script_path=$(dirname $(readlink -f $0))

# define default path for MASTER.CURRENT
master_path=$(pwd)

# define default verbose option
verbose=false

# parsing input arguments using getopts
while getopts ':hp:k:t:o:m:v' option; do
  case "$option" in
    h) echo "$usage"
       exit
       ;;
    p) pdb_folder=$(readlink -m $OPTARG); 
       ;;
    k) vek_folder=$(readlink -m $OPTARG);
         ;;
    t) patterns_folder=$(readlink -m $OPTARG);  
         ;;
    o) results_folder=$(readlink -m $OPTARG);
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
    echo "ERROR: The input directory with pdbs does not exist"
    echo "$usage"
    exit 1
fi

if [ ! -d $patterns_folder ]; then
    echo "ERROR: The input directory with patterns to be compared does not exist"
    echo "$usage"
    exit 1
fi

if [ ! -d $vek_folder ]; then
    echo "ERROR: The input directory with vek files does not exist"
    echo "$usage"
    exit 1
fi

if [ ! -d $master_path ]; then
    echo "ERROR: The input directory with MASTER.CURRENT does not exist"
    echo "$usage"
    exit 1
fi

# Check if the output directory exists and create it if not
if [ ! -d $results_folder ]; then
    mkdir $results_folder
fi 

#set environment
if command -v export &> /dev/null
then
    export ASS_EXE="$script_path/bin"
    export MAS_LIST="$master_path/MASTER.CURRENT"
    export VEC_D="$vek_folder"
    export VMC_D="$(pwd)/VMC_D"
    export PIK_D="$(pwd)/PIK_D"
    export PDB_D="$pdb_folder"
    export NEW_D="$(pwd)/NEW_D"
else
    setenv ASS_EXE="$script_path/bin"
    setenv MAS_LIST="$master_path/MASTER.CURRENT"
    setenv VEC_D="$vek_folder"
    setenv VMC_D="$(pwd)/VMC_D"
    setenv PIK_D="$(pwd)/PIK_D"
    setenv PDB_D="$pdb_folder"
    setenv NEW_D "$(pwd)/NEW_D"
fi

#run ASSAM
for file in $patterns_folder/*.pdb; do
    echo "Processing $file"
    # copy file to directory
    bname=$(basename "$file")
    filenamewoext="${bname%.*}"
    cwd=$(pwd)
    newdir="$results_folder/results/${bname%.*}"
    newfilename="${newdir}/hebe.pdb"
    mkdir -p $newdir
    cp ${file} ${newfilename}

    #run ASSAM
    cd $newdir
    $ASS_EXE/aamc_pat < $script_path/bin/aamc_pat.INP
    $ASS_EXE/aspsettings
    (ulimit -t 900; $ASS_EXE/aspreg2f < AUTOASP.INP) | head -c 1073741824
    #right-handed superposition (output: AUTO.LPS)
    $ASS_EXE/aspsup < $script_path/bin/aspsup_R.INP
    #left-handed superposition (output: AUTO.LPL)
    $ASS_EXE/aspsup < $script_path/bin/aspsup_L.INP
    python $script_path/output.py ${filenamewoext} R ${pdb_folder}/ ${script_path}/bin
    python $script_path/output.py ${filenamewoext} L ${pdb_folder}/ ${script_path}/bin
    
    # execute Rhebpat.sh and Lhebpat.sh only if they exist
    # if they do not exist, it means that no matches were found
    if [ -f Rhebpat.sh ]; then
        if [ $verbose = true ]; then
            sh Rhebpat.sh
        else
            sh Rhebpat.sh >/dev/null 2>&1
        fi

        #get transformed files
        for file in transformr/Rhebe*.pdb; do
            newfile="${file/hebe/pat}"
            pattfile="${file/hebe/patt}"
            echo "MODEL        1\n" >> $newfile
            cat $file >> $newfile
            echo "END\nMODEL        2\n" >> $newfile
            cat $pattfile >> $newfile
            rm $file
            rm $pattfile
        done

    fi
    
    if [ -f Lhebpat.sh ]; then
        if [ $verbose = true ]; then
            sh Lhebpat.sh
        else
            sh Lhebpat.sh >/dev/null 2>&1
        fi

        #get transformed files
        for file in transforml/Lhebe*.pdb; do
            newfile="${file/hebe/pat}"
            pattfile="${file/hebe/patt}"
            echo "MODEL        1\n" >> $newfile
            cat $file >> $newfile
            echo "END\nMODEL        2\n" >> $newfile
            cat $pattfile >> $newfile
            rm $file
            rm $pattfile
        done
    fi

    # back to original directory
    cd $cwd

done