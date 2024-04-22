#!/bin/bash

set -e

# This script divides the databse folders into folders of 10000 files maximum.
pdbdir=/mnt/data/lorenzo/00-work/CBP/db/20230309/pdb_clean
vekdir=/mnt/data/lorenzo/00-work/CBP/db/20230309/vek_clean
masdir=/mnt/data/lorenzo/00-work/CBP/db/20230309/mas_clean
outdir=/mnt/data/lorenzo/00-work/CBP/db/20230309/blocks3

# base name of the input folders
pdbdir_base=$(basename $pdbdir)
vekdir_base=$(basename $vekdir)

# count number of files in the original folder
nfiles=$(ls $pdbdir | wc -l)

# dimension of the blocks
blocksize=6000

# decide number of blocks
nblocks=$(($nfiles/$blocksize+1))

# counter
cnt=1

# create output folders
for dir in $(seq -w 01 $nblocks); do

    echo "Creating folder $outdir/$dir-block$dir.."

    mkdir -p $outdir/$dir-block$dir/pdb/ 
    mkdir -p $outdir/$dir-block$dir/vek/ 

    # divide the files in the original folder into the new ones
    head_files=$(($blocksize*$cnt))

    if [ $cnt -eq $nblocks ]; then
        tail_files=$(($nfiles - $head_files + $blocksize))
    else
        tail_files=$blocksize
    fi
    
    echo "Copying pdb and vek files for block $cnt .."

    # cycle over files comprises between head_files and tail_files in the original folder
    for pdb in $(ls $pdbdir | head -n $head_files | tail -n $tail_files); do
        ln -s ../../../$pdbdir_base/$pdb $outdir/$dir-block$dir/pdb/
        ln -s ../../../$vekdir_base/${pdb%.*}.vek $outdir/$dir-block$dir/vek/
    done

    cat $masdir/MASTER.CURRENT | head -n 1 >$outdir/$dir-block$dir/MASTER.CURRENT
    cat $masdir/MASTER.CURRENT | tail -n +2 | head -n $head_files | tail -n $tail_files >>$outdir/$dir-block$dir/MASTER.CURRENT
    
    cnt=$(($cnt+1))

    echo "Done"

done        