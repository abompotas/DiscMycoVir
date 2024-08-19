#!/bin/bash

discovery_dir="virus_discovery_output"

usage() {
  echo "Usage (Single End): ./discovery_blast.sh [options] -g reference_genome -s file"
  echo "Usage (Paired End): ./discovery_blast.sh [options] -g reference_genome -f forward_file -r reverse_file"
  echo "OPTIONS:"
  echo "-d	Discovery directory"
  exit 1
}


while getopts ":d:h:" option; do
  case $option in
    d)
      discovery_dir=$OPTARG
      ;;
    h)
      usage
      ;;
    *)
      usage
      ;;
  esac
done


blastn -db nt_viruses -max_target_seqs 5 -outfmt 5 \
  -query "$discovery_dir/output_unmapped.fasta" \
  -out "$discovery_dir/output_blast.xml" \
  -remote

touch "$discovery_dir/blast_finished"