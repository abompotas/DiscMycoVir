#!/bin/bash

threads=1
single_paired="pair"
sample_name=""
forward_file=""
reverse_file=""
output_dir="virus_discovery_output"

usage() {
  echo "Usage (Single End): ./analysis.sh [options] -s file"
  echo "Usage (Paired End): ./analysis.sh [options] -f forward_file -r reverse_file"
  echo "OPTIONS:"
  echo "-t	Number of threads for parallel execution"
  echo "-n	Sample name"
  echo "-s	Single end input file"
  echo "-f	Paired forward input file"
  echo "-r	Paired reverse input file"
  echo "-o	Output directory"
  exit 1
}


while getopts ":t:n:g:s:f:r:o:h:" option; do
  case $option in
    t)
      threads=$OPTARG
      ;;
    n)
      sample_name=$OPTARG
      ;;
    g)
      ref_genome=$OPTARG
      ;;
    s)
      single_paired="single"
      forward_file=$OPTARG
      ;;
    f)
      single_paired="pair"
      forward_file=$OPTARG
      ;;
    r)
      single_paired="pair"
      reverse_file=$OPTARG
      ;;
    o)
      output_dir=$OPTARG
      ;;
    h)
      usage
      ;;
    *)
      usage
      ;;
  esac
done


if [ "${forward_file}" = "" ]; then
  usage
fi
if [ "${single_paired}" = "pair" ]; then
  if [ "${reverse_file}" = "" ]; then
    usage
  fi
fi

current_dir=$(pwd)
mkdir -p "/tmp/${sample_name}"
cd "/tmp/${sample_name}"


#FastQC
tmp=(${forward_file//./ })
fastqc_out_f="${tmp[0]}_fastqc"
if [ "${single_paired}" = "single" ]; then
  fastqc -o ./fastqc_analysis "${forward_file}" -t $threads -q
else
  tmp=(${reverse_file//./ })
  fastqc_out_r="${tmp[0]}_fastqc"
  fastqc -o ./fastqc_analysis "${forward_file}" "${reverse_file}" -t $threads -q
fi


# Copy files and cleanup
mkdir -p "${output_dir}"
mv fastqc_analysis "$output_dir"
cd "${current_dir}"
rm -r "/tmp/${sample_name}"