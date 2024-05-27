#!/bin/bash

threads=1
single_paired="pair"
adapter="NexteraPE-PE.fa:2:30:10"
sliding_window="5:20"
min_len="50"
sample_name=""
forward_file=""
reverse_file=""
output_dir="virus_discovery_output"

usage() {
  echo "Usage (Single End): ./trimming.sh [options] -s file"
  echo "Usage (Paired End): ./trimming.sh [options] -f forward_file -r reverse_file"
  echo "OPTIONS:"
  echo "-t	Number of threads for parallel execution"
  echo "-a	Adapter for Trimmomatic"
  echo "-w	Sliding window for Trimmomatic"
  echo "-l	Minimum length for Trimmomatic"
  echo "-n	Sample name"
  echo "-s	Single end input file"
  echo "-f	Paired forward input file"
  echo "-r	Paired reverse input file"
  echo "-o	Output directory"
  exit 1
}


while getopts ":t:a:w:l:n:s:f:r:o:h:" option; do
  case $option in
    t)
      threads=$OPTARG
      ;;
    a)
      adapter=$OPTARG
      ;;
    w)
      sliding_window=$OPTARG
      ;;
    l)
      min_len=$OPTARG
      ;;
    n)
      sample_name=$OPTARG
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


#Trimmomatic
tmp=(${forward_file//./ })
trimmomatic_out_f="${tmp[0]}.trimmed.fastq"
trimmomatic_out_f_un="${tmp[0]}.untrimmed.fastq"
trimmomatic_out_r=""
trimmomatic_out_r_un=""
if [ "${single_paired}" = "single" ]; then
  TrimmomaticSE \
    -threads $threads \
    "${forward_file}" \
    "${trimmomatic_out_f}" \
    $adapter \
    SLIDINGWINDOW:$sliding_window \
    MINLEN:$min_len
else
  tmp=(${reverse_file//./ })
  trimmomatic_out_r="${tmp[0]}.trimmed.fastq"
  trimmomatic_out_r_un="${tmp[0]}.untrimmed.fastq"
  TrimmomaticPE \
     -threads $threads \
    "${forward_file}" "${reverse_file}" \
    "${trimmomatic_out_f}" "${trimmomatic_out_f_un}" "${trimmomatic_out_r}" "${trimmomatic_out_r_un}" \
    $adapter \
    SLIDINGWINDOW:$sliding_window \
    MINLEN:$min_len
fi


# Copy files and cleanup
mkdir -p "$output_dir/trimming"
mv "${trimmomatic_out_f}" "$output_dir/trimming"
if [ "${single_paired}" = "pair" ]; then
  mv "${trimmomatic_out_f_un}" "${trimmomatic_out_r}" "${trimmomatic_out_r_un}" "$output_dir/trimming"
fi
cd "${current_dir}"
rm -r "/tmp/${sample_name}"