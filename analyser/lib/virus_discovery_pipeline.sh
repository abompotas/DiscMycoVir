#!/bin/bash

threads=1
single_paired="pair"
adapter="NexteraPE-PE.fa:2:30:10"
sliding_window="5:20"
min_len="50"
max_memory="8G"
seq_type="fq"
sample_name=""
ref_genome=""
forward_file=""
reverse_file=""
output_dir="virus_discovery_output"

usage() {
  echo "Usage (Single End): ./virus_discovery_pipeline.sh [options] -g reference_genome -s file"
  echo "Usage (Paired End): ./virus_discovery_pipeline.sh [options] -g reference_genome -f forward_file -r reverse_file"
  echo "OPTIONS:"
  echo "-t	Number of threads for parallel execution"
  echo "-a	Adapter for Trimmomatic"
  echo "-w	Sliding window for Trimmomatic"
  echo "-l	Minimum length for Trimmomatic"
  echo "-q	Sequence type to be used by Trinity"
  echo "-m	Max memory to be used by Trinity"
  echo "-n	Sample name"
  echo "-g	Reference genome"
  echo "-s	Single end input file"
  echo "-f	Paired forward input file"
  echo "-r	Paired reverse input file"
  echo "-o	Output directory"
  exit 1
}


while getopts ":t:i:w:l:m:q:n:c:s:f:r:h:" option; do
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
    q)
      seq_type=$OPTARG
      ;;
    m)
      max_memory=$OPTARG
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


if [ "${ref_genome}" = "" ]; then
  usage
fi


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
  # rm forward
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
    # rm forward reverse
fi


if [ "${seq_type}" = "fa" ]; then
  tmp=(${trimmomatic_out_f//./ })
  trinity_in_f="${tmp[0]}.trimmed.fasta"
  sed -n '1~4s/^@/>/p;2~4p' ${trimmomatic_out_f} > "${trinity_in_f}"
  sed -i 's#\.1\s#/1 #g' ${trinity_in_f}
  rm "${trimmomatic_out_f}"
  if [ "${single_paired}" = "pair" ]; then
    tmp=(${trimmomatic_out_f_un//./ })
    trinity_in_f_un="${tmp[0]}.untrimmed.fasta"
    sed -n '1~4s/^@/>/p;2~4p' ${trimmomatic_out_f_un} > "${trinity_in_f_un}"
    sed -i 's#\.1\s#/1 #g' ${trinity_in_f_un}
    rm "${trimmomatic_out_f_un}"
    tmp=(${trimmomatic_out_r//./ })
    trinity_in_r="${tmp[0]}.trimmed.fasta"
    sed -n '1~4s/^@/>/p;2~4p' ${trimmomatic_out_r} > "${trinity_in_r}"
    sed -i 's#\.2\s#/2 #g' ${trinity_in_r}
    tmp=(${trimmomatic_out_r_un//./ })
    trinity_in_r_un="${tmp[0]}.untrimmed.fasta"
    sed -n '1~4s/^@/>/p;2~4p' ${trimmomatic_out_r_un} > "${trinity_in_r_un}"
    sed -i 's#\.2\s#/2 #g' ${trinity_in_r_un}
    rm "${trimmomatic_out_r}" "${trimmomatic_out_r_un}"
  fi
else
  trinity_in_f="${trimmomatic_out_f}"
  trinity_in_f_un="${trimmomatic_out_f_un}"
  trinity_in_r="${trimmomatic_out_r}"
  trinity_in_r_un="${trimmomatic_out_r_un}"
fi


#Trinity
trinity_dir="./${sample_name}_trinity"
if [ "${single_paired}" = "single" ]; then
  Trinity --seqType $seq_type \
    --max_memory $max_memory \
    --single "${trinity_in_f}" \
    --no_normalize_reads \
    --CPU $threads \
    --output "${trinity_dir}"
else
  Trinity --seqType $seq_type \
    --max_memory $max_memory \
    --left "${trinity_in_f},${trinity_in_f_un}" \
    --right "${trinity_in_r},${trinity_in_r_un}" \
    --CPU $threads \
    --output "${trinity_dir}"
fi


#BWA-MEM
bwa index $ref_genome
bwa mem $ref_genome "${trinity_dir}.Trinity.fasta" > ./output.sam

#SAMTOOLS
#Getting the unmapped reads from a sam file:
samtools view -f 4 output.sam > output_unmapped.sam
#Getting only the mapped reads from a sam file:
samtools view -b -F 4 output.sam> output_mapped.sam
#Unmapped sequences from SAM to FASTA
samtools fasta output_unmapped.sam > output_unmapped.fasta

#BLAST
blastn -db nt -query output_unmapped.fasta -out output_unmapped_BLASTn.txt -max_target_seqs 5 -remote

# Copy files and cleanup
mkdir -p "$output_dir"
mv output_unmapped.sam output_mapped.sam output_unmapped.fasta output_unmapped_BLASTn.txt "$output_dir"
mv "${trinity_in_f}" "${trinity_dir}.Trinity.fasta" "$output_dir"
if [ "${single_paired}" = "pair" ]; then
  mv "${trinity_in_f_un}" "${trinity_in_r}" "${trinity_in_r_un}" "$output_dir"
fi
cd "${current_dir}"
rm -r "/tmp/${sample_name}"