#!/bin/bash

threads=1
single_paired="pair"
adapter="NexteraPE-PE.fa:2:30:10"
sliding_window="5:20"
min_len="50"
max_memory="8G"
seq_type="fa"
sample_name=""
ref_genome=""
forward_file=""
reverse_file=""
output_dir="$(pwd)"

usage() {
  echo "Usage (Single End): ./virus_discovery_pipeline.sh [options] -c reference_genome -s file"
  echo "Usage (Paired End): ./virus_discovery_pipeline.sh [options] -c reference_genome -f forward_file -r reverse_file"
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


mkdir -p "$output_dir"
cd "$output_dir"


#Trimmomatic
tmp=(${forward_file//./ })
fastq_files_f="${tmp[0]}.trimmed.fastq"
fastq_files_f_un="${tmp[0]}.untrimmed.fastq"
fastq_files_r=""
fastq_files_r_un=""
if [ "${single_paired}" = "single" ]; then
  TrimmomaticSE \
    -threads $threads \
    "${forward_file}" \
    "${fastq_files_f}" \
    $adapter \
    SLIDINGWINDOW:$sliding_window \
    MINLEN:$min_len
else
  tmp=(${reverse_file//./ })
  fastq_files_r="${tmp[0]}.trimmed.fastq"
  fastq_files_r_un="${tmp[0]}.untrimmed.fastq"
  TrimmomaticPE \
     -threads $threads \
    "${forward_file}" "${reverse_file}" \
    "${fastq_files_f}" "${fastq_files_f_un}" "${fastq_files_r}" "${fastq_files_r_un}" \
    $adapter \
    SLIDINGWINDOW:$sliding_window \
    MINLEN:$min_len
fi


#fastq -> fasta
tmp=(${fastq_files_f//./ })
fasta_files_f="${tmp[0]}.trimmed.fasta"
sed -n '1~4s/^@/>/p;2~4p' ${fastq_files_f} > "${fasta_files_f}"
sed -i 's#\.1\s#/1 #g' ${fasta_files_f}
rm "${fastq_files_f}" 
if [ "${single_paired}" = "pair" ]; then
  tmp=(${fastq_files_f_un//./ })
  fasta_files_f_un="${tmp[0]}.untrimmed.fasta"
  sed -n '1~4s/^@/>/p;2~4p' ${fastq_files_f_un} > "${fasta_files_f_un}"
  sed -i 's#\.1\s#/1 #g' ${fasta_files_f_un}
  rm "${fastq_files_f_un}"
  tmp=(${fastq_files_r//./ })
  fasta_files_r="${tmp[0]}.trimmed.fasta"
  sed -n '1~4s/^@/>/p;2~4p' ${fastq_files_r} > "${fasta_files_r}"
  sed -i 's#\.2\s#/2 #g' ${fasta_files_r}
  tmp=(${fastq_files_r_un//./ })
  fasta_files_r_un="${tmp[0]}.untrimmed.fasta"
  sed -n '1~4s/^@/>/p;2~4p' ${fastq_files_r_un} > "${fasta_files_r_un}"
  sed -i 's#\.2\s#/2 #g' ${fasta_files_r_un}
  rm "${fastq_files_r}" "${fastq_files_r_un}"
fi


#Trinity
trinity_dir="./${sample_name}_trinity"
if [ "${single_paired}" = "single" ]; then
  Trinity --seqType $seq_type \
    --max_memory $max_memory \
    --single "${fasta_files_f}" \
    --no_normalize_reads \
    --CPU $threads \
    --output "${trinity_dir}"
else
  Trinity --seqType $seq_type \
    --max_memory $max_memory \
    --left "${fasta_files_f},${fasta_files_f_un}" \
    --right "${fasta_files_r},${fasta_files_r_un}" \
    --CPU $threads \
    --output "${trinity_dir}"
fi


#BWA-MEM
bwa index $ref_genome
bwa mem $ref_genome "${trinity_dir}.Trinity.fasta" > "./output.sam"


#SAMTOOLS
#Getting the unmapped reads from a sam file:
samtools view -f 4 "./output.sam" > "./output_unmapped.sam"
#Getting only the mapped reads from a sam file:
samtools view -b -F 4 "./output.sam" > "./output_mapped.sam"
#Unmapped sequences from SAM to FASTA
samtools fasta "./output_unmapped.sam" > "./output_unmapped.fasta"

#BLAST
blastn -db nt -query "./output_unmapped.fasta" -out "./output_unmapped_BLASTn.txt" -max_target_seqs 5 -remote