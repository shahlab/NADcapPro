#!/bin/bash

#SBATCH --partition=main          # Partition (job queue)
#SBATCH --requeue                 # Return job to the queue if preempted
#SBATCH --job-name=fastp          # Assign an short name to your job
#SBATCH --cpus-per-task=8         # Cores per task (>1 if multithread tasks)
#SBATCH --mem-per-cpu=1G          # Real memory (RAM) required (MB)
#SBATCH --time=03:00:00           # Total run time limit (HH:MM:SS)
#SBATCH --error=/scratch/jsf149/kiledjian_2/errs/slurm.%N.%j.err
#SBATCH --output=/dev/null
#SBATCH --array=0-19

#SBATCH --mail-type=all
#SBATCH --mail-user=john.favate@rutgers.edu

source /home/jsf149/miniconda3/bin/activate

conda activate fastp_env

# define the array of files
ARRAY=( $(
# for each file
for file in /scratch/jsf149/kiledjian_2/pcr_seqdata/1-original/*.gz; do
  # get the file name sans path
  PARTS=`basename $file`
  
  # remove the part that denotes the pair
  PARTS=${PARTS//_R[12]_00[12].fastq.gz/}
  
  # get and push only unique elements to the array
  echo $PARTS
done | sort -u
) )

# assign each file stem to a slurm array number
FILE=${ARRAY[$SLURM_ARRAY_TASK_ID]}

# define the pairs
PAIR1=/scratch/jsf149/kiledjian_2/pcr_seqdata/1-original/$FILE\_R1_001.fastq.gz
PAIR2=/scratch/jsf149/kiledjian_2/pcr_seqdata/1-original/$FILE\_R2_001.fastq.gz
  
OUT1=/scratch/jsf149/kiledjian_2/pcr_seqdata/2-cleaned/$FILE\_R1_001.fastq.gz
OUT2=/scratch/jsf149/kiledjian_2/pcr_seqdata/2-cleaned/$FILE\_R2_001.fastq.gz
  
# reports files, these need to be named fastp.json and fastp.html for
# multiqc to see them but they cant all end up in the same place so they each
# get dirnames that are the sample names which contain fastp.json/html

# send the job ID to file to mathch error reports and file names
echo $SLURM_ARRAY_JOB_ID\_$SLURM_ARRAY_TASK_ID $FILE >> /scratch/jsf149/kiledjian_2/errs/errs_files.txt

# make a dir for each sample pair
mkdir -p /scratch/jsf149/kiledjian_2/reports/$FILE

# the locations of the jsons/htmls
JSON=/scratch/jsf149/kiledjian_2/reports/$FILE/fastp.json
HTML=/scratch/jsf149/kiledjian_2/reports/$FILE/fastp.html

# adapters
A1=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
A2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT

fastp -i $PAIR1 -I $PAIR2 -o $OUT1 -O $OUT2 --adapter_sequence $A1 --adapter_sequence_r2 $A2 -w 16 -j $JSON -h $HTML -R "$FILE"