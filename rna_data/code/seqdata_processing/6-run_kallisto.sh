#!/bin/bash

#SBATCH --partition=main    # Partition (job queue)
#SBATCH --requeue           # Return job to the queue if preempted
#SBATCH --job-name=kallisto # Assign an short name to your job
#SBATCH --cpus-per-task=8   # Cores per task (>1 if multithread tasks)
#SBATCH --mem-per-cpu=1000  # Real memory (RAM) required (MB)
#SBATCH --time=00:30:00     # Total run time limit (HH:MM:SS)
#SBATCH --array=0-19        # slurm array members
#SBATCH --error=/scratch/jsf149/kiledjian_2/rna_data/errs/slurm.%N.%j.err
#SBATCH --output=/dev/null
              
# define the array of files
ARRAY=( $(
# for each file
for file in /scratch/jsf149/kiledjian_2/rna_data/seqdata/3-depleted/*.gz; do
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

# send the job ID to file to mathch error reports and file names
echo $SLURM_ARRAY_JOB_ID\_$SLURM_ARRAY_TASK_ID $FILE >> /scratch/jsf149/kiledjian_2/rna_data/errs/errs_files.txt

# make the output directory, one per pair
OUTPUT=/scratch/jsf149/kiledjian_2/rna_data/alignment/kallisto/output/$FILE
  
mkdir -p $OUTPUT
  
# define the read pairs
PAIR1=/scratch/jsf149/kiledjian_2/rna_data/seqdata/3-depleted/$FILE\_R1_001.fastq.gz
PAIR2=/scratch/jsf149/kiledjian_2/rna_data/seqdata/3-depleted/$FILE\_R2_001.fastq.gz
  
# denote the index
INDEX=/scratch/jsf149/kiledjian_2/rna_data/alignment/kallisto/indices/GCF_000146045.2_R64_rna.kidx

# stderr from kallisto
ERR=/scratch/jsf149/kiledjian_2/rna_data/reports/$FILE\/$FILE\.txt.kallisto
  
/home/jsf149/bin/kallisto quant -i $INDEX -o $OUTPUT -t 8 --bias $PAIR1 $PAIR2 2>$ERR
