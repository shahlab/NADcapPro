#!/bin/bash

#SBATCH --partition=main          # Partition (job queue)
#SBATCH --requeue                 # Return job to the queue if preempted
#SBATCH --job-name=hisat2         # Assign an short name to your job
#SBATCH --cpus-per-task=16        # Cores per task (>1 if multithread tasks)
#SBATCH --mem-per-cpu=500M        # Real memory (RAM) required
#SBATCH --time=01:00:00           # Total run time limit (HH:MM:SS)
#SBATCH --array=0-19              # n array members
#SBATCH --error=/scratch/jsf149/kiledjian_2/errs/slurm.%N.%j.err
#SBATCH --out=/dev/null

module load gcc/5.4

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

# index to use
INDEX=/scratch/jsf149/kiledjian_2/alignment/hisat2/indices/GCF_000146045.2_R64_rna

# define the input read pairs
PAIR1=/scratch/jsf149/kiledjian_2/seqdata/3-depleted/$FILE\_R1_001.fastq.gz
PAIR2=/scratch/jsf149/kiledjian_2/seqdata/3-depleted/$FILE\_R2_001.fastq.gz

# SAM output
SAM=/scratch/jsf149/kiledjian_2/alignment/hisat2/output/$FILE\.sam

# stderr output
ERR=/scratch/jsf149/kiledjian_2/reports/$FILE\/$FILE\.txt.hisat2

# run it
/home/jsf149/bin/hisat2-2.2.1/hisat2 -x $INDEX -1 $PAIR1 -2 $PAIR2 -p 16 -S $SAM \
  --no-unal --no-mixed --no-discordant --no-spliced-alignment --rna-strandness FR 2>$ERR
