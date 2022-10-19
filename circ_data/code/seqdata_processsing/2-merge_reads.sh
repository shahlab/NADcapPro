#!/bin/bash

#SBATCH --partition=main          # Partition (job queue)
#SBATCH --requeue                 # Return job to the queue if preempted
#SBATCH --job-name=merge          # Assign an short name to your job
#SBATCH --cpus-per-task=2         # Cores per task (>1 if multithread tasks)
#SBATCH --mem-per-cpu=4G          # Real memory (RAM) required (MB)
#SBATCH --time=03:00:00           # Total run time limit (HH:MM:SS)
#SBATCH --error=/scratch/jsf149/kiledjian_2/circ_data/errs/slurm.%N.%j.err
#SBATCH --output=/dev/null
#SBATCH --array=0-19
#SBATCH --mail-type=all
#SBATCH --mail-user=john.favate@rutgers.edu

# the plan is to operate on read 1, read 2, and merges separately. This merges reads
# and saves them to the processed directory

source /home/jsf149/miniconda3/bin/activate && conda activate bbmap_env

module load java/14.0.1

# echo a list of packages used in this environment to a file
conda list > /scratch/jsf149/kiledjian_2/conda_envs/bbmap_env.txt

# define the array of files
ARRAY=( $(
# for each file
for file in /scratch/jsf149/kiledjian_2/circ_data/seqdata/2-cleaned/*.gz; do
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
echo $SLURM_ARRAY_JOB_ID\_$SLURM_ARRAY_TASK_ID $FILE >> /scratch/jsf149/kiledjian_2/circ_data/errs/errs_files.txt

# define the pairs
PAIR1=/scratch/jsf149/kiledjian_2/circ_data/seqdata/2-cleaned/$FILE\_R1_001.fastq.gz
PAIR2=/scratch/jsf149/kiledjian_2/circ_data/seqdata/2-cleaned/$FILE\_R2_001.fastq.gz
  
# Merged output
OUT=/scratch/jsf149/kiledjian_2/circ_data/seqdata/3-processed/$FILE\.fq.gz
  
# merge stats
STATS=/scratch/jsf149/kiledjian_2/circ_data/reports/$FILE/$FILE\.txt.merge

# insert size histogram location
IHIST=/scratch/jsf149/kiledjian_2/circ_data/reports/$FILE/$FILE\.txt.ihist

bbmerge.sh in=$PAIR1 in2=$PAIR2 out=$OUT ihist=$IHIST t=2 -Xmx8G 2>$STATS