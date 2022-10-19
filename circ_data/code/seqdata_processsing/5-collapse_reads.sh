#!/bin/bash

#SBATCH --partition=main   # Partition (job queue)
#SBATCH --requeue                 # Return job to the queue if preempted
#SBATCH --job-name=dumb_counts    # Assign an short name to your job
#SBATCH --cpus-per-task=1         # Cores per task (>1 if multithread tasks)
#SBATCH --mem-per-cpu=4000M       # Real memory (RAM) required
#SBATCH --time=01:00:00           # Total run time limit (HH:MM:SS)
#SBATCH --array=0-59              # n array members
#SBATCH --error=/scratch/jsf149/kiledjian_2/circ_data/errs/slurm.%N.%j.err
#SBATCH --out=/dev/null
#SBATCH --mail-type=all
#SBATCH --mail-user=john.favate@rutgers.edu

# Make each file a member of an array
ARRAY=($(find /scratch/jsf149/kiledjian_2/circ_data/seqdata/4-selected -name "*.gz"))

# assign each file to a slurm array number
FILE=${ARRAY[$SLURM_ARRAY_TASK_ID]}

FNAME=`basename $FILE | cut -d '.' -f 1`

# output location
OUT=/scratch/jsf149/kiledjian_2/circ_data/unique_counts/$FNAME\.txt

# count the unique reads in a file
zcat $FILE | awk 'NR%4==2' | sort | uniq -c | sort -rg > $OUT