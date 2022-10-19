#!/bin/bash

# count reads at each step

#SBATCH --partition=main          # Partition (job queue)
#SBATCH --requeue                 # Return job to the queue if preempted
#SBATCH --job-name=counts         # Assign an short name to your job
#SBATCH --cpus-per-task=1         # Cores per task (>1 if multithread tasks)
#SBATCH --mem-per-cpu=4000M       # Real memory (RAM) required
#SBATCH --time=01:00:00           # Total run time limit (HH:MM:SS)
#SBATCH --error=/scratch/jsf149/kiledjian_2/circ_data/errs/slurm.%N.%j.err
#SBATCH --out=/dev/null
#SBATCH --mail-type=all
#SBATCH --mail-user=john.favate@rutgers.edu

# output location
OUT=/scratch/jsf149/kiledjian_2/circ_data/read_counts_per_step.tsv

# count reads for mate 1 
for FILE in /scratch/jsf149/kiledjian_2/circ_data/seqdata/1-original/*R1*.gz; do
  echo $FILE `zcat $FILE | awk 'NR%4==0' | wc -l` >> $OUT
done

for FILE in /scratch/jsf149/kiledjian_2/circ_data/seqdata/2-cleaned/*R1*.gz; do
  echo $FILE `zcat $FILE | awk 'NR%4==0' | wc -l` >> $OUT
done

for FILE in /scratch/jsf149/kiledjian_2/circ_data/seqdata/3-merged/*.gz; do
  echo $FILE `zcat $FILE | awk 'NR%4==0' | wc -l` >> $OUT
done

for FILE in /scratch/jsf149/kiledjian_2/circ_data/seqdata/4-selected/*.gz; do
  echo $FILE `zcat $FILE | awk 'NR%4==0' | wc -l` >> $OUT
done