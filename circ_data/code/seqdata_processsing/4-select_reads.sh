#!/bin/bash

#SBATCH --partition=main          # Partition (job queue)
#SBATCH --requeue                 # Return job to the queue if preempted
#SBATCH --job-name=selection      # Assign an short name to your job
#SBATCH --cpus-per-task=4         # Cores per task (>1 if multithread tasks)
#SBATCH --mem-per-cpu=2G          # Real memory (RAM) required (MB)
#SBATCH --time=01:00:00           # Total run time limit (HH:MM:SS)
#SBATCH --error=/scratch/jsf149/kiledjian_2/circ_data/errs/slurm.%N.%j.err
#SBATCH --output=/dev/null
#SBATCH --mail-type=all
#SBATCH --mail-user=john.favate@rutgers.edu

# get only reads that have an AAAAAAAA in them 
source /home/jsf149/miniconda3/bin/activate && conda activate cutadapt_env

# echo a list of packages used in this environment to a file
conda list > /scratch/jsf149/kiledjian_2/conda_envs/cutadapt_env.txt

for FILE in /scratch/jsf149/kiledjian_2/circ_data/seqdata/3-processed/*.gz; do
    # output location
    OUT=/scratch/jsf149/kiledjian_2/circ_data/seqdata/4-selected/`basename $FILE`
    
    # dont do mito genes
    if [[ "$FILE" == *"COX"* ]]; then
        cp $FILE $OUT
    elif [[ "$FILE" == *"21S"* ]]; then
        cp $FILE $OUT
    else
        # rev comp the file
        cutadapt -j 4 -a AAAAAAAA --discard-untrimmed --action=none -o $OUT $FILE
    fi
done
