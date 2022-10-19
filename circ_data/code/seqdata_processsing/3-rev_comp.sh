#!/bin/bash

#SBATCH --partition=main          # Partition (job queue)
#SBATCH --requeue                 # Return job to the queue if preempted
#SBATCH --job-name=revcomp        # Assign an short name to your job
#SBATCH --cpus-per-task=1         # Cores per task (>1 if multithread tasks)
#SBATCH --mem-per-cpu=8G          # Real memory (RAM) required (MB)
#SBATCH --time=01:00:00           # Total run time limit (HH:MM:SS)
#SBATCH --error=/scratch/jsf149/kiledjian_2/circ_data/errs/slurm.%N.%j.err
#SBATCH --output=/dev/null
#SBATCH --mail-type=all
#SBATCH --mail-user=john.favate@rutgers.edu

# read 2 can get reverse complemented because we're going to be treating them as single
# end reads and it's easier this way, activate bbmap conda env
module load java/14.0.1

source /home/jsf149/miniconda3/bin/activate && conda activate bbmap_env

# reverse complement these files 
for FILE in /scratch/jsf149/kiledjian_2/circ_data/seqdata/2-cleaned/*R2*.gz; do
    # output location
    OUT=/scratch/jsf149/kiledjian_2/circ_data/seqdata/3-processed/`basename $FILE`
    
    # rev comp the file
    reformat.sh -Xmx8G rcomp=t in=$FILE out=$OUT
done

# deactivate conda env
conda deactivate

# also copy the read 1 files over, untouched
cp /scratch/jsf149/kiledjian_2/circ_data/seqdata/2-cleaned/*R1* /scratch/jsf149/kiledjian_2/circ_data/seqdata/3-processed/