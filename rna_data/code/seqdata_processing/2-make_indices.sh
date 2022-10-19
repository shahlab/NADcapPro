#!/bin/bash

# Given the correct fastas, this makes the rRNA/tRNA, and transcript indices for 
# hisat2/kallisto

#SBATCH --partition=main          # Partition (job queue)
#SBATCH --requeue                 # Return job to the queue if preempted
#SBATCH --job-name=indexing       # Assign an short name to your job
#SBATCH --cpus-per-task=2         # Cores per task (>1 if multithread tasks)
#SBATCH --mem-per-cpu=4000M       # Real memory (RAM) required
#SBATCH --time=01:00:00           # Total run time limit (HH:MM:SS)
#SBATCH --error=/dev/null
#SBATCH --out=/dev/null

module load gcc/5.4

source /home/jsf149/miniconda3/bin/activate && conda activate hisat2_env

# hisat2 indices
for FILE in /scratch/jsf149/kiledjian_2/fastas/*.fa; do
    INAME=`basename $FILE | cut -d '.' -f 1`

    OUT=/scratch/jsf149/kiledjian_2/rna_data/alignment/hisat2/indices/$INAME

    hisat2-build -p 2 $FILE $OUT
done

# kallisto indices, only the transcripts
conda deactivate && conda activate kallisto_env

# echo a list of packages used in this environment to a file
conda list > /scratch/jsf149/kiledjian_2/conda_envs/kallisto_env.txt

kallisto index -i /scratch/jsf149/kiledjian_2/rna_data/alignment/kallisto/indices/GCF_000146045.2_R64_rna.kidx /scratch/jsf149/kiledjian_2/fastas/GCF_000146045.2_R64_rna.fa