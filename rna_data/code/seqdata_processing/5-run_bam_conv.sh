#!/bin/bash

#SBATCH --partition=main   # Partition (job queue)
#SBATCH --requeue                 # Return job to the queue if preempted
#SBATCH --job-name=bam_conv       # Assign an short name to your job
#SBATCH --cpus-per-task=4         # Cores per task (>1 if multithread tasks)
#SBATCH --mem-per-cpu=3G          # Real memory (RAM) required
#SBATCH --time=01:00:00           # Total run time limit (HH:MM:SS)
#SBATCH --array=0-8              # n array members
#SBATCH --error=/scratch/jsf149/kiledjian_2/rna_data/errs/slurm.%N.%j.err
#SBATCH --out=/dev/null

source /home/jsf149/miniconda3/bin/activate

conda activate samtools_env

# Make each SAM a member of an array
ARRAY=($(find /scratch/jsf149/kiledjian_2/rna_data/alignment/hisat2/output/ -name "*.sam"))

# assign each file to a slurm array number
FILE=${ARRAY[$SLURM_ARRAY_TASK_ID]}

# file name sans extension
FNAME=`basename $FILE | cut -d '.' -f 1`

# send the job ID to file to mathch error reports and file names
echo $SLURM_ARRAY_JOB_ID\_$SLURM_ARRAY_TASK_ID $FILE >> /scratch/jsf149/kiledjian_2/rna_data/errs/errs_files.txt

# new BAM file location
BAM=/scratch/jsf149/kiledjian_2/rna_data/alignment/hisat2/output/$FNAME\.bam

# sort and convert to bam, -F 16 excludes - strand reads
samtools view -F 16 -b $FILE | \
    samtools sort \
    -T "/scratch/jsf149/kiledjian_2/rna_data/alignment/hisat2/output/$FNAME" -@ 4 -o $BAM -

# index it too
INDEX=/scratch/jsf149/kiledjian_2/rna_data/alignment/hisat2/output/$FNAME\.bai

samtools index $BAM $INDEX

# also do flagstat
STATS=/scratch/jsf149/kiledjian_2/rna_data/reports/$FNAME\/$FNAME\.txt.flagstats

samtools flagstat $BAM > $STATS

# remove the sam
rm $FILE