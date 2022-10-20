#!/bin/bash

#SBATCH --partition=p_shahlab_1   # Partition (job queue)
#SBATCH --requeue                 # Return job to the queue if preempted
#SBATCH --job-name=rrna_dep       # Assign an short name to your job
#SBATCH --cpus-per-task=4         # Cores per task (>1 if multithread tasks)
#SBATCH --mem-per-cpu=2000        # Real memory (RAM) required (MB)
#SBATCH --time=02:00:00           # Total run time limit (HH:MM:SS)
#SBATCH --error=/scratch/jsf149/kiledjian_2/errs/slurm.%N.%j.err
#SBATCH --output=/dev/null
#SBATCH --array=0-19

source /home/jsf149/miniconda3/bin/activate && conda activate hisat2_env

module load gcc/5.4

# define the array of files
ARRAY=( $(
# for each file
for file in /scratch/jsf149/kiledjian_2/rna_data/seqdata/2-ar/*.gz; do
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
INDEX=/scratch/jsf149/kiledjian_2/rna_data/alignment/hisat2/indices/rtrna

# define the input read pairs
PAIR1=/scratch/jsf149/kiledjian_2/rna_data/seqdata/2-ar/$FILE\_R1_001.fastq.gz
PAIR2=/scratch/jsf149/kiledjian_2/rna_data/seqdata/2-ar/$FILE\_R2_001.fastq.gz

# define the unaligned i.e. rRNA depleted file locations. The % symbol
# will cause insertion of a 1 or 2 depending on which file the read
# originates from
UN=/scratch/jsf149/kiledjian_2/rna_data/seqdata/3-depleted/$FILE\_R%_001.fastq.gz

# run it
hisat2 -x $INDEX -1 $PAIR1 -2 $PAIR2 -p 4 \
  --un-conc-gz $UN --no-unal --no-mixed --no-discordant --no-spliced-alignment \
  --rna-strandness FR -S /dev/null