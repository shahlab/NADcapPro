This repository contains the code associated with the article XXXXXX.
It can be divided into two parts, code that deals with the standard RNAseq data analysis and code that deals with the NAD based circularization assay.
If you're attempting to reproduce analyses, the first thing you should do is clone this repository.
Then, you'll need to acquire the raw sequencing data from XXXXXX and place it in the correct folders.
Processing of the sequencing data was done on a HPC system running CentOS and SLURM. 
If this is not the configuration of your system, you'll need to modify the code to suit your needs.

# RNAseq data

## Generate an index for alignment

After being cleaned, the data will be aligned with kallisto.
This requires an index, you need to acquire the fastas for this first.
Instructions for doing so can be found in the script `rna/code/make_indices.Rmd`.
They should end up in the correct directories on their own.

## Process the RNAseq data

Various conda environments were used to process the sequencing data, the details of these environments can be seen in the files in the `conda_envs` directory.
Setting up identical environments is recommended.

1. Download the data from XXXXXX and place it in the folder at `/rna/seqdata/1-original`.
2. Run the scripts in `/rna/code/seqdat_processing` in the indicated order. These will
   1. Remove the adapters and clean the data
   2. Make indices if you have not already.
   3. Run an *in silico* rRNA/tRNA depletion, this reduces file sizes by removing rRNA and tRNA reads, which we're not interested in.
   4. Run hisat2 to generate SAM files
   5. Convert said SAMs to BAMs
   6. Run kallisto to get transcript quantifications.
   
The result should be a bunch of directories in `/alignment/kallisto/output` which contain the read counts.
You should then run the code in `rna/code/kallisto_cleanup.Rmd` to collect and clean these files.

## Differential expression

After you get the counts, you can run the code in `/rna/code/deseq2.Rmd`, this will get you the differential expression results.
The HTML of the same name shows the versions of R packages used in this and other code.

## Figures

You can then run the code in `rna/code/figures.Rmd` or `figures_no_labs.Rmd` (same code, just produces figures with no labels), and `supp_figs.Rmd` to make the figures.
They should end up in the figures folder.

# NAD based circularization assay

1. As before, download the raw sequencing data from XXXXXX and place it in `circ_data/seqdata/1-original`.
2. As before, run the scripts in `circ_data/coded/seqdata_processing` in the indicated order
   1. Remove the adapters and clean the data
   2. Merge the PE reads to a single read, but keep the pairs as well
   3. reverse complement read 2 so it's easier to deal with
   4. filter reads based on some criteria
   5. count the unique reads
   6. count reads at each step so you can see how many reads you're left with after step
   
After that, you can run the `circ_data/code/msas.Rmd` script to make MSAs that guided the manual analysis.
