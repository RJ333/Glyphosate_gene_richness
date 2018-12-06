# Preparing and analysing 16S data using mothur and the R phyloseq-package 

## mothur (V 1.39.5 in conda environment)

In this workflow the goal is to start from raw, demultiplexed amplicon reads to primer- and adapterfree, qualitytrimmed, chimera free amplicon sequences, which are counted and taxonomically annotated. We can do all this within mothur.

A possible way to set up directories, download required data bases and link your data to the work directory is described in `00_set_up_folder_and_start_mothur.sh`

To test mothur and for references, have a look into `00_katis_mothur_skript.h`. But note that she used the primer set 515f - 806r with V2-500 cycle reads.

The script for my glyphosate data including three separate sequencing runs is located under `01_mothur_workflow_1395.h`. I also stored the summary.seqs() information, which allows to track the development of your sequences for each processes (V1.39.5 and 1.41.0). This information can be found in `summary_seqs_1395_collected.txt`

