# Preparing and analysing 16S data using mothur and the R phyloseq-package 

## mothur (V 1.39.5 in conda environment)

In this workflow the goal is to start from raw, demultiplexed amplicon reads to primer- and adapterfree, qualitytrimmed, chimera free amplicon sequences, which are counted and taxonomically annotated. We can do all this within mothur.

A possible way to set up directories, download required data bases and link your data to the work directory is described in `00_set_up_folder_and_start_mothur.sh`

To test mothur and for references, have a look into `00_katis_mothur_skript.h`. But note that she used the primer set 515f - 806r with V2-500 cycle reads.

The script for my glyphosate data including three separate sequencing runs is located under `01_mothur_workflow_1395.h`. I also stored the summary.seqs() information, which allows to track the development of your sequences for each processes (V1.39.5 and 1.41.0). This information can be found in `summary_seqs_1395_collected.txt`

## phyloseq

Later parts of the workflow are taken from the dada2 workflow (which we also ran) described here: https://f1000research.com/articles/5-1492/v2

After we prepared, counted and annotated our amplicon data with mothur, it is time to get the data into R. For this we use the phyloseq package, as it allows us to start the exactly same workflow as we would with e.g. dada2 or qiime data. phyloseq transforms all data into a standardized format, so once we have it there, we can use the same commands for amplicon data from different sources.

How to import the data to phyloseq is described in `02_import_mothur_to_phyloseq.r`. It will make use of mothur output files such as the .shared-file and the .constax-file. We will then add a file with the meta or sample data and with the representative sequence per OTU. Based on this data we can calculate a tree, which completes the phyloseq-object. This script also contains the calls to generate relative abundances. We will generate a phyloseq-melt-object, which contains all meta data, taxonomy and OTU counts/abundance and is suited for ggplot2.

The following scripts are based on the phyloseq-objects or the phyloseq-melt-objects to target specific questions:

* evaluate the distribution of OTUs per genus, their abundance and diversity `02_overview_on_abundant_OTUs_and_genera.r`

* bar-plot the microbial composition on order level (water column, DNA and RNA) over all time points `03_community_overview_plot.r`

* plot the abundance of a specific OTU over time for water column and biofilm, treatment and control, DNA and RNA `04_*_OTU_abundance_plots.r`

* ordination plots `missing`

* DESeq2-Test for differentially abundant OTUs before and after glyphosate addition `missing`

* Heatmap showing phylogenetic relation and reaction to glyphosate `missing`