# Abstract

The widespread herbicide glyphosate has been detected in aquatic coastal zones of the southern Baltic Sea recently. We monitored community dynamics in glyphosate-impacted chemostats for a total of 20 weeks to evaluate its potential impact on free-living and biofilm-associated bacterial community assemblages in brackish ecosystems. The system was analyzed up to twice a week concerning glyphosate concentrations as well as prokaryotic succession by total cell counts and next generation 16S rRNA (gene) amplicon sequencing, respectively. Shotgun metagenomics revealed new insights into the glyphosate degradation potential of the microbial communities by evaluating abundance and sequence similarity of phn genes. Glyphosate caused temporally increases of total cell counts, bacterial diversity and abundances of distinct bacterial operational taxonomic units in the water column. Biofilm communities proved to be less affected than pelagic ones, but their responses were longer lasting. The increase of phn-operon gene abundance indicated glyphosate degradation by the sarcosine pathway, which was indirectly supported by the absence of aminomethylphosphonic acid. However, despite the fact that glyphosate concentrations were reduced by 99%, quantities of 1 µM remained till the end of the experiment. This indicates that at such low concentrations degradation in the ecosystem can be disadvantageous for bacteria; thus, glyphosate entering the Baltic Sea could also resist bacterial degradation at low concentrations.

## Scripts

The scripts can be found in the respective repos, to run the code please refer to the `installation_guides` folder

## 16S amplicon analysis

The program `mothur` (1.39.5) was used to process amplicon reads, which were then read into `R` package `phyloseq`. For processing we also tested the R package `dada2` (1.8) (together with `cutadapt` (1.8.3)), but did not analyze further on.

## Metagenomic analysis

We checked for glyphosate degradation relevant genes (`gox`, `soxABGD`, `phnC-P`) and their abundance and sequence similarity


