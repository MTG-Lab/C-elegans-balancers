# C. elegans Balancers

Scripts and data for used for figures on the manuscript.

karyoploteR.R: R script to reproduce the Figures 2A, 5B and 5D. 

karyoploteR_map.R: R script to reproduce the Figure 1A. Uses the data in mart_export_map.tsv (genes genomic position extarcted from Biomart).

coverage_analysis.sh: bash script implementing samtools dept to estimate the average of coverage by window of 1 kb, normalized by the avrge of coverage of the entire genome. 
