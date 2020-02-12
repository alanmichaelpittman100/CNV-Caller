# CNV-Caller
## Exome CNV caller in R

Automated script to run ExomeDepth in R to call CNVs in 2 or more Exomes

Requires the following dependencies/libraries installed under bioconductor:

library("ExomeDepth") \
library("GenomeInfoDb") \
library("Rsamtools") 

Requres pre-aligned .bam files.

Requires config file with header listing the .bam files you want to analyse \
save as a config.csv file 

list_of_bam_files \
1_realigned.bam \
2_realigned.bam \
3_realigned.bam \
3_realigned.bam

run script in directory with BAM files (and BAI) and setwd() of the directory in which you want to perform the analysis.

ExomeDepth_withConfig_V2.R : is depreciated and not used

CNV-Calling-With-Exome-Depth.Rmd : is active and for use in in interactice session and for troubleshooting. In RMarkdown Format

ExomeDepth_withConfig_V3.R : is active and runs atomatically with config.csv file ; Rscript ExomeDepth_withConfig_V3.R



Includes python script to join outputs from multiple samples into a single excel file for ease of viewing 


