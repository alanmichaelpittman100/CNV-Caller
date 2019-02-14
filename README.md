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

Includes python script to join outputs from multiple samples into a single excel file for easy of viewing 
