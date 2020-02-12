### ExomeDepth - Alan Pittman -  - reusable_ Rscript - CNV caller Febbruary 2020

#requires index .bam files to run
#Use with the congifguration file (list number of samples basically)
#calls CNVs on sample set(~4-10) using rest as reference files

## R Vession 3.6.1
## ExomeDepth 1.1.10
## Genomic Ranges 1.32.7


#example config file with header (save as a config.csv file) 
################################################
#requires config file with header listing the .bam files you want to analyse
#save as a config.csv file 

#list_of_bam_files
#1_realigned.bam
#2_realigned.bam
#3_realigned.bam
#3_realigned.bam


library(ExomeDepth)

setwd("/homedirs-porthos/sgul/shares/incc/porthos/Genetics_Centre_Bioinformatics/CNV_calling")

data(exons.hg19)
print(head(exons.hg19))

analysisConfig <- read.csv('config.csv', 
							              header = TRUE, 
							              fill = TRUE)

list_of_bam_files <- as.vector(analysisConfig$list_of_bam_files)

list_of_bam_files

my.counts <- getBamCounts(bed.frame = exons.hg19,
                          bam.files = list_of_bam_files,
                          include.chr = FALSE,)

print(head(my.counts))

# Create dataframe

ExomeCount.dafr <- as(my.counts[, colnames(my.counts)], 'data.frame')

print(head(ExomeCount.dafr))

# Create matrix of the bam counts

ExomeCount.mat <- as.matrix(ExomeCount.dafr[, grep(names(ExomeCount.dafr), 
							pattern = '*.bam')])

print(head(ExomeCount.mat))
							
test <- new('ExomeDepth',
          test = ExomeCount.dafr$S28449_sorted_unique_recalibrated.bam,
          reference = ExomeCount.dafr$S32519_sorted_unique_recalibrated.bam,
          formula = 'cbind(test, reference) ~ 1',
          subset.for.speed = seq(1, nrow(ExomeCount.dafr), 100))

message('Now looping over all the samples innit')

nsamples <- ncol(ExomeCount.mat)

print(head(nsamples))

for (i in 1:nsamples) {


my.test.data <- as.matrix(ExomeCount.mat[, i])

my.reference.set <- as.matrix(ExomeCount.mat[, -i])

head(my.test.data)
head(my.reference.set)

my.choice <- select.reference.set(test.counts = my.test.data,
                                    reference.counts = my.reference.set,
                                    bin.length = (ExomeCount.dafr$end - ExomeCount.dafr$start)/1000,
                                    n.bins.reduced = 10000)


head(my.choice)

my.matrix <- as.matrix( ExomeCount.dafr[, my.choice$reference.choice, drop = FALSE])

my.reference.set <- apply(X = my.matrix,
                                MAR = 1,
                                FUN = sum)


## CNV calling

all.exons <- new('ExomeDepth',
                  test = ExomeCount.dafr$S28449_sorted_unique_recalibrated.bam,
                  reference = my.reference.set,
                  formula = 'cbind(test, reference) ~ 1')

#We can now call the CNV by running the underlying hidden Markov model:

all.exons <- CallCNVs(x = all.exons,
                      transition.probability = 10^-4,
                      chromosome = ExomeCount.dafr$space,
                      start = ExomeCount.dafr$start,
                      end = ExomeCount.dafr$end,
                      name = ExomeCount.dafr$names)

#check output
head(all.exons@CNV.calls)

#Annotating with Conrad Common CNVs

data(Conrad.hg19)

head(Conrad.hg19.common.CNVs)

all.exons <- AnnotateExtra(x = all.exons,
                          reference.annotation = Conrad.hg19.common.CNVs,
                          min.overlap = 0.5,
                          column.name = 'Conrad.hg19')



print(head(all.exons@CNV.calls))

#Now save it in an easily readable format

output.file <- paste('Exome_', 
					          i, 
					          '_', 	
					          list_of_bam_files[i], 
					          '.csv', sep = '')

write.csv(file = output.file, 
		      x = all.exons@CNV.calls, 
		      row.names = FALSE)

}


q()