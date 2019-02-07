### ExomeDepth - Alan Pittman - Lucas Cortes - reusable_ Rscript - CNV caller April 2018 

#requires index .bam files to run
#Use with the congifguration file (list number of samples basically)
#calls CNVs on sample set(~4-10) using rest as reference files


library("ExomeDepth")
library("GenomeInfoDb")
library("Rsamtools")

data(exons.hg19)

################

#example config file with header (save as a config.csv file) 
################################################
#requires config file with header listing the .bam files you want to analyse
#save as a config.csv file 

#list_of_bam_files
#1_realigned.bam
#2_realigned.bam
#3_realigned.bam
#3_realigned.bam

message('reading config file')

analysisConfig <- read.csv('config.csv', 
							header = TRUE, 
							fill = TRUE)

list_of_bam_files <- as.vector(analysisConfig$list_of_bam_files)

print(head(list_of_bam_files))

#generate counts of out list of bam files
my.counts <- getBamCounts(bed.frame = exons.hg19, 
						bam.files = list_of_bam_files, 
						include.chr = FALSE)

#check the counting worked:
print(head(my.counts))

# Create dataframe
ExomeCount.dafr <- as(my.counts[, colnames(my.counts)], 'data.frame')
ExomeCount.dafr$chromosome <- gsub(as.character(ExomeCount.dafr$space), 
									pattern = 'chr',replacement = '')

print(head(ExomeCount.dafr))

# Create matrix of the bam counts
ExomeCount.mat <- as.matrix(ExomeCount.dafr[, grep(names(ExomeCount.dafr), 
							pattern = '*.bam')])

print(head(ExomeCount.mat))

nsamples <- ncol(ExomeCount.mat)

print(head(nsamples))

#now loop over each sample

message('Now looping over all the samples innit')

for (i in 1:nsamples) {

### Create the aggregate reference set for the test samples:

my.choice <- select.reference.set (test.counts =  ExomeCount.mat[,i], reference.counts = ExomeCount.mat[,-i], 
									bin.length = (ExomeCount.dafr$end - ExomeCount.dafr$start)/1000, 
									n.bins.reduced = 10000)

my.reference.selected <- apply(X = ExomeCount.mat[, my.choice$reference.choice, drop = FALSE], 
							MAR = 1, FUN = sum)

message('Now creating the ExomeDepth object')
all.exons <- new('ExomeDepth',
					test = ExomeCount.mat[,i], 
					reference = my.reference.selected, 
					formula = 'cbind(test, reference) ~ 1')			

############### Now call the CNVs
all.exons <- CallCNVs(x = all.exons,
                      transition.probability = 10^-4,
                      chromosome = ExomeCount.dafr$chromosome,
                      start = ExomeCount.dafr$start,
                      end = ExomeCount.dafr$end,
                      name = ExomeCount.dafr$names)				   
							   
#Lets add the annotation of common CNV dataset:

data(Conrad.hg19)

#Then one can use this information to annotate our CNV calls with the function AnnotateExtra from GenomicRanges
#The CNVs must overlap by at least 50% to get annotated.

all.exons <- AnnotateExtra(x = all.exons, 
							reference.annotation = Conrad.hg19.common.CNVs, 
							min.overlap = 0.5, 
							column.name = 'Conrad.hg19')
 
#now annotating with exon/gene level information. 

data(exons.hg19)
 
exons.hg19.GRanges <- GenomicRanges::GRanges(seqnames = exons.hg19$chromosome, 
											IRanges::IRanges(start=exons.hg19$start,end=exons.hg19$end), 
											names = exons.hg19$name)

#here the minimum overlap should be very close to 0  
all.exons <- AnnotateExtra(x = all.exons, 
						reference.annotation = exons.hg19.GRanges, 
						min.overlap = 0.0001, 
						column.name = 'exons.hg19')
 
#writing results to .csv format:

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
