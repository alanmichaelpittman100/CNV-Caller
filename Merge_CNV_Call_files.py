#!/usr/bin/env python

#Alan Pittman Julu 2018
#generates a single .xlsx sheet containg all the CNV data labelled by sample

import sys
import glob


orig_stdout = sys.stdout
f = open('All_cnv_calls.csv', 'w')
sys.stdout = f

for file in glob.glob('*bam.csv'):

    for line in open(file):
	    print(file + "," + line,)
		
sys.stdout = orig_stdout
f.close()			


# remove first line..
# leftstrip ('
# rightstip \n',)