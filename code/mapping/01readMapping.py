#%%#############################################################################
# readMapping.py
# Copyright (c) 2017, Joshua J Hamilton and Katherine D McMahon
# Affiliation: Department of Bacteriology
#              University of Wisconsin-Madison, Madison, Wisconsin, USA
# URL: http://http://mcmahonlab.wisc.edu/
# All rights reserved.
################################################################################
# Use BBMap to map metranscriptomic reads to our reference genome collection
################################################################################

#%%#############################################################################
### Import packages
################################################################################

import os
import pandas as pd
import subprocess

#%%#############################################################################
### Static folder structure
################################################################################

# Define fixed input and output files
genomeFolder = '../../data/refGenomes/concat'
sampleFolder = '../../data/sequences'
mapFolder = '../../data/mapping'
bamFolder = '../../data/mapping/bamFiles'

# Check that the new output directory exists and create if it doesn't

if not os.path.exists(mapFolder):
        os.makedirs(mapFolder)

if not os.path.exists(bamFolder):
        os.makedirs(bamFolder)

#%%#############################################################################
### Read in sample and genome lists.
################################################################################

sampleList = []
for sample in os.listdir(sampleFolder):
    if sample.endswith('.fastq'):
       sampleList.append(sample)
sampleList = [sample.replace('.fastq', '') for sample in sampleList]

genomeList = []
for genome in os.listdir(genomeFolder):
    if genome.endswith('.fna'):
       genomeList.append(genome)

genomeList = [genome.replace('.fna', '') for genome in genomeList]

##%%#############################################################################
#### Run BBMap
#################################################################################
#    
#i = 1
#for sample in sampleList:
#    for genome in genomeList:
#        print('Mapping pair '+str(i)+' of '+str(len(sampleList)*len(genomeList)))
#        
#        subprocess.call(['bbmap.sh', 'ref='+genomeFolder+'/'+genome+'.fna',
#                         'in='+sampleFolder+'/'+sample+'.fastq',
#                         'outm='+bamFolder+'/'+sample+'-'+genome+'.sam',
#                         'minid=0.80', 
#                         'ambig=random', 
#                         'nodisk', 
#                         'sam=1.3'
#                ])
#        i = i+1
        
#%%#############################################################################
### Count reads
################################################################################

readCountDF = pd.DataFrame(index=genomeList, columns=sampleList)

i = 1
for sample in sampleList:
    for genome in genomeList:
        print('Counting pair '+str(i)+' of '+str(len(sampleList)*len(genomeList)))

        readCountDF[sample][genome] = subprocess.check_output(
            'samtools view -F 0x4 '+bamFolder+'/'+sample+'-'+genome+'.sam'+'| cut -f 1 | sort | uniq | wc -l', shell=True).strip().decode('ascii')
        i = i+1
        
readCountDF.to_csv(mapFolder+'/concatMappedReads.csv')