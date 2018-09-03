#%%#############################################################################
# filterAndAggregate.py
# Copyright (c) 2017, Joshua J Hamilton and Katherine D McMahon
# Affiliation: Department of Bacteriology
#              University of Wisconsin-Madison, Madison, Wisconsin, USA
# URL: http://http://mcmahonlab.wisc.edu/
# All rights reserved.
################################################################################
# Count reads which map to each (genome, gene) pair and to each (clade, COG)
# pair.
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
concatFolder = '../../data/refGenomes/concat'
genomeFolder = '../../data/refGenomes/fna'
sampleFolder = '../../data/sequences'
mapFolder = '../../data/mapping'
bamFolder = '../../data/mapping/bamFiles'
coverageFolder = '../../data/mapping/coverage-pooled'
countFolder = '../../data/mapping/htseq'

cogTable = '../../data/orthoMCL/cogTable.csv'
taxonFile = '../../data/externalData/taxonomy.csv'

# Check that the new output directory exists and create if it doesn't
if not os.path.exists(countFolder):
        print("Creating output directory\n")
        os.makedirs(countFolder)

#%%#############################################################################
### Read in sample and genome lists. Create DF to store read countDict.
################################################################################
# Read in list of samples
sampleList = []
for sample in os.listdir(sampleFolder):
    if sample.endswith('.fastq'):
       sampleList.append(sample)
sampleList = [sample.replace('.fastq', '') for sample in sampleList]

# Read in list of genomes.
genomeList = []
for genome in os.listdir(genomeFolder):
    if genome.endswith('.fna'):
       genomeList.append(genome)

genomeList = [genome.replace('.fna', '') for genome in genomeList]

# Read in list of genomes.
concatList = []
for concat in os.listdir(concatFolder):
    if concat.endswith('.fna'):
       concatList.append(concat)

concatList = [concat.replace('.fna', '') for concat in concatList]

#%%#############################################################################
### Count the reads which align to each CDS
################################################################################

# Define parameters for HTSeq-Count script
minQual = 0
featureType = 'CDS'
idAttr = 'locus_tag'
overlapMode = 'intersection-strict'

for sample in sampleList:
    for concat in concatList:
        samFile = bamFolder+'/'+sample+'-'+concat+'.sam'
        gffFile = concatFolder+'/'+concat+'.gff' 
        outFile = countFolder+'/'+sample+'-'+concat+'.CDS.out'
        
        subprocess.call('htseq-count -f sam -r pos -s no -a 0 -t CDS -i '+
                        'locus_tag -m intersection-strict '+samFile+' '+
                        gffFile+' > '+outFile, shell=True)
                        
#%%#############################################################################
### Filtering. In this section, filter out all coding sequences which do not
### recruit at least 50 reads
################################################################################

# First, read in the read counts for each CDS

# Create empty dataframe to merge into
tempDF = pd.read_csv(countFolder+'/'+sampleList[0]+'-'+concatList[0]+'.CDS.out', sep='\t', index_col=0, names=[sampleList[0]])
readCountDF = pd.DataFrame(index=tempDF.index)

# And read in the counts
for sample in sampleList:
    for concat in concatList:
        tempDF = pd.read_csv(countFolder+'/'+sample+'-'+concat+'.CDS.out', sep='\t', index_col=0, names=[sample])
        tempDF = tempDF[:-5]
        # Merge with readCountDF
        readCountDF = pd.concat([readCountDF, tempDF], axis=1, join='outer')

## Drop stats from the readCountsDF
readCountDF = readCountDF.drop(['__alignment_not_unique', '__ambiguous', '__no_feature', '__not_aligned', '__too_low_aQual'], axis=0)
readCountDF = readCountDF.sum(axis=1)

readCountDF.to_csv(countFolder+'/readCounts.csv', sep=',')


# Filter the results by dropping all genes which don't recruit at least ten reads
readCutoff = 0
readCountDF = readCountDF.loc[readCountDF >= readCutoff]
readCountDF.to_csv(countFolder+'/filteredReadCounts.csv', sep=',')
    
#%%#############################################################################
### Integrate data from taxonomy and cog tables into a single data
###  structure
################################################################################

# Read in the taxonomy table and create a list of genomes for each clade
cladeToGenomeDict = {}
cladeList = []

for concat in concatList:
    taxonClass = pd.DataFrame.from_csv(taxonFile, sep=',')
    taxonClass = taxonClass.dropna()
    
# Extract the unique clades
    taxonClass = taxonClass.drop(taxonClass[taxonClass['Lineage'] != concat].index)
    innerCladeList = pd.unique(taxonClass['Clade'].values)
    
    for clade in innerCladeList:
        innerconcatList = taxonClass[taxonClass['Clade'] == clade].index.tolist()
        cladeToGenomeDict[clade] = innerconcatList

    cladeList = cladeList + innerCladeList.tolist()

# Read in the COG table
cogTableDF = pd.read_csv(cogTable, index_col=0)

# Create and populate the dataframe, indexed by clade and group
cladeCogToCdsIndex = pd.MultiIndex.from_product([cladeList, cogTableDF.index.tolist()], names=['Clade', 'COG'])
cladeCogToCdsDF = pd.DataFrame(index=cladeCogToCdsIndex, columns=['CDS'])

for index in cladeCogToCdsDF.index:
    clade = index[0]
    cog = index[1]
    innerconcatList = cladeToGenomeDict[clade]
    cdsList = []

    for innerGenome in innerconcatList:
        if not pd.isnull(cogTableDF.loc[cog][innerGenome]):
            tempList = cogTableDF.loc[cog][innerGenome].split(';')        
            cdsList = cdsList + tempList

    cladeCogToCdsDF.loc[index] = ','.join(cdsList)

# Sort by the multi-index and write to file
cladeCogToCdsDF = cladeCogToCdsDF.sort_index(axis=0)
cladeCogToCdsDF.to_csv(mapFolder+'/cladesCogsToCDS.csv')

# Read in singly-indexed, drop empty rows, and write to file
cladeCogToCdsDF = pd.read_csv(mapFolder+'/cladesCogsToCDS.csv', index_col=0)
cladeCogToCdsDF = cladeCogToCdsDF.dropna(axis=0, how='any')

cladeCogToCdsDF.to_csv(mapFolder+'/cladesCogsToCDS.csv')
