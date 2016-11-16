#%%#############################################################################
# filterAndAggregate.py
# Copyright (c) 2016, Joshua J Hamilton and Katherine D McMahon
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
genomeFolder = '../../data/refGenomes/concat'
sampleFolder = '../../data/sequences'
mapFolder = '../../data/mapping'
bamFolder = '../../data/mapping/bamFiles'
countFolder = '../../data/mapping/htseq'

cogTable = '../../data/orthoMCL/cogTable.csv'
taxonFile = '../../data/externalData/taxonomy.csv'

# Check that the new output directory exists and create if it doesn't
if not os.path.exists(countFolder):
        print "Creating output directory\n"
        os.makedirs(countFolder)

#%%#############################################################################
### Read in MT and genome lists. Create DF to store read countDict.
################################################################################
# Read in list of MTs
sampleList = []
for sample in os.listdir(sampleFolder):
    if sample.endswith('.fastq'):
       sampleList.append(sample)

sampleList = [sample.replace('.fastq', '') for sample in sampleList]

# Read in list of genomes. Ignore internal standard genome.
genomeList = []
for genome in os.listdir(genomeFolder):
    if genome.endswith('.fna'):
       genomeList.append(genome)

genomeList = [genome.replace('.fna', '') for genome in genomeList]

#%%#############################################################################
### Count the reads which align to each CDS
################################################################################

# Define parameters for HTSeq-Count script, which we reimplement below.
minQual = 0
featureType = 'CDS'
idAttr = 'locus_tag'
overlapMode = 'intersection-strict'

for sample in sampleList:
    for genome in genomeList:
        samFile = bamFolder+'/'+sample+'-'+genome+'.sam'
        gffFile = genomeFolder+'/'+genome+'.gff' 
        outFile = countFolder+'/'+sample+'-'+genome+'.CDS.out'
        
        subprocess.call('htseq-count -f sam -r pos -s no -a 0 -t CDS -i '+
                        'locus_tag -m intersection-strict '+samFile+' '+
                        gffFile+' > '+outFile, shell=True)
                        
#%%#############################################################################
### Filtering. In this section, filter out all coding sequences which do not 
### have at least one read in all samples.
################################################################################

# First, read in the read counts for each CDS

# Create empty dataframe to merge into
tempDF = pd.read_csv(countFolder+'/'+sampleList[0]+'-'+genomeList[0]+'.CDS.out', sep='\t', index_col=0, names=[sampleList[0]])
readCountDF = pd.DataFrame(index=tempDF.index)

# And read in the counts
for sample in sampleList:
    for genome in genomeList:
        tempDF = pd.read_csv(countFolder+'/'+sample+'-'+genome+'.CDS.out', sep='\t', index_col=0, names=[sample])
        tempDF = tempDF[:-5]
        # Merge with readCountDF
        readCountDF = pd.concat([readCountDF, tempDF], axis=1, join='outer')

# Drop stats from the readCountsDF
readCountDF = readCountDF[:-5]

# Filter the results
cutoff = 1
for sample in sampleList:
    readCountDF = readCountDF.loc[readCountDF[sample] >= cutoff ]

# Write new CDS files
for sample in sampleList:
    tempDF = readCountDF[sample]
    tempDF.to_csv(countFolder+'/'+sample+'-'+genome+'.CDS.out', sep='\t')
    
#%%#############################################################################
### Count total and unique reads which map to each (clade, group) pairing
### Requires integrating data from taxonomy and cog tables into a single data
###  structure
################################################################################

# Read in the taxonomy table and create a list of genomes for each clade
cladeToGenomeDict = {}
cladeList = []

for genome in genomeList:
    taxonClass = pd.DataFrame.from_csv(taxonFile, sep=',')
    taxonClass = taxonClass.dropna()
    
# Extract the unique clades
    taxonClass = taxonClass.drop(taxonClass[taxonClass['Lineage'] != genome].index)
    innerCladeList = pd.unique(taxonClass['Clade'].values)
    
    for clade in innerCladeList:
        innerGenomeList = taxonClass[taxonClass['Clade'] == clade].index.tolist()
        cladeToGenomeDict[clade] = innerGenomeList

    cladeList = cladeList + innerCladeList.tolist()

# Read in the COG table
cogTableDF = pd.read_csv(cogTable, index_col=0)

# Create and populate the dataframe, indexed by clade and group
cladeCogToCdsIndex = pd.MultiIndex.from_product([cladeList, cogTableDF.index.tolist()], names=['Clade', 'COG'])
cladeCogToCdsDF = pd.DataFrame(index=cladeCogToCdsIndex, columns=['CDS'])

for index in cladeCogToCdsDF.index:
    clade = index[0]
    cog = index[1]
    innerGenomeList = cladeToGenomeDict[clade]
    cdsList = []

    for innerGenome in innerGenomeList:
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

#%%#############################################################################
### Count total reads which map to each (clade, group) pairing
### Requires integrating data from taxonomy and cog tables into a single data
###  structure
################################################################################

# Convert the cladeCogToCdsDF to a multi-index so we can construct the 
# count table for COGs

cladeCogToCdsDF = cladeCogToCdsDF.set_index([cladeCogToCdsDF.index, 'COG'])

# For each MT, construct the count table and write to file
for sample in sampleList:
    for genome in genomeList:    
    # Reset columns for total and unique reads
        cladeCogToCdsDF['Total'] = 0

        for index in cladeCogToCdsDF.index:
            cdsList = cladeCogToCdsDF.loc[index, 'CDS'].split(',')
            readList = []
            for cds in cdsList:
                if cds in readCountDF.index:
                    cladeCogToCdsDF = cladeCogToCdsDF.set_value(index, 'Total', cladeCogToCdsDF.loc[index, 'Total'] + readCountDF.loc[cds, sample])
        
        cladeCogToCdsDF = cladeCogToCdsDF.loc[cladeCogToCdsDF['Total'] > 0 ]
        cladeCogToCdsDF.to_csv(countFolder+'/'+sample+'-'+genome+'.COG.out')
