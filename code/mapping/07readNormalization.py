#%%#############################################################################
# readNormalization.py
# Copyright (c) 2017, Joshua J Hamilton and Katherine D McMahon
# Affiliation: Department of Bacteriology
#              University of Wisconsin-Madison, Madison, Wisconsin, USA
# URL: http://http://mcmahonlab.wisc.edu/
# All rights reserved.
################################################################################
# Normalize read count data. Given read counts for (genome, gene) pairings, 
# compute RPKM-normalized counts for (clade, COG) pairs.
################################################################################

#%%#############################################################################
### Import packages
################################################################################

from Bio import SeqIO
from collections import Counter
import os
import pandas as pd
import re
import subprocess

#%%#############################################################################
### Static folder structure
################################################################################
# Define fixed input and output files
genomeFolder = '../../data/refGenomes/concat'
ffnFolder = '../../data/refGenomes/ffn'
sampleFolder = '../../data/sequences'
mapFolder = '../../data/mapping'
bamFolder = '../../data/mapping/bamFiles'
countFolder = '../../data/mapping/htseq'
outputFolder = '../../results/expression'

cogTable = '../../data/orthoMCL/cogTable.csv'
taxonFile = '../../data/externalData/taxonomy.csv'

cladesCogsToCDSTable = mapFolder+'/cladesCogsToCDS.csv'
annotTable = '../../data/orthoMCL/annotTable.csv'

# Check that the new output directory exists and create if it doesn't
if not os.path.exists(outputFolder):
        os.makedirs(outputFolder)

#%%#############################################################################
### Read in MT and genome lists. Create DF to store read countDict.
################################################################################
# Read in list of MTs
sampleList = []
for sample in os.listdir(sampleFolder):
    if sample.endswith('.fastq'):
       sampleList.append(sample)

sampleList = [sample.replace('.fastq', '') for sample in sampleList]

# Read in list of individual genomes
genomeList = []
for genome in os.listdir(ffnFolder):
    if genome.endswith('.ffn'):
       genomeList.append(genome)
       
genomeList = [genome.replace('.ffn', '') for genome in genomeList]
       
# Read in list of concatenated genomes
concatList = []
for genome in os.listdir(genomeFolder):
    if genome.endswith('.fna'):
       concatList.append(genome)

concatList = [genome.replace('.fna', '') for genome in concatList]

# Use the taxonomy file to obtain the correct genomes
taxonClass = pd.DataFrame.from_csv(taxonFile, sep=',')
taxonClass = taxonClass.dropna()
        
#%%#############################################################################
### Read in dataframe showing CDS associated with each (clade, COG)
################################################################################

cladeCogToCdsDF = pd.read_csv(cladesCogsToCDSTable, index_col=[0, 1])

#%%#############################################################################
### Establish data structures necessary for normalization:
### totalMappedReads, giving genome and total number of mapped reads (per
###  million basis), for 'M' in RPKM
### geneLength, giving gene length of each gene (per kilobase basis), for 'K'
###  in RPKM
################################################################################

totalMappedReadsDF = pd.DataFrame(index=genomeList,columns=sampleList)
for concat in concatList:
    for sample in sampleList:
        for genome in genomeList:
            totalMappedReadsDF.loc[genome][sample] = float(subprocess.check_output('samtools view -F 0x4 '+bamFolder+'/'+sample+'-'+genome+'.bam | cut -f 1 | sort | uniq | wc -l', shell=True).strip())

totalMappedReadsDF.to_csv(mapFolder+'/'+'genomeMappedReads.csv')

totalMappedReadsDF = totalMappedReadsDF.sum(axis=1)

# Create an empty gene length dictionary. Then read in the individual fasta 
# files, and store the length of each gene in the dictionary. Write to file.
geneLengthDict = {}
for genome in genomeList:
    for seqRecord in SeqIO.parse(ffnFolder+'/'+genome+'.ffn', "fasta"):
        geneLengthDict[seqRecord.description.split(' ')[0]] = float(len(seqRecord.seq))
        
#%%#############################################################################
### Count total and unique reads which map to each (clade, group) pairing
################################################################################

for concat in concatList:    

    # Create a dataframe to store results
    cladeCogNormDF = cladeCogToCdsDF.copy()
    cladeCogNormDF = cladeCogNormDF.drop('CDS', axis=1)
    cladeCogNormDF['RPKM'] = float(0)

    # Read in the DF of gene counts
    geneCountDF = pd.read_csv(countFolder+'/filteredReadCounts.csv', sep=',', index_col=0, header=None, names=['Count'])
    # For each (clade, COG) pairing, read in the list of genes
    for index in cladeCogToCdsDF.index:
        totRPKM = float(0)
        if not pd.isnull(cladeCogToCdsDF.loc[index, 'CDS']):
            cdsList = cladeCogToCdsDF.loc[index, 'CDS'].split(',')
            # For each gene, comput the RPKM and update for the COG
            for cds in cdsList:
                if cds in geneCountDF.index:
                    genome = cds.split('.')[0]
                    # Sometimes low abundance genomes don't map any reads, so check for this before computing RPKM
                    if totalMappedReadsDF.loc[genome] > 0:
                        curRPKM = float(geneCountDF.loc[cds]['Count']) / ((geneLengthDict[cds] / 1000)*(totalMappedReadsDF[genome] / 1000000))
                        totRPKM = totRPKM + curRPKM
        cladeCogNormDF.loc[index]['RPKM'] = totRPKM

    cladeCogNormDF = cladeCogNormDF.loc[cladeCogNormDF['RPKM'] > 0]
    cladeCogNormDF.to_csv(countFolder+'/'+concat+'.COG.norm')

#%%#############################################################################
### Now average the RPKMs for the indicated sets of samples and normalize
### Also compute the percentile rank for each clade, and extract the majority
### annotation
################################################################################

# Create additional columns
rpkmDF = cladeCogNormDF.copy()
rpkmDF['Median'] = 0
rpkmDF['Rel to Med'] = 0
rpkmDF['Annotation'] = 'None'

# Compute the Log 2 RPKM and replace infinity with zero
#rpkmDF['Log2 RPKM'] = np.log2(rpkmDF['RPKM'])
#rpkmDF = rpkmDF.replace(to_replace='-inf', value='0')

# Compute the median RPKM
for clade in rpkmDF.index.levels[0]:
    median = rpkmDF.loc[clade, 'RPKM'].median()
    rpkmDF.loc[(clade), 'Median'] = median

rpkmDF['Rel to Med'] = rpkmDF['RPKM'] / rpkmDF['Median']

# Read in the annotation table
annotDF = pd.read_csv(annotTable, index_col=0)

# Split across samples and compute the majority annotation and percentile rank for each COG
splitDF = rpkmDF.groupby(level=0)
for (clade, COGs) in splitDF:

    # Compute majority annotation
    # First subset the dataframe, keep genomes for that clade and dropping 
    genomeList = taxonClass.loc[taxonClass['Clade'] == clade].index.tolist()
    tempDF = annotDF[genomeList].dropna(axis=0, how='all')
    
    for cog in tempDF.index:
        annotList = []
        for genome in tempDF.columns:
            if not pd.isnull(tempDF.loc[cog][genome]):
                innerString = tempDF.loc[cog][genome]
                # Dataframe element is a string enclosed in brackets with a comma separating elements
                innerString = re.sub('[\[\]]' , '', innerString)
                innerList = re.split('\', \'|\", \"', innerString)
                innerList = [re.sub('\"|\'', '', string) for string in innerList]
                annotList = annotList + innerList
        # Find the most common 
        annotCounter = Counter(annotList)
        majorityAnnot = annotCounter.most_common(1)[0][0]
            
        # Assign the Annotation
        COGs.set_value((clade, cog), 'Annotation', majorityAnnot)
        
    # Drop empty rows
    COGs = COGs.dropna(axis=0, how='any')
    
    COGs.to_csv(outputFolder+'/'+clade+'.norm')
