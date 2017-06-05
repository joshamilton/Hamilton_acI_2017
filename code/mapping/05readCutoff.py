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

from Bio import SeqIO
import os
import pandas as pd
import subprocess

#%%#############################################################################
### Static folder structure
################################################################################
# Define fixed input and output files
concatFolder = '../../data/refGenomes/concat'
genomeFolder = '../../data/refGenomes/fna'
ffnFolder = '../../data/refGenomes/ffn'
sampleFolder = '../../data/sequences'
mapFolder = '../../data/mapping'
bamFolder = '../../data/mapping/bamFiles'
coverageFolder = '../data/mapping/coverage-pooled'
countFolder = '../../data/mapping/cutoff'

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

# Read in list of concatenated genomes
concatList = []
for genome in os.listdir(concatFolder):
    if genome.endswith('.fna'):
       concatList.append(genome)

concatList = [genome.replace('.fna', '') for genome in concatList]

# Read in list of genomes.
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
    for concat in concatList:
        samFile = bamFolder+'/'+sample+'-'+concat+'.sam'
        gffFile = concatFolder+'/'+concat+'.gff' 
        outFile = countFolder+'/'+sample+'-'+concat+'.CDS.out'
        
        subprocess.call('htseq-count -f sam -r pos -s no -a 0 -t CDS -i '+
                        'locus_tag -m intersection-strict '+samFile+' '+
                        gffFile+' > '+outFile, shell=True)
                        
#%%#############################################################################
### Aggregate read counts
#%%#############################################################################

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

# Drop stats from the readCountsDF
readCountDF = readCountDF[:-5]

# Compute total reads across all four samples
readCountDF = readCountDF.sum(axis=1)
readCountDF.to_csv(countFolder+'/readCounts.csv')

#%%#############################################################################f
### Compute RPKM for each gene
#%%#############################################################################

# Read in the dataframe of readCounts
readCountDF = pd.read_csv(countFolder+'/readCounts.csv', names = ['Counts'], index_col=0)

for genome in genomeList:

    # Subset to include the genome of interest
    rpkmCountDF = readCountDF[readCountDF.index.str.contains(genome)]
    
    # Compute total mapped reads
    totalReads = rpkmCountDF.sum()[0]
    
    # Create an empty gene length dictionary. Then read in the individual fasta 
    # files, and store the length of each gene in the dictionary. Write to file.
    geneLengthDict = {}
    for seqRecord in SeqIO.parse(ffnFolder+'/'+genome+'.ffn', "fasta"):
        geneLengthDict[seqRecord.description.split(' ')[0]] = float(len(seqRecord.seq))

    # Expand the readCounts DF with an RPKM column
    for gene in rpkmCountDF.index:
        readCountDF.loc[gene, 'RPKM'] = readCountDF.loc[gene, 'Counts'] / ((geneLengthDict[gene] / 1000) * (totalReads / 1000000))
    
readCountDF.to_csv(countFolder+'/rpkmCounts.csv')
    
#%%#############################################################################
### Visualize RPKM values for each unique read count as a boxplot
#%%#############################################################################

rpkmCountDF = pd.read_csv(countFolder+'/rpkmCounts.csv', index_col=0)

axes = rpkmCountDF.boxplot(column='RPKM', by='Counts')
figure = axes.get_figure()
figure.savefig(countFolder+'/readCountDF-boxplot.png')

# This plot shows there is a lot of noise at low count values, as evidenced
# by the numnber of outliers.

# Let's calculate the number of outliers for each count value
outliersDF = pd.DataFrame(0, index=rpkmCountDF['Counts'].unique(), columns=['Outliers'])

for count in rpkmCountDF['Counts'].unique():
    _, bp = pd.DataFrame.boxplot(rpkmCountDF[rpkmCountDF['Counts'] == count], return_type='both')

    outliers = [flier.get_ydata() for flier in bp["fliers"]]
    outliersDF.loc[count, 'Outliers'] = len(outliers[1])
    
# Plot the results
outliersDF.index.name = 'Count'
outliersDF = outliersDF.reset_index()
axes = outliersDF.plot(x='Count', y='Outliers', kind='scatter', legend=False)
axes.set_xlim(0, 100)
axes.set_ylim(0, 100)
figure = axes.get_figure()
figure.savefig(countFolder+'/outliers.png')

# As can be seen from this plot, the number of outliers levels off after 
# around 45 reads. Therefore, I will use a cutoff of 50 reads when calculating
# RPKMs.
