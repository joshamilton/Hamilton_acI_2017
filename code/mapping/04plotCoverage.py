#%%#############################################################################
# plotCoverage.py
# Copyright (c) 2017, Joshua J Hamilton and Katherine D McMahon
# Affiliation: Department of Bacteriology
#              University of Wisconsin-Madison, Madison, Wisconsin, USA
# URL: http://http://mcmahonlab.wisc.edu/
# All rights reserved.
################################################################################
# This function makes coverage plots of each gene for inclusion to ensure
# sufficient depth/evenness for inclusion in RPKM calculations
################################################################################

#%%#############################################################################
### Import packages
################################################################################

from Bio import SeqIO
import matplotlib.pyplot as plt
import os
import pandas as pd
import subprocess

#%%#############################################################################
### Static folder structure
################################################################################

# Define fixed input and output files
concatFolder = '../../data/refGenomes/concat'
externalDataFolder = '../../data/externalData'
genomeFolder = '../../data/refGenomes/fnaKBase'
gffFolder = '../../data/refGenomes/gff'
sampleFolder = '../../data/sequences'
mapFolder = '../../data/mapping'
bamFolder = '../../data/mapping/bamFiles'
coverageFolder = '../../data/mapping/coverage-pooled'
plotFolder = '../../data/mapping/plots'

# Check that the new output directory exists and create if it doesn't
if not os.path.exists(plotFolder):
        os.makedirs(plotFolder)
        
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

concatList = []
for concat in os.listdir(concatFolder):
    if concat.endswith('.fna'):
       concatList.append(concat)

concatList = [concat.replace('.fna', '') for concat in concatList]

##%%#############################################################################
#### For each genome, plot the depth of each base along each contig
#### Part 1 - Compute depth pooled across samples
#################################################################################
#
## Create dataframe to store genome coverage on a per-sample basis
#coverageDF = pd.DataFrame(index = genomeList, columns=['% Covered', 'Coverage'])
#
#for concat in concatList:
#    
#    for genome in genomeList:
#
#        print('Computing pooled depth profile for genome '+genome)
#
#        # Create a dataframe to store the % covered - must track each base on each contig
#        coverContigList = []
#        coverPositionList = []
#        for curSeq in SeqIO.parse(genomeFolder+'/'+genome+'.fna', 'fasta'):
#            coverContigList = coverContigList + ([curSeq.id] * len(curSeq))
#            coverPositionList = coverPositionList + list(range(1,len(curSeq)+1))
#
#        coverDF = pd.DataFrame(0, index = pd.MultiIndex.from_tuples(list(zip(*[coverContigList, coverPositionList]))), columns=['Depth'])
#        
#        for sample in sampleList:
#    
#            # Make a depth file for the genome
#            subprocess.call('grep '+genome+' '+coverageFolder+'/'+sample+'-'+concat+'.depth > '+coverageFolder+'/'+sample+'-'+genome+'.depth', shell=True)
#            
#            # Import the depth file (if it exists)
#            if os.stat(coverageFolder+'/'+sample+'-'+genome+'.depth').st_size > 0:
#                depthDF = pd.read_csv(coverageFolder+'/'+sample+'-'+genome+'.depth', index_col = [0,1], names = ['Depth'], sep='\t')
#                
#                # Store depth of each base across samples
#                coverDF = coverDF.add(depthDF, axis=1, fill_value=0)
#        
#        coverDF.to_csv(coverageFolder+'/'+genome+'.depth')
#        
##%%#############################################################################
#### For each genome, plot the depth of each base along each contig
#### Part 2 - Construct plots
#################################################################################
#
#for genome in genomeList:
#    
#    # Read in depth profile
#    sampleDF = pd.read_csv(coverageFolder+'/'+genome+'.depth', index_col = [0,1], sep=',')
#    
#    # For each contig, plot depth at each position
#    for contig in sampleDF.index.levels[0]:
#        
#        contigsampleDf = sampleDF.loc[contig]
#        contigsampleDf.index.name = 'Position'
#        contigsampleDf = contigsampleDf.reset_index()
#        axes = contigsampleDf.plot(x='Position', y='Depth', kind='line', legend=False, title='sample: '+contig)
#        figure = axes.get_figure()
#        figure.savefig(plotFolder+'/'+contig+'.png')
#        plt.close('all')
        
#%%#############################################################################
### For each genome, plot the depth of each base along each gene
### Part 1 - Construct plots
################################################################################

geneDict = {}
with open(externalDataFolder+'/transporters.csv') as inFile:
    for line in inFile.readlines():   
        line = line.strip()
        key = line.split('.')[0]
        if key in geneDict.keys():
            geneDict[key].append(line)
        else:
            geneDict[key] = [line]

for genome in geneDict.keys():
#for genome in genomeList:
    
    # Read in depth profile
    sampleDF = pd.read_csv(coverageFolder+'/'+genome+'.depth', index_col = [0,1], sep=',')
    
    # Create a dataframe to store gene depth, length, and coverage
    geneDF = pd.DataFrame(columns=['Contig', 'Start', 'Stop'])
    
    # Import the GFF file and use its values to populate the geneDF
    with open(gffFolder+'/'+genome+'.gff', 'r') as gffFile:
        next(gffFile)
        for line in gffFile:
            gffArray = line.split('\t')
            contig = gffArray[0]
            start = int(gffArray[3])
            stop = int(gffArray[4])
            locusArray = gffArray[8].split(';')
            locus = locusArray[0].split('=')[1]

            # Populate the geneDF dataframe using these values
            geneDF.loc[locus, 'Contig'] = contig
            geneDF.loc[locus, 'Start'] = start
            geneDF.loc[locus, 'Stop'] = stop
            
    # For each gene, extract the relevant depth profile
    for gene in geneDict[genome]:
#    for gene in geneDF.index:
        contig = geneDF.loc[gene, 'Contig']
        start = geneDF.loc[gene, 'Start']
        stop = geneDF.loc[gene, 'Stop']
        
        genesampleDf = sampleDF.loc[(contig, start):(contig, stop)]
        genesampleDf = genesampleDf.loc[contig]
        genesampleDf.index.name = 'Position'
        genesampleDf = genesampleDf.reset_index()
        axes = genesampleDf.plot(x='Position', y='Depth', kind='line', legend=False, title='sample: '+gene)
        axes.set_ylim(ymin = 0)
        figure = axes.get_figure()
        figure.savefig(plotFolder+'/'+gene+'.png')
        plt.close('all')