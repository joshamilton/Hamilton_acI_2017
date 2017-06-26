#%%#############################################################################
# coverage.py
# Copyright (c) 2017, Joshua J Hamilton and Katherine D McMahon
# Affiliation: Department of Bacteriology
#              University of Wisconsin-Madison, Madison, Wisconsin, USA
# URL: http://http://mcmahonlab.wisc.edu/
# All rights reserved.
################################################################################
# This function plots coverage of each contig and genome
################################################################################

#%%#############################################################################
### Import packages
################################################################################

from Bio import SeqIO
import math
import os
import pandas as pd
import subprocess

#%%#############################################################################
### Static folder structure
################################################################################

# Define fixed input and output files
concatFolder = '../../data/refGenomes/concat'
genomeFolder = '../../data/refGenomes/fnaKBase'
gffFolder = '../../data/refGenomes/gff'
sampleFolder = '../../data/sequences'
mapFolder = '../../data/mapping'
bamFolder = '../../data/mapping/bamFiles'
coverageFolder = '../../data/mapping/coverage-pooled'

# Check that the new output directory exists and create if it doesn't

if not os.path.exists(mapFolder):
        os.makedirs(mapFolder)

if not os.path.exists(bamFolder):
        os.makedirs(bamFolder)

if not os.path.exists(coverageFolder):
        os.makedirs(coverageFolder)

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

#%%#############################################################################
### Compute coverage of each individual genome in each sample
### Part 1: Compute SAM to BAM
################################################################################

# Using bamtools, compute the depth of each position along the chromosome
## First convert sam to indexed bam
for sample in sampleList:
    print('Indexing sample '+str(sample))
    for concat in concatList:
        subprocess.call('samtools view -bS '+bamFolder+'/'+sample+'-'+concat+'.sam | samtools sort -o '+bamFolder+'/'+sample+'-'+concat+'.bam', shell=True)

#%%#############################################################################
### Compute coverage of each individual genome in each sample
### Part 2: Compute depth of each base
################################################################################

## Compute depth for each sample
for sample in sampleList:
    print('Depth calculations for sample '+str(sample))
    for concat in concatList:
        subprocess.call('samtools depth '+bamFolder+'/'+sample+'-'+concat+'.bam > '+coverageFolder+'/'+sample+'-'+concat+'.depth', shell=True)

#%%#############################################################################
### Compute coverage of each individual genome
### Part 3: Compute coverage of each contig and genome
################################################################################

# Using FASTA files (contigs in each genome), compute average coverage of each genome and contig

## Create dataframe to store genome coverage on a per-sample basis
coverageDF = pd.DataFrame(index = genomeList, columns=['% Covered', 'Coverage'])

for concat in concatList:
    for genome in genomeList:

        print('Contig coverage for genome '+str(genome))

        # Create a dataframe to store contig depth, length, and coverage
        # Create a dataframe to store the % covered - must track each base on each contig
        contigList = []
        coverContigList = []
        coverPositionList = []
        for curSeq in SeqIO.parse(genomeFolder+'/'+genome+'.fna', 'fasta'):
            contigList.append(curSeq.id)
            coverContigList = coverContigList + ([curSeq.id] * len(curSeq))
            coverPositionList = coverPositionList + list(range(1,len(curSeq)+1))
        contigDF = pd.DataFrame(0, index = contigList, columns=['Covered', 'Depth', 'Length', '% Covered', 'Coverage'])
        coverDF = pd.DataFrame(0, index = pd.MultiIndex.from_tuples(list(zip(*[coverContigList, coverPositionList]))), columns=['Depth'])
        
        for sample in sampleList:
    
            # Make a depth file for the genome
            subprocess.call('grep '+genome+' '+coverageFolder+'/'+sample+'-'+concat+'.depth > '+coverageFolder+'/'+sample+'-'+genome+'.depth', shell=True)
            
            # Import the depth file (if it exists)
            if os.stat(coverageFolder+'/'+sample+'-'+genome+'.depth').st_size > 0:
                depthDF = pd.read_csv(coverageFolder+'/'+sample+'-'+genome+'.depth', index_col = [0,1], names = ['Depth'], sep='\t')
                
                # Store depth of each base across samples
                #coverDF['Depth'] = coverDF['Depth'] + depthDF['Depth']
                coverDF = coverDF.add(depthDF, axis=1, fill_value=0)

                # Compute the depth of each contig
                for curIndex in depthDF.index.levels[0]:
                    contigDF.loc[curIndex, 'Depth'] = contigDF.loc[curIndex, 'Depth'] + depthDF.loc[curIndex].sum()[0]
        
        # Store the length of each contig
        for curSeq in SeqIO.parse(genomeFolder+'/'+genome+'.fna', 'fasta'):
            contigDF.loc[curSeq.id, 'Length'] = len(curSeq)
            
        # Compute the coverage of each contig across all samples
        contigDF['Coverage'] = contigDF['Depth'] / contigDF['Length']
                
        # Compute the % covered of each contig across all samples
        for curIndex in coverDF.index.levels[0]:
            # Subset the coverDF belonging to this contig
            subsetCoverDF = coverDF.loc[curIndex]
            # Subset the coverDF having nonzero depth
            subsetCoverDF = subsetCoverDF.loc[subsetCoverDF['Depth']>0]

            contigDF.loc[curIndex, 'Covered'] = len(subsetCoverDF)
            
            # Check to see if contig is covered. Update coverage appropriately.
            if len(subsetCoverDF) > 0:
                contigDF.loc[curIndex, '% Covered'] = contigDF.loc[curIndex, 'Covered'] / contigDF.loc[curIndex, 'Length']
            else:
                contigDF.loc[curIndex, '% Covered'] = 0

        contigDF.to_csv(coverageFolder+'/'+genome+'.contig.coverage')        

        # Update coverage and % covered for the entire genome
        coverageDF.loc[genome, 'Coverage'] = float(contigDF['Depth'].sum()) / contigDF['Length'].sum()
        coverageDF.loc[genome, '% Covered'] = float(contigDF['Covered'].sum()) / contigDF['Length'].sum()

    coverageDF.to_csv(coverageFolder+'/coverage.csv')

#%%#############################################################################
### Compute coverage of each individual genome
### Part 4: Compute coverage of each individual gene in each genome
################################################################################

# Using GFF files (contigs in each genome), compute average coverage of each gene

for concat in concatList:
    for genome in genomeList:
        print('Gene coverage for genome '+str(genome))

        # Create a dataframe to store gene depth, length, and coverage
        geneDF = pd.DataFrame(columns=['Contig', 'Start', 'Stop', 'Length', 'Covered', 'Depth', '% Covered', 'Coverage', 'Evenness'])
        
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
                geneDF.loc[locus, 'Length'] = stop - start + 1

        # Create a dataframe to store contig depth at each position
        contigList = []
        coverContigList = []
        coverPositionList = []
        for curSeq in SeqIO.parse(genomeFolder+'/'+genome+'.fna', 'fasta'):
            contigList.append(curSeq.id)
            coverContigList = coverContigList + ([curSeq.id] * len(curSeq))
            coverPositionList = coverPositionList + list(range(1,len(curSeq)+1))

        coverDF = pd.DataFrame(0, index = pd.MultiIndex.from_tuples(list(zip(*[coverContigList, coverPositionList]))), columns=['Depth'])

        for sample in sampleList:
    
            # Make a depth file for the genome
            subprocess.call('grep '+genome+' '+coverageFolder+'/'+sample+'-'+concat+'.depth > '+coverageFolder+'/'+sample+'-'+genome+'.depth', shell=True)
            
            # Import the depth file (if it exists)
            if os.stat(coverageFolder+'/'+sample+'-'+genome+'.depth').st_size > 0:
                depthDF = pd.read_csv(coverageFolder+'/'+sample+'-'+genome+'.depth', index_col = [0,1], names = ['Depth'], sep='\t')
                
                # Store depth of each base across samples
                #coverDF['Depth'] = coverDF['Depth'] + depthDF['Depth']
                coverDF = coverDF.add(depthDF, axis=1, fill_value=0)

        # Compute the depth and coverage of each gene
        for curIndex in geneDF.index:
            contig = geneDF.loc[curIndex]['Contig']
            start = geneDF.loc[curIndex]['Start']
            stop = geneDF.loc[curIndex]['Stop']
            
            # If the depth file exists:
            if os.stat(coverageFolder+'/'+sample+'-'+genome+'.depth').st_size > 0:
                # Subset the depthDF belonging to this contig and range
                subsetCoverDF = coverDF.loc[(contig, start):(contig, stop)]
                # Subset the depthDF having nonzero depth
                subsetCoverDF = subsetCoverDF.loc[subsetCoverDF['Depth']>0]
                
                # Update the geneDF with coverage and depth of each gene
                geneDF.loc[curIndex, 'Covered'] = len(subsetCoverDF)
                geneDF.loc[curIndex, 'Depth'] = subsetCoverDF.sum()[0]
    
                # Compute the coverage of each gene
                geneDF.loc[curIndex, '% Covered'] = float(geneDF.loc[curIndex, 'Covered']) / float(geneDF.loc[curIndex, 'Length'])
                geneDF.loc[curIndex, 'Coverage'] = float(geneDF.loc[curIndex, 'Depth']) / float(geneDF.loc[curIndex, 'Length'])
        
                # If the gene is covered...
                # Compute the evenness of coverage using Pielou's eveness
                if len(subsetCoverDF) > 0:
                    countList = subsetCoverDF['Depth'].value_counts().tolist()
                    # Account for unmapped loci
                    if (geneDF.loc[curIndex]['Stop'] - geneDF.loc[curIndex]['Start'] + 1) > sum(countList):
                        countList.append((geneDF.loc[curIndex]['Stop'] - geneDF.loc[curIndex]['Start'] + 1) - sum(countList))
                
                    freqDist = [float(x) / sum(countList) for x in countList]
                    Hprime = 0
                    for freq in freqDist:
                        Hprime = Hprime + freq*math.log(freq)
                    # If only one frequency, assign an evenness of 1
                    if len(freqDist) == 1:
                        geneDF.loc[curIndex, 'Evenness'] = 1
                    else:
                        Hprime_max = math.log(len(freqDist))
                        geneDF.loc[curIndex, 'Evenness'] = - Hprime / Hprime_max

        # Simplify the dataframe and write to file
        geneDF = geneDF.drop(['Contig', 'Start', 'Stop', 'Length', 'Covered', 'Depth'], 1)
        geneDF.to_csv(coverageFolder+'/'+genome+'.gene.coverage')
