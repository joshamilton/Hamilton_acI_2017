#%%#############################################################################
# aggregateUniqueReads.py
# Copyright (c) 2016, Joshua J Hamilton and Katherine D McMahon
# Affiliation: Department of Bacteriology
#              University of Wisconsin-Madison, Madison, Wisconsin, USA
# URL: http://http://mcmahonlab.wisc.edu/
# All rights reserved.
################################################################################
# Process read count data from mapping of OMD-TOIL MT reads to our reference
# genomes.
################################################################################

#%%#############################################################################
### Import packages
################################################################################

import HTSeq
import os
import pandas as pd

#%%#############################################################################
### Static folder structure
################################################################################
# Define fixed input and output files
genomeFolder = '../../data/refGenomes/concat'
sampleFolder = '../../data/metatranscriptomes'
mapFolder = '../../data/mapping/bamFiles'
countFolder = '../../data/mapping/htseq'
readsDir = '../../data/mapping/reads'

cogTable = '../../data/orthoMCL/cogTable.csv'
taxonFile = '../../data/externalData/taxonomy.csv'

# Check that the new output directory exists and create if it doesn't
if not os.path.exists(countFolder):
        print "Creating output directory\n"
        os.makedirs(countFolder)

#if os.path.exists(readsDir):
#    shutil.rmtree(readsDir)
        
if not os.path.exists(readsDir):
        print "Creating output directory\n"
        os.makedirs(readsDir) 

#%%#############################################################################
### Read in MT and genome lists. Create DF to store read countDict.
################################################################################
# Read in list of MTs
mtList = []
for mt in os.listdir(sampleFolder):
    if mt.endswith('.fastq'):
       mtList.append(mt)

mtList = [mt.replace('.fastq', '') for mt in mtList]

# Read in list of genomes. Ignore internal standard genome.
genomeList = []
for genome in os.listdir(genomeFolder):
    if genome.endswith('.fna'):
       genomeList.append(genome)

genomeList = [genome.replace('.fna', '') for genome in genomeList]

#%%#############################################################################
### Obtain the same information for each (genome, COG) pairing
################################################################################

# Define parameters for HTSeq-Count script, which we reimplement below.
minQual = 0
featureType = 'CDS'
idAttr = 'locus_tag'
overlapMode = 'intersection-strict'

for mt in mtList:
    for genome in genomeList:
        samFile = mapFolder+'/'+mt+'-'+genome+'.sam'
        gffFile = genomeFolder+'/'+genome+'.gff'    

        # CReate feature array to store GFF data        
        features = HTSeq.GenomicArrayOfSets("auto", stranded = False)     
        
        # Create dictionary to store count data
        countDict = {}

        # Read in the GFF file
        gffObj = HTSeq.GFF_Reader(gffFile)   
        i = 0
        for gffLine in gffObj:
            if gffLine.type == featureType:
                featureID = gffLine.attr[idAttr]
                features[gffLine.iv] += featureID
            # Create dictionary entry
                countDict[gffLine.attr[idAttr]] = 0
            i = i + 1
        print(str(i)+' GFF lines processed')
            
        # Read in the SAM file
        readSeqFile = HTSeq.SAM_Reader(samFile)
        
        # Variables to store counts of read types
        alignOutCDS = 0
        nonUnique = 0
        i = 0   

        # Count reads
        for read in readSeqFile:
            i += 1
            # Skip ambiguous reads
            if read.optional_field("NH") > 1:
                nonUnique += 1
                continue
            iv_seq = (co.ref_iv for co in read.cigar if co.type == "M" and co.size > 0)
            fs = None
            for iv in iv_seq:
                for iv2, fs2 in features[iv].steps():
                    if len(fs2) > 0 or overlapMode == "intersection-strict":
                        if fs is None:
                            fs = fs2.copy()
                        else:
                            fs = fs.intersection(fs2)
            if fs is None or len(fs) == 0:
                alignOutCDS += 1
            else:
                countDict[list(fs)[0]] += 1
                # Write the read to file
                with open(readsDir+'/'+mt+'-'+list(fs)[0]+".reads", "a") as myfile:
                   myfile.write(read.original_sam_line)
            
        print(str(i)+' SAM alignments processed')      
            
        with open(countFolder+'/'+mt+'-'+genome+'.CDS.out', "w") as outFile:
            for fn in sorted(countDict.keys()):   
                outFile.write("%s\t%d\n" % (fn, countDict[fn]))
            outFile.write("__align_outside_CDS\t%d\n" % alignOutCDS)
            outFile.write("__alignment_not_unique\t%d" % nonUnique)
                        
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
cladeCogToCdsDF = cladeCogToCdsDF.sort(columns=None, axis=0)
cladeCogToCdsDF.to_csv(countFolder+'/cladesCogsToCDS.csv')

# Read in singly-indexed, drop empty rows, and write to file
cladeCogToCdsDF = pd.read_csv(countFolder+'/cladesCogsToCDS.csv', index_col=0)
cladeCogToCdsDF = cladeCogToCdsDF.dropna(axis=0, how='any')

cladeCogToCdsDF.to_csv(countFolder+'/cladesCogsToCDS.csv')

#%%#############################################################################
### Count total and unique reads which map to each (clade, group) pairing
### Requires integrating data from taxonomy and cog tables into a single data
###  structure
################################################################################

# Convert the cladeCogToCdsDF to a multi-index so we can construct the 
# count table for COGs

cladeCogToCdsDF = cladeCogToCdsDF.set_index([cladeCogToCdsDF.index, 'COG'])

# For each MT, construct the count table and write to file
for mt in mtList:
    for genome in genomeList:    
    # Reset columns for total and unique reads
        cladeCogToCdsDF['Total'] = 0

        for index in cladeCogToCdsDF.index:
            cdsList = cladeCogToCdsDF.loc[index, 'CDS'].split(',')
            readList = []
            for cds in cdsList:
                if os.path.isfile(readsDir+'/'+mt+'-'+cds+'.reads'):
                    readSeqFile = HTSeq.SAM_Reader(readsDir+'/'+mt+'-'+cds+'.reads')
                    for read in readSeqFile:
                        readList.append(read.read.name)
            cladeCogToCdsDF.loc[index, 'Total'] = len(readList)
        
        cladeCogToCdsDF.to_csv(countFolder+'/'+mt+'-'+genome+'.COG.out')
