###############################################################################
# parseCOGs.py
# Copyright (c) 2016, Joshua J Hamilton and Katherine D McMahon
# Affiliation: Department of Bacteriology
#              University of Wisconsin-Madison, Madison, Wisconsin, USA
# URL: http://http://mcmahonlab.wisc.edu/
# All rights reserved.
################################################################################
# Rearrange groups.txt output from OrthoMCL into a more usable format
################################################################################

#%%#############################################################################
### Import packages
################################################################################

from Bio import SeqIO
import pandas as pd

#%%#############################################################################
### User-defined files and folder structure
################################################################################
genomeFolder = '../genomes/faa'
resultsFolder = '../results'

#%%#############################################################################
### Create a hash for mapping names to taxon IDs
### Create a dataFrame to store the results
################################################################################

taxonDict = {}
with open('taxonMapping.txt') as dictFile:
    for line in dictFile:
       (val, key) = line.split() 
       taxonDict[key] = val
       
inCogDict = {}
with open(resultsFolder+'/groups.txt') as dictFile:
    for line in dictFile:
       (key, val) = line.split(':') 
       inCogDict[key] = val

groupNum = len(inCogDict)

with open(resultsFolder+'/singletons.txt') as singletonFile:
    for line in singletonFile:
        key = 'group'+str(groupNum).zfill(5)
        inCogDict[key] = line
        groupNum = groupNum + 1

outCogFrame = pd.DataFrame(index=inCogDict.keys(), columns=taxonDict.values())

#%%#############################################################################
### Parse the inCogDict
### For each key, split the line on whitespace
### For each element, split along the '|'
### Use the prefix to look up the genome in taxonDict
### Write the results to the appropriate dataFrame element
################################################################################

for key in inCogDict.keys():
    cogList = inCogDict[key].split()
    for cogLocus in cogList:
        code = cogLocus.split('|')[0]
        locus = cogLocus.split('|')[1]
#        locus = locus.split('.')[2]+'.'+locus.split('.')[3]
# Check if field already exists
        if pd.isnull(outCogFrame.loc[key, taxonDict[code]]):
            outCogFrame.loc[key, taxonDict[code]] = locus
        else:
            outCogFrame.loc[key, taxonDict[code]] = outCogFrame.loc[key, taxonDict[code]]+';'+locus

outCogFrame.to_csv(resultsFolder+'/cogTable.csv')

#%%#############################################################################
### Assign annotations to COGs
### For each genome (column of outCogFrame), read the annotation information
### into a hash.
### For each cog in that genome, look up the annotation and assign it to
### cogAnnotFrame
################################################################################

cogAnnotFrame = outCogFrame.copy()

for genome in outCogFrame.columns:
# Create the annotation hash
    annotHash = {}        
    inFile = open(genomeFolder+'/'+genome+'.faa', 'r')
    for record in SeqIO.parse(inFile, 'fasta'):
        locus = record.description.split()[1]
#        locus = locus.split('.')[2]+'.'+locus.split('.')[3]
        annotation = record.description.split()[2:]
        annotation = ' '.join(annotation)
        if not annotation:
            annotation = 'None Provided'
        annotHash[locus] = annotation
    
    for index in cogAnnotFrame.index:
        print('Processing genome: '+genome+' and locus: '+index)
        if not pd.isnull(cogAnnotFrame.loc[index, genome]):
            locusList = cogAnnotFrame.loc[index, genome].split(';')
            annotList = []
            for locus in locusList:
                annotList.append(annotHash[locus])
            cogAnnotFrame.loc[index, genome] = annotList
    inFile.close()
    
cogAnnotFrame.to_csv(resultsFolder+'/annotTable.csv')

#%%#############################################################################
### Extract uniqe annotations
### Create a empty DF indexed by groups
### For each group, extract unique annotations and drop 'nan'
### Add to dataframe
################################################################################

def flatten(x):
    result = []
    for el in x:
        if hasattr(el, "__iter__") and not isinstance(el, basestring):
            result.extend(flatten(el))
        else:
            result.append(el)
    return result

annotSummary = pd.DataFrame(index=inCogDict.keys(), columns=['Annotations'])

for group in cogAnnotFrame.index:
    annotMess = cogAnnotFrame.loc[group].tolist()
    annotList = flatten(annotMess)
    annotList = [annot for annot in annotList if str(annot) != 'nan']
    annotList = list(set(annotList))
    annotList = '; '.join(annotList)
    annotSummary.loc[group] = annotList.strip(';')

annotSummary.to_csv(resultsFolder+'/annotSummary.csv')
