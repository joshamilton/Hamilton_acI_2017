###############################################################################
# findGenesAnnotExpr.py
# Copyright (c) 2016, Joshua J Hamilton and Katherine D McMahon
# Affiliation: Department of Bacteriology
#              University of Wisconsin-Madison, Madison, Wisconsin, USA
# URL: http://http://mcmahonlab.wisc.edu/
# All rights reserved.
################################################################################
# Identify transport reactions and their corresponding genes. Map genes to COGs
# and extract their expression profiles.
################################################################################

#%%#############################################################################
### Import packages
################################################################################

from collections import Counter
import cobra
import copy
import os
import pandas as pd
import re

#%%#############################################################################
### Define folder structure
################################################################################
    
modelDir = '../../models/rawModels'
genomeDir = '../../data/transporters/modelGenes'
cladeDir = '../../data/transporters/cladeCOGs'
cogDir = '../../data/orthoMCL/genomeCOGs'
exprDir = '../../results/expression'
resultsDir = '../../results/transporters'

taxonFile = '../../data/externalData/taxonomy.csv'
annotTable = '../../data/orthoMCL/annotTable.csv'
    
#%%#############################################################################
### Model pre-processing
################################################################################

# Check that the output directory exists. If not, create it.
if not os.path.exists(genomeDir):
    os.makedirs(genomeDir)

if not os.path.exists(cladeDir):
    os.makedirs(cladeDir)

if not os.path.exists(resultsDir):
    os.makedirs(resultsDir)
    
# Import the list of models
modelList = []
for model in os.listdir(modelDir):
    if not model.startswith('.'):
        modelList.append(model)
        
#%%#############################################################################
### Identify transporter genes
################################################################################
        
for curModel in modelList:
    
    # Read in model from SBML and create dict to store stuff
    model = cobra.io.read_sbml_model(modelDir+'/'+curModel+'/'+curModel+'.xml')
    transDict = {}    
    
    for curRxn in model.reactions:
        
        # Transport reactions, based on keywords
        if re.search('transport', curRxn.name) or re.search('permease', curRxn.name) or re.search('symport', curRxn.name) or re.search('diffusion', curRxn.name) or re.search('excretion', curRxn.name) or re.search('export', curRxn.name) or re.search('secretion', curRxn.name) or re.search('uptake', curRxn.name) or re.search('antiport', curRxn.name):
            cdsList = []
            for gene in curRxn.genes:
                if gene.id != 'Unknown':
                    cdsList = cdsList + [gene.id]
            transDict[curRxn.id] = cdsList
            
        # Transport reactions which don't get picked up based on keywords
        elif curRxn.id == 'rxn05226_c0' or curRxn.id == 'rxn05292_c0' or curRxn.id == 'rxn05305_c0' or curRxn.id == 'rxn05312_c0' or curRxn.id == 'rxn05315_c0' or curRxn.id == 'rxn10945_c0' or curRxn.id == 'rxn10116_c0':
            cdsList = []
            for gene in curRxn.genes:
                if gene.id != 'Unknown':
                    cdsList = cdsList + [gene.id]
            transDict[curRxn.id] = cdsList
            

    with open(genomeDir+'/'+curModel+'.txt', 'w') as outFile:
        for key in transDict.keys():
            outFile.write(key+';')
            for cds in transDict[key]:
                outFile.write(cds+',')
            outFile.write('\n')

#%%#############################################################################
### For each clade, identify the COGs associated with each reaction
### For each (rxn, cog) pairing, identify the expression data
################################################################################

# Read in the taxonomy table and create a list of genomes for each clade
cladeToGenomeDict = {}
cladeList = []

taxonClass = pd.DataFrame.from_csv(taxonFile, sep=',')
taxonClass = taxonClass.dropna()

# Extract the unique clades
cladeList = pd.unique(taxonClass['Clade'].values)

for clade in cladeList:
    genomeList = taxonClass[taxonClass['Clade'] == clade].index.tolist()
    cladeToGenomeDict[clade] = genomeList

# Read in the annotation table
annotDF = pd.read_csv(annotTable, index_col=0)

# Identity all the unique reactions within the clade and their associated COGs.
for clade in cladeList:
    
    geneDict = {}   
    cogDict = {}    
    
    modelList = cladeToGenomeDict[clade]

    for model in modelList:

        # Create a dictionary associating CDS with each reaction
        with open(genomeDir+'/'+model+'.txt', 'r') as inFile:
            for line in inFile:
                [gene, cdsArray] = line.strip().split(';')
                cdsArray = cdsArray.split(',')
                cdsArray = filter(None, cdsArray)
                if len(cdsArray) > 0:
                    if gene in geneDict.keys():
                        geneDict[gene] = geneDict[gene] + cdsArray
                    else:
                        geneDict[gene] = cdsArray                
                    for cds in cdsArray:
                        cogDict[cds] = None

        # Create a dictionary associating with COGs with each CDS
        # Temporary dict to store all associations for that genome
        tempDict = {}
        with open(cogDir+'/'+model+'COGs.txt', 'r') as inFile:
            for line in inFile:
                [cds, cog] = line.strip().split(',')
                tempDict[cds] = cog

        # Populate the cogDict using this info
        for cds in cogDict.keys():
            if cds.replace('_CDS_', '.genome.CDS.') in tempDict.keys():
                cogDict[cds] = tempDict[cds.replace('_CDS_', '.genome.CDS.')]
    
    with open(cladeDir+'/'+clade+'.CDS.txt', 'w') as outFile:
        for key in geneDict.keys():
            outFile.write(key+';')
            for cds in geneDict[key]:
                outFile.write(cds+',')
            outFile.write('\n')
                
    # Now, we need to map the CDS for each reaction to its COG.
    rxnCogDict = copy.deepcopy(geneDict)
    
    for rxn in rxnCogDict.keys():
        for pos, cds in enumerate(rxnCogDict[rxn]):
           rxnCogDict[rxn][pos] = cogDict[cds]
       
    # Some CDS map to the same COG, so update the lists to only include
    # unique entries
    for rxn in rxnCogDict.keys():
        rxnCogDict[rxn] = list(set(rxnCogDict[rxn]))

    with open(cladeDir+'/'+clade+'.COG.txt', 'w') as outFile:
        for key in sorted(rxnCogDict.keys()):
            for cds in sorted(rxnCogDict[key]):
                outFile.write(key+','+str(cds)+'\n')
                
    # Now, read in the expression data for that clade
    exprDataFrame = pd.read_csv(exprDir+'/'+clade+'.norm', index_col=1)
    exprDataFrame = exprDataFrame.drop('Clade', axis=1)

    # Create an empty dataframe
    rxnCogExprMultiIndex = pd.MultiIndex(levels=[[],[]],
                             labels=[[],[]],
                             names=['Reaction', 'COG'])
                             
    rxnCogExprDataFrame = pd.DataFrame(index=rxnCogExprMultiIndex, columns=exprDataFrame.columns)
                             
    # Iterate over the rxnCogDict and look up expression values in the exprDataFrame
    # Use these to populate the rxnCogExprDataFrame
    for rxn in sorted(rxnCogDict.keys()):
            for cds in sorted(rxnCogDict[rxn]):
                # If CDS IS in the genome AND expressed
                if cds in exprDataFrame.index:
                    rxnCogExprDataFrame.loc[(rxn, cds),:] = exprDataFrame.loc[cds]
                # If CDS IS in the genome AND NOT expressed
                elif cds in cogDict.values():
                    rxnCogExprDataFrame.loc[(rxn, cds),:] = 0
                # If CDS IS NOT in the genome
                else:
                    rxnCogExprDataFrame.loc[(rxn, cds),:] = None

    # The genes which are not expressed will not have consensus annotations
    # Rerun that piece of code

    # Compute majority annotation
    # First subset the dataframe, keep genomes for that clade and dropping 
    
    tempDF = annotDF[modelList]
    
    for curIndex in rxnCogExprDataFrame.index:
        cog = curIndex[1]
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
        rxnCogExprDataFrame.loc[curIndex,'Annotation'] = majorityAnnot
        
    # Write the results to file
    rxnCogExprDataFrame.to_csv(cladeDir+'/'+clade+'.COG.norm')
                    
#%%#############################################################################
### Aggregate it all into a single dataframe
################################################################################

# Create a master dataframe
masterMultiIndex = pd.MultiIndex(levels=[[],[]],
                             labels=[[],[]],
                             names=['Reaction', 'COG'])
                             
masterDataFrame = pd.DataFrame(index=rxnCogExprMultiIndex)

for clade in cladeList:
    
    # Read in the expression data for that clade
    rxnCogExprDataFrame = pd.read_csv(cladeDir+'/'+clade+'.COG.norm', index_col=[0,1])
    
    # Rename the columns
    for column in rxnCogExprDataFrame.columns:
        rxnCogExprDataFrame.rename(columns={column:column+' ('+clade+')'}, inplace=True)
    
    # Merge into the masterDataFrame
    masterDataFrame = pd.concat([masterDataFrame, rxnCogExprDataFrame], axis=1, join='outer')

# Write to file
masterDataFrame.to_csv(resultsDir+'/transAnnotExpr.csv')
