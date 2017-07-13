###############################################################################
# auxotrophyCOGs.py
# Copyright (c) 2017, Joshua J Hamilton and Katherine D McMahon
# Affiliation: Department of Bacteriology
#              University of Wisconsin-Madison, Madison, Wisconsin, USA
# URL: http://http://mcmahonlab.wisc.edu/
################################################################################
# Identify genes and COGs associated with a given set of reactions.
################################################################################

#%%#############################################################################
### Import packages
################################################################################

#from collections import Counter
import cobra
#import copy
import os
#import pandas as pd
import re

#%%#############################################################################
### Define folder structure
################################################################################
    
modelDir = '../../models/rawModels'
rxnFile = '../../data/externalData/auxotrophies.csv'
resultsDir = '../../data/auxotrophies'
genomeDir = '../../data/auxotrophies/modelGenes'
cogDir = '../../data/orthoMCL/genomeCOGs'

# Check that the output directory exists. If not, create it.
if not os.path.exists(resultsDir):
    os.makedirs(resultsDir)

if not os.path.exists(genomeDir):
    os.makedirs(genomeDir)
    
# Import the list of models
modelList = []
for model in os.listdir(modelDir):
    if not model.startswith('.'):
        modelList.append(model)

##%%#############################################################################
#### Identify genes associated with auxotrophies
#################################################################################
#
#rxnList = []
#with open(rxnFile, 'r') as inFile:
#    for line in inFile:
#        rxnList.append(line.strip())
#        
#for curModel in modelList:
#    
#    # Read in model from SBML and create dict to store stuff
#    model = cobra.io.read_sbml_model(modelDir+'/'+curModel+'/'+curModel+'.xml')
#    transDict = {}    
#    
#    for curRxn in model.reactions:
#    
#        # If current reaction in the reaction list
#        rxnName = re.sub('_\w0', '', curRxn.id)
#        if rxnName in rxnList:
#            cdsList = []
#            for gene in curRxn.genes:
#                if gene.id != 'Unknown':
#                    cdsList = cdsList + [gene.id]
#            transDict[rxnName] = cdsList
#
#    with open(genomeDir+'/'+curModel+'.txt', 'w') as outFile:
#        for key in transDict.keys():
#            outFile.write(key+';')
#            for cds in transDict[key]:
#                outFile.write(cds+',')
#            outFile.write('\n')
    
#%%#############################################################################
### Identify the CDS and COGs associated with each reaction
################################################################################

# Identity all the unique reactions within the clade and their associated COGs.
    
geneDict = {}
cogDict = {}

for model in modelList:

    # Create a dictionary associating with COGs with each CDS
    # Temporary dict to store all associations for that genome
    tempDict = {}
    with open(cogDir+'/'+model+'COGs.txt', 'r') as inFile:
        for line in inFile:
            [cds, cog] = line.strip().split(',')
            tempDict[cds] = cog

    # Create a dictionary associating CDS with each reaction
    with open(genomeDir+'/'+model+'.txt', 'r') as inFile:
        for line in inFile:
            [rxn, cdsArray] = line.strip().split(';')
            cdsArray = cdsArray.split(',')
            
            # Remove empty strings and update CDS string format
            cdsArray = list(filter(None, cdsArray))
            cdsArray = [cds.replace('_CDS_', '.genome.CDS.') for cds in cdsArray]
            
            if len(cdsArray) > 0:
                if rxn in geneDict.keys():
                    geneDict[rxn] = geneDict[rxn] + cdsArray
                else:
                    geneDict[rxn] = cdsArray 
                # Transform cdsArray into a cogArray
                cogArray = [tempDict[cds] for cds in cdsArray]
                if rxn in cogDict.keys():
                    cogDict[rxn] = cogDict[rxn] + cogArray
                else:
                    cogDict[rxn] = cogArray 
                
with open(resultsDir+'/rxn-to-CDS.txt', 'w') as outFile:
    for key in geneDict.keys():
        outFile.write(key+';')
        for cds in geneDict[key]:
            outFile.write(cds+',')
        outFile.write('\n')
        
with open(resultsDir+'/rxn-to-COG.txt', 'w') as outFile:
    for key in cogDict.keys():
        outFile.write(key+';')
        for cds in cogDict[key]:
            outFile.write(cds+',')
        outFile.write('\n')
