###############################################################################
# kBaseGenbankToGFF.py
# Copyright (c) 2016, Joshua J Hamilton and Katherine D McMahon
# Affiliation: Department of Bacteriology
#              University of Wisconsin-Madison, Madison, Wisconsin, USA
# URL: http://http://mcmahonlab.wisc.edu/
# All rights reserved.
################################################################################
# Convert Genbank files downloaded from KBase into GFF files.
################################################################################

#%%#############################################################################
### Import packages
################################################################################

from BCBio import GFF
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature
import os

#%%#############################################################################
### Define folder structure
################################################################################
gbkDir = '../gbk'
gffDir = '../gff'


#%%#############################################################################
### Create list of genomes to process
################################################################################

genomeList = []
for genome in os.listdir(gbkDir):
    if genome.endswith('.gbk'):
        genomeList.append(genome)

genomeList = [genome.replace('.gbk', '') for genome in genomeList]

#%%#############################################################################
### Convert Genbank files to GFF files
### These files are necessary for counting mapped readsconti
################################################################################

#    seqid -- SeqRecord ID
#    source -- Feature qualifier with key "source"
#    type -- Feature type attribute
#    start, end -- The Feature Location
#    score -- Feature qualifier with key "score"
#    strand -- Feature strand attribute
#    phase -- Feature qualifier with key "phase" 
    
for genome in genomeList:
    inFile = open(gbkDir+'/'+genome+'.gbk')
    outFile = open(gffDir+'/'+genome+'.gff', "w")
    
    for gbkRecord in SeqIO.parse(inFile, "genbank"):
        for contigFeature in gbkRecord.features:
# The Genbank record is the entire contig
# Check the SeqFeature type to see if we want to procedue
            if contigFeature.type in ('CDS', 'misc_RNA', 'tRNA'):
            
                gffSeq = gbkRecord.seq
                gffRecord = SeqRecord(gffSeq, gbkRecord.id)
                
                gffQualifiers = {"source": "feature", "ID": contigFeature.qualifiers['gene'], "locus_tag": contigFeature.qualifiers['gene'],"Product": contigFeature.qualifiers['function']}
                
                gffFeature = SeqFeature(contigFeature.location, type=contigFeature.type,
                             qualifiers=gffQualifiers)

                gffRecord.features = [gffFeature]
                                    
                GFF.write([gffRecord], outFile)
            else:   
                next
    
    inFile.close()
    outFile.close()
