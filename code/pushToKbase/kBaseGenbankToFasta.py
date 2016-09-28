###############################################################################
# kBaseGenbankToFasta.py
# Copyright (c) 2016, Joshua J Hamilton and Katherine D McMahon
# Affiliation: Department of Bacteriology
#              University of Wisconsin-Madison, Madison, Wisconsin, USA
# URL: http://http://mcmahonlab.wisc.edu/
# All rights reserved.
################################################################################
# Code description.
################################################################################

#%%#############################################################################
### Import packages
################################################################################

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
import os

#%%#############################################################################
### Define folder structure
################################################################################
faaDir = '../ffn'
ffnDir = '../ffn'
gbkDir = '../gbk'

#%%#############################################################################
### Create list of genomes to process
################################################################################

genomeList = []
for genome in os.listdir(gbkDir):
    if genome.endswith('.gbk'):
        genomeList.append(genome)

genomeList = [genome.replace('.gbk', '') for genome in genomeList]

#%%#############################################################################
### Convert Genbank files to FASTA nucleotide files
### These are necessary to run clustering algorithms which require input gene
### sequences.
################################################################################

for genome in genomeList:
    inFile = open(gbkDir+'/'+genome+'.gbk')
    outFile = open(ffnDir+'/'+genome+'.ffn', "w")
    
    for gbkRecord in SeqIO.parse(inFile, "genbank"):
        for contigFeature in gbkRecord.features:
            if contigFeature.type in ('CDS', 'misc_RNA', 'tRNA'):
# Print the FASTA header file, with contig name, locus tag, and gene product
                print >>outFile, ">"+gbkRecord.id+" "+contigFeature.qualifiers['gene'][0]+" "+contigFeature.qualifiers['function'][0]
# Extract the gene sequence based on its position. Use 'strand' information to check if we need to retrieve the reverse-complement.
                geneStart = contigFeature.location.start
                geneStop = contigFeature.location.end
                strand = contigFeature.location.strand
                if strand == -1:
                    geneSequence = gbkRecord.seq[geneStart:geneStop].reverse_complement()
                else :
                    geneSequence = gbkRecord.seq[geneStart:geneStop]
                print >>outFile, geneSequence+"\n"
            next;
 
    inFile.close()
    outFile.close()

for genome in genomeList:
    inFile = open(gbkDir+'/'+genome+'.gbk')
    outFile = open(faaDir+'/'+genome+'.faa', "w")
    
    for gbkRecord in SeqIO.parse(inFile, "genbank"):
        for contigFeature in gbkRecord.features:
            if "translation" in contigFeature.qualifiers:
# Print the FASTA header file, with contig name, locus tag, and gene product
                print >>outFile, ">"+gbkRecord.id+" "+contigFeature.qualifiers['gene'][0]+" "+contigFeature.qualifiers['function'][0]
# Extract the gene sequence based on its position. Use 'strand' information to check if we need to retrieve the reverse-complement.
                print >>outFile, contigFeature.qualifiers["translation"][0]+"\n"
 
    inFile.close()
    outFile.close()
