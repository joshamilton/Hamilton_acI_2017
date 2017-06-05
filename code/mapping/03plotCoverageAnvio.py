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
import os
import subprocess

#%%#############################################################################
### Static folder structure
################################################################################

# Define fixed input and output files
concatFolder = '../../data/refGenomes/concat'
genomeFolder = '../../data/refGenomes/fnaKBase'
genomeForAnvioFolder = '../../data/refGenomes/fnaAnvio'
gffFolder = '../../data/refGenomes/gff'
gffForAnvioFolder = '../../data/refGenomes/gffAnvio'
sampleFolder = '../../data/sequences'
sampleFolder = '../../data/sequences'
bamFolder = '../../data/mapping/bamFiles'
anvioFolder = '../../data/mapping/anvio'

# Check that the new output directory exists and create if it doesn't
if not os.path.exists(genomeForAnvioFolder):
        os.makedirs(genomeForAnvioFolder)

# Check that the new output directory exists and create if it doesn't
if not os.path.exists(gffForAnvioFolder):
        os.makedirs(gffForAnvioFolder)

# Check that the new output directory exists and create if it doesn't
if not os.path.exists(anvioFolder):
        os.makedirs(anvioFolder)
        
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
### Create FASTA files with properly-formatted deflines
################################################################################

print('Formatting FASTA files ...')

for genome in genomeList:
    with open(genomeForAnvioFolder+'/'+genome+'.fna', 'w') as outFile:
        for curSeq in SeqIO.parse(genomeFolder+'/'+genome+'.fna', 'fasta'):
            outFile.write('>'+curSeq.id.replace('.', '_')+'\n')
            outFile.write(str(curSeq.seq)+'\n')

#%%#############################################################################
### Create the '--external-gene-calls' files for anvio
### These are built from my GFF files
################################################################################

print('Formatting GFF files ...')

for genome in genomeList:
    curID = 1

    with open(gffForAnvioFolder+'/'+genome+'.gff', 'w') as outFile:
        outFile.write('gene_callers_id\tcontig\tstart\tstop\tdirection\tpartial\tsource\tversion\n')
        with open(gffFolder+'/'+genome+'.gff', 'r') as inFile:
            next(inFile)
            for line in inFile:
                gffArray = line.split('\t')
    
                # Obtain values for the fields
                gene_callers_id = curID
                contig = gffArray[0].replace('.contig.', '_contig_')
                start = int(gffArray[3]) - 1
                stop = int(gffArray[4])
                if gffArray[6] == '+':
                    direction = 'f'
                elif gffArray[6] == '-':
                    direction = 'r'
                if gffArray[2] != 'CDS':
                    partial = '1'
                elif (stop - start) % 3 == 0:
                    partial = '0'
                else:
                    partial = '1'
                source = 'Kbase'
                version = '1'
                
                # Write to file
                outFile.write(str(gene_callers_id)+'\t'+contig+'\t'+str(start)+'\t'+str(stop)+'\t'+direction+'\t'+partial+'\t'+source+'\t'+version+'\n')
                curID = curID + 1

#%%#############################################################################
### Create bam files for use with anvio
################################################################################

# The following samtools commands will convert the sample-acI.bam file to
# sorted and indexed bam files for each genome
# grep 'MEint885' ME150263.sam > ME150263-genome.sam
# samtools view -H ME150263-acI.bam > ME150263-header.sam
# cat ME150263-header.sam ME150263-genome.sam > foo.sam
# mv foo.sam ME150263-genome.sam
# samtools view -bS ME150263-genome.sam | samtools sort -o ME150263-genome.bam

for concat in concatList:
    for sample in sampleList:
        
        print('Formatting BAM files for sample '+sample+' ...')
            
        # samtools view ME150263-acI.bam > ME150263.sam
        # samtools view -H ME150263-acI.bam > ME150263-header.sam
        subprocess.call('samtools view '+bamFolder+'/'+sample+'-'+concat+'.bam >'+bamFolder+'/'+sample+'.sam', shell=True)
        subprocess.call('samtools view -H '+bamFolder+'/'+sample+'-'+concat+'.bam >'+bamFolder+'/'+sample+'-header.sam', shell=True)

        for genome in genomeList:
            # The header file needs to get cleaned up - for anvio, should only include
            with open(bamFolder+'/'+sample+'-'+genome+'-header.sam', 'w') as outFile:
                with open(bamFolder+'/'+sample+'-header.sam', 'r') as inFile:
                    fileArray = inFile.readlines()
                    outFile.write(fileArray[0]) # write first line
                    # Check if line contains a contig for our genome and write to file.
                    for curLine in fileArray:
                        if curLine.find(genome) != -1:
                            outFile.write(curLine)
                    outFile.write(fileArray[len(fileArray) - 1]) # write last line

            subprocess.call('grep \''+genome+'\' '+bamFolder+'/'+sample+'.sam >'+bamFolder+'/'+sample+'-'+genome+'.sam', shell=True)

            # Alternative to 'cat'
            catFiles = [bamFolder+'/'+sample+'-'+genome+'-header.sam', bamFolder+'/'+sample+'-'+genome+'.sam']
            with open(bamFolder+'/foo.sam', 'w') as outFile:
                for curFile in catFiles:
                    with open(curFile) as inFile:
                        for line in inFile:
                            outFile.write(line)

            # Alternative to 'mv'
            os.rename(bamFolder+'/foo.sam', bamFolder+'/'+sample+'-'+genome+'.sam')
            
            # Before converting to BAM, update contig names
            # Rename .contig. to _contig_
            
            with open(bamFolder+'/'+sample+'-'+genome+'.sam', 'r') as inFile:
                with open(bamFolder+'/foo.sam', 'w') as outFile:
                    for line in inFile:
                        outFile.write(line.replace('.contig.', '_contig_'))
            
            os.rename(bamFolder+'/foo.sam', bamFolder+'/'+sample+'-'+genome+'.sam')
            
            subprocess.call('samtools view -bS '+bamFolder+'/'+sample+'-'+genome+'.sam | samtools sort -o '+bamFolder+'/'+sample+'-'+genome+'.bam', shell=True)
            
            # Finally, index the bam file
            subprocess.call('samtools index '+bamFolder+'/'+sample+'-'+genome+'.bam '+bamFolder+'/'+sample+'-'+genome+'.bai', shell=True)
            os.remove(bamFolder+'/'+sample+'-'+genome+'.sam')

            os.remove(bamFolder+'/'+sample+'-'+genome+'-header.sam')

        os.remove(bamFolder+'/'+sample+'.sam')
        os.remove(bamFolder+'/'+sample+'-header.sam')
        
    
#%%#############################################################################
### Make the anvio visualization
################################################################################

for genome in genomeList:
    
    print('Running anvio for genome '+genome+' ...')
            
    # Make a folder
    contigFolder = anvioFolder+'/'+genome
    if not os.path.exists(contigFolder):
        os.makedirs(contigFolder)

    # Make the contigs database
    subprocess.call(['anvi-gen-contigs-database',
                     '-f', genomeForAnvioFolder+'/'+genome+'.fna',
                     '-o', contigFolder+'/'+genome+'.db',
                     '--external-gene-calls', gffForAnvioFolder+'/'+genome+'.gff',
                     '--split-length', '-1'])

    # Profile each sample
    for sample in sampleList:
        # Make a folder
        profileFolder = contigFolder+'/'+sample
        if not os.path.exists(profileFolder):
            os.makedirs(profileFolder)

        subprocess.call(['anvi-profile',
                         '-i', bamFolder+'/'+sample+'-'+genome+'.bam',
                         '--sample-name', sample+'_'+genome,
                         '-c', contigFolder+'/'+genome+'.db',
                         '--output-dir', profileFolder,
                         '--overwrite-output-destinations',
                         '--skip-SNV-profiling',
                         '--min-contig-length', '0'])
    
    # Make a folder
    mergeFolder = contigFolder+'/merged'
    if not os.path.exists(mergeFolder):
        os.makedirs(mergeFolder)
    
    mergeList = []
    for sample in sampleList:
        mergeList.append(contigFolder+'/'+sample+'/PROFILE.db')
    mergeString = ' '.join(mergeList)

    subprocess.call('anvi-merge ' + mergeString + 
                     ' -o ' + mergeFolder + 
                     ' --sample-name ' + genome + 
                     ' -c ' + contigFolder+'/'+genome+'.db' +
                     ' --skip-concoct-binning' +
                     ' --overwrite-output-destinations', shell=True)

