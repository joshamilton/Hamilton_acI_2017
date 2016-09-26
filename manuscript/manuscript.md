---
title: Metabolic Network Analysis and Metatranscriptomics of a Cosmopolitan and Streamlined Freshwater Lineage
author: Joshua J. Hamilton^1,\*^, Sarahi L. Garcia^2^, Brittany S. Brown^1^, Jeffrey R. Dwulit-Smith^1^, Francisco Moya^4^, Ben O. Oyserman^4^, Sarah L.R. Stevens^1^, Stefan Bertilsson^2^, Katrina T. Forest^1^, Susannah G, Tringe^3^, Tanja Woyke^3^, and Katherine D. McMahon^1,4,\*^
abstract: \* Correspondence&#58; Joshua J. Hamilton, jjhamilton2@wisc.edu and Katherine D. McMahon, trina.mcmahon@wisc.edu
date: ^1^ Department of Bacteriology, University of Wisconsin-Madison, Madison, WI, USA; ^2^ Department of Ecology and Genetics, Uppsala University, Uppsala, Sweden; ^3^ United States Department of Energy Joint Genome Institute, Walnut Creek, CA, USA; ^4^ Department of Civil and Environmental Engineering, University of Wisconsin-Madison, Madison, WI, USA
---

# Abstract

# Introduction

# Materials and Methods

## Single-Cell Genome Generation, Selection, and Sequencing

Single-cell genomes were collected from the top of the water column (depth <1m) from each of two lakes, Mendota (Madison, WI, USA) and Damariscotta (Lincoln County, ME USA), in 2009. Samples were cryopreserved and sent to the Single Cell Genomics Center at the Bigelow Laboratory for Ocean Sciences for sorting, as previously described [@Martinez-Garcia2012, @Garcia2013]. Partial 16S rRNA genes amplified previously [@Martinez-Garcia2012] were phylogenetically classified using a controlled nomenclature for freshwater bacteria [@Newton2011a] by insertion into references trees created in the ARB software package [@Ludwig2004].

Actinobacterial SAGs used in this study were then sent to the JGI for sequencing and assembly, also as previously described [@Ghylin2014]. Briefly, shotgun libraries were constructed for each of the SAGs from re-amplified MDA products and sequenced on an Illumina HiSeq2000. All general aspects of and detailed protocols for library construction and sequencing can be found on the JGI website (http://www.jgi.doe.gov/).

For assembly, raw sequence data was first passed through a filtering program developed at JGI to eliminate known sequencing and library preparation artifacts. Assembly was then performed using Velvet [@Zerbino2008 and ALLPATHS-LG [@Butler2008]. Additional details of the assembly process have been previously described [@Ghylin2014] and are available through the JGI Genome Portal (http://genome.jgi.doe.gov) Genome sequences are available through IMG (https://img.jgi.doe.gov/cgi-bin/mer/main.cgi). Genome-specific information can be accessed in both databases by searching for the IMG Taxon OIDs given in Table 1.

## Metagenome Sampling, Sequencing, Assembly, and Binning

Sample collection, DNA sequencing, metabgenomic assembly, and genomic binning for the Trout Bog samples have been described previously [@Bendall2016], and similar procedures were followed for Lake Mendota samples. A summary is provided here.

For Lake Mendota, Depth-integrated water samples were collected from the top 12 meters at 96 time points during ice-free periods from 2008 to 2011. For Trout Bog, depth-integrated water samples were collected from the epilimnion (44 samples) and hypolimnion (45 samples) layers during ice-free periods from 2007 to 2009. All samples were filtered on 0.2 μm polyethersulfone filters (Supor, Pall Corp) prior to storage at -80°C, as described previously [@Bendall2016].. DNA was extracted from these filters using the FastDNA kit (MP Biomedicals) and sent to the JGI for sequencing, as described previously [@Bendall2016].

Shotgun libraries were constructed for each of the samples and sequenced on an Illumina GA IIx (four Trout Bog samples) or an Illumina HiSeq2000 (all other samples), following a 2x150 indexed run recipe as previously described [@Bendall2016]. All general aspects of and detailed protocols for library construction and sequencing can be found on the JGI website (http://www.jgi.doe.gov/). Metagenomic sequence reads are publicly available on the JGI Genome Portal (http://genome.jgi.doe.gov/) under Proposal ID 394.

Raw sequence data was passed through a filtering program developed at JGI to eliminate known sequencing and library preparation artifacts. Prior to assembly, reads were merged with FLASH [@Magoc2011], as previously described [@Bendall2016]. Merged reads were pooled by lake and layer into three co-assemblies using SOAPdenovo [@Luo2012a], and contigs from the resulting assemblies were assembled into a final assembly using Minimus [@Sommer2007], as previously described [@Bendall2016]. Additional details of the assembly process and metagenomic sequence reads are available through the JGI Genome Portal (http://genome.jgi.doe.gov) under Proposal ID 394.

Genomes were binned from each metagenomic co-assembly using MetaBat [@Kang2015], as described previously [@Bendall2016]. Briefly, contigs were classified into bins using tetranucleotide frequency and coverage patterns across the time-series and then manually curated, as previously described [@Bendall2016]. Genome sequences are available through IMG (https://img.jgi.doe.gov/cgi-bin/mer/main.cgi) by searching for the IMG Taxon OIDs given in Table 1. Genomes were classified using taxonomic assignments from a set of 37 highly-conserved single-copy marker genes using Phylosift [@Darling2014], as previously described [@Bendall2016]. Final bin size and number of contigs are reported in Table 1.

# Results

# Discussion

# Acknowledgements

# Conflict of Interest

The authors declare no conflict of interest.

# References
