---
title: Metabolic Network Analysis and Metatranscriptomics of a Cosmopolitan and Streamlined Freshwater Lineage
author: Joshua J. Hamilton^1,\*^, Sarahi L. Garcia^2^, Brittany S. Brown^1^, Jeffrey R. Dwulit-Smith^1^, Francisco Moya^4^, Ben O. Oyserman^4^, Sarah L.R. Stevens^1^, Stefan Bertilsson^2^, Katrina T. Forest^1^, Susannah G, Tringe^3^, Tanja Woyke^3^, and Katherine D. McMahon^1,4,\*^
abstract: \* Correspondence&#58; Joshua J. Hamilton, jjhamilton2@wisc.edu and Katherine D. McMahon, trina.mcmahon@wisc.edu
date: ^1^ Department of Bacteriology, University of Wisconsin-Madison, Madison, WI, USA; ^2^ Department of Ecology and Genetics, Uppsala University, Uppsala, Sweden; ^3^ United States Department of Energy Joint Genome Institute, Walnut Creek, CA, USA; ^4^ Department of Civil and Environmental Engineering, University of Wisconsin-Madison, Madison, WI, USA
---

# Abstract

# Introduction

## acI

### What we know from tag data
* persistent and abundant through time and across lakes
* three clades and thirteen tribes

### What we know from SAG data
* small and streamlined genome
* conclusions on substrate utilization from SAG studies, incl. rhodopsin
* conclusions on substrate utilization from FISH studies
* no genomes from acI-C

## Reverse Ecology
* genome-gazing is a bottleneck, FISH studies are laborious - need an alternative
* reverse-ecology is an emerging paradigm to infer ecology from genomic information w/ a focus on systems-level signatures
* metabolism one aspect of ecology, long tradition of metabolic network analysis
* example of reverse ecology - seed compounds and environment

## Overview

* assembled largest catalog of acI genomes, including acI-C
* use existing reverse ecology to predict seed compounds for acI genomes
  * high-throughput approach to leverage incomplete genome data
* confirm conserved gene functions and revealed clade-specific differences
* integrate with metatranscriptomic data
* many diverse transporters, all highly expressed. possible key to success of passive & planktonic lifestyle

# Materials and Methods

## Single-Cell Genome Generation, Selection, and Sequencing

Single-cell genomes were collected from the top of the water column (depth <1m) from each of two lakes, Mendota (Madison, WI, USA) and Damariscotta (Lincoln County, ME USA), in 2009. Samples were cryopreserved and sent to the Single Cell Genomics Center at the Bigelow Laboratory for Ocean Sciences for sorting, as previously described [@Martinez-Garcia2012, @Garcia2013]. Partial 16S rRNA genes amplified previously [@Martinez-Garcia2012] were phylogenetically classified using a controlled nomenclature for freshwater bacteria [@Newton2011a] by insertion into references trees created in the ARB software package [@Ludwig2004].

Actinobacterial SAGs used in this study were then sent to the JGI for sequencing and assembly, also as previously described [@Ghylin2014]. Briefly, shotgun libraries were constructed for each of the SAGs from re-amplified MDA products and sequenced on an Illumina HiSeq2000. All general aspects of and detailed protocols for library construction and sequencing can be found on the JGI website (http://www.jgi.doe.gov/).

For assembly, raw sequence data was first passed through a filtering program developed at JGI to eliminate known sequencing and library preparation artifacts. Assembly was then performed using Velvet [@Zerbino2008 and ALLPATHS-LG [@Butler2008]. Additional details of the assembly process have been previously described [@Ghylin2014] and are available through the JGI Genome Portal (http://genome.jgi.doe.gov) Genome sequences are available through IMG (https://img.jgi.doe.gov/cgi-bin/mer/main.cgi). Genome-specific information can be accessed in both databases by searching for the IMG Taxon OIDs given in Table 1.

## Metagenome Sampling, Sequencing, Assembly, and Binning

Sample collection, DNA sequencing, metabgenomic assembly, and genomic binning for the Trout Bog samples have been described previously [@Bendall2016], and similar procedures were followed for Lake Mendota samples. A summary is provided here.

For Lake Mendota, depth-integrated water samples were collected from the top 12 meters at 96 time points during ice-free periods from 2008 to 2011. For Trout Bog, depth-integrated water samples were collected from the epilimnion (44 samples) and hypolimnion (45 samples) layers during ice-free periods from 2007 to 2009. All samples were filtered on 0.2 μm polyethersulfone filters (Supor, Pall Corp) prior to storage at -80°C, as described previously [@Bendall2016]. DNA was extracted from these filters using the FastDNA kit (MP Biomedicals) and sent to the JGI for sequencing, as described previously [@Bendall2016].

Shotgun libraries were constructed for each of the samples and sequenced on an Illumina GA IIx (four Trout Bog samples) or an Illumina HiSeq2000 (all other samples), following a 2x150 indexed run recipe as previously described [@Bendall2016]. All general aspects of and detailed protocols for library construction and sequencing can be found on the JGI website (http://www.jgi.doe.gov/). Metagenomic sequence reads are publicly available on the JGI Genome Portal (http://genome.jgi.doe.gov/) under Proposal ID 394.

Raw sequence data was passed through a filtering program developed at JGI to eliminate known sequencing and library preparation artifacts. Prior to assembly, reads were merged with FLASH [@Magoc2011], as previously described [@Bendall2016]. Merged reads were pooled by lake and layer into three co-assemblies using SOAPdenovo [@Luo2012a], and contigs from the resulting assemblies were assembled into a final assembly using Minimus [@Sommer2007], as previously described [@Bendall2016]. Additional details of the assembly process and metagenomic sequence reads are available through the JGI Genome Portal (http://genome.jgi.doe.gov) under Proposal ID 394.

Genomes were binned from each metagenomic co-assembly using MetaBat [@Kang2015], as described previously [@Bendall2016]. Briefly, contigs were classified into bins using tetranucleotide frequency and coverage patterns across the time-series and then manually curated, as previously described [@Bendall2016]. Genome sequences are available through IMG (https://img.jgi.doe.gov/cgi-bin/mer/main.cgi) by searching for the IMG Taxon OIDs given in Table 1. Genomes were classified using taxonomic assignments from a set of 37 highly-conserved single-copy marker genes using Phylosift [@Darling2014], as previously described [@Bendall2016]. acI genomes (class Actinobacteria) were further classified using a defined ontology for freshwater bacteria [@Newton2011a] as described below. Final bin size and number of contigs are reported in Table 1.

## Metatranscriptome Sampling and Sequencing

Four samples were collected from the top of the water column (depth <1m) from Lake Mendota (Madison, WI, USA) over a twenty-four hour period on August 20 and 21, 2015. For each sample, between 200 and 400 mL lake water was filtered onto a 0.2 μm polyethersulfone filter (Supor, Pall Corp) and stored at -80°C until extraction.

Prior to extraction, three samples were spiked with an internal standard to enable quantification of total transcript abundance, following an established protocol [@Satinsky2014a]. Briefly, a 970-nucleotide-long mRNA standard was synthesized using a a T7 RNA polymerase and the Riboprobe In Vitro Transcription System (Promega, Madison, WI), according to the manufacturer’s protocol. A fixed quantity of each standard (1.172 x 10^10 copies) was added independently to each lysis tube immediately prior to the addition of the sample filter.

Samples were subject to TRIzol-based RNA extraction (Thermo Fisher Scientific, Waltham, MA) followed by on-column DNAse digestion and RNA purification using an RNeasy Mini Kit (Qiagen, Venlo, Netherlands). RNA was then sent to the University of Wisconsin-Madison Biotechnology Center (https://www.biotech.wisc.edu) for sequencing. There, samples were prepared for sequencing using the TruSeq RNA Library Prep Kit v2 (Illumina, San Diego, CA), with the addition of a step for selective ribosomal RNA depletion using the Ribo-Zero rRNA Removal Kit (Bacteria) (Illumina). The resulting cDNA libraries were pooled in an equimolar ratio, and sequenced on an Illumina HiSeq2500.

Raw paired-end reads were then merged using FLASH [@Magoc2011] using default parameters. Finally, additional rRNA and ncRNA sequences were removed using SortMeRNA [@Kopylova2012] using default parameters. SortMeRNA was run using eight built-in databases for bacterial, archaeal, and eukaryotic small and large ribosomal subunits and ncRNAs, derived from the SILVA 119 [@Quast2013] and RFAM [@Nawrocki2015] databases.

Additional information, including all protocols and scripts for RNA analysis, can be found on Github (https://github.com/McMahonLab/OMD-TOILv2). Raw RNA sequences can be found on the National Center for Biotechnology Information (NCBI) website under BioProject PRJNA######.

## Genome Completeness and Phylogenetic Relationships

CheckM [@Parks2015] was used to estimate genome completeness based on 204 single-copy marker genes conserved across the phylum Actinobacteria. Phylogenetic analysis of Actinobacterial SAGs and MAGs was performed using a concatenated alignment of single-copy marker genes obtained via Phylosift [@Darling2014]. Maximum likelihood trees were generated using RAxML [@Stamatakis2014] using the automatic protein model assignment option (PROTGAMMAAUTO) and 100 bootstraps.

## Metabolic Network Reconstruction and Reverse Ecology

### Genome Annotation and Model Processing (including Genome Merging)

Genome annotations and metabolic network reconstructions were performed using KBase (http://kbase.us/). Unannotated contigs for each genome were pushed to KBase and annotated using the "Annotate Microbial Contigs" method using default options, which uses components of the RAST toolkit [@Brettin2015, @Overbeek2014] for genome annotation. Genome-scale metabolic network reconstructions were performed using the "Build Metabolic Model" app using default parameters, which relies on the Model SEED framework [@Henry2010a] to build a draft metabolic model without gapfilling. As acI genomes are known to be small with numerous auxotrophies [@Garcia2015], we opted to forgo the gapfilling step and instead leverage multiple genomes to fill gaps, as described below.

Metabolic models were then downloaded via KBase, pruned, and converted to metabolic network graphs. In particular, biomass, exchange, transport, spontaneous, and DNA/RNA biosynthesis reactions were removed from the model, and reactions were mass- and charge-balanced. Next, currency metabolites (compounds used to carry electrons and functional groups) were removed following an established procedure  [@Ma2003]. Finally, the network model was converted to a metabolic network graph, in which nodes denote compounds and edges denote reactions. A directed edge from compound _a_ to compound _b_ means that _a_ is a substrate in a reaction which produces _b_. The metabolic network graph representation of glycolysis is shown in Supplemental Figure S2.

Finally, all genome-level metabolic network graphs for a single acI clade were combined to generate a composite clade-level metabolic network graph for that clade. Beginning with two genomes, nodes and edges unique to the second genome are identified and appended to the network graph for the first genome, giving a composite metabolic network graph. The process is repeated for each genome belonging to the clade, until all of the network graphs have been incorporated into the composite.

### Calculation and Evaluation of Seed Compounds

Seed compounds for each composite clade-level metabolic network graph were calculated using established methods [@Borenstein2008]. Briefly, the metabolic network is decomposed into its strongly connected components (SCCs), sets of nodes where each node in the set is reachable from every other node. Seed compounds can then be found by identifying source components (components with no incoming edges) on the condensation of the original graph: each source component represents a collection of seed compounds. The process is illustrated for an artificial network in Supplemental Figure S2. Finally, all predicted seed compounds were manually evaluated to identify those which may be biologically meaningful. Examples are given in the Results section.

### Re-annotation of Peptidases and Glycoside Hydrolases

Many seed compounds were associated with reactions catalyzed by peptidases or glycoside hydrolases, and genes associated with these reactions were re-annotated. Peptidase sequences were annotated using the MEROPS batch BLAST interface using default parameters [@Rawlings2015]. Glycoside hydrolases were first annotated using dbCAN [@Yin2012] to assign these genes to glycoside hydrolase families, as defined in the Carbohydrate-Active enZYmes Database CAZY [@Lombard2014]. Hidden Markov Models for these sub-families were then downloaded from dbCAN, and HMMER3 [@Eddy2011] was used to assign these genes to individual sub-families using default parameters. Re-annotated peptidases and glycoside hydrolases are given in Supplementary Tables S2 and S3, respectively.

## Integrating Reverse Ecology with Metatranscriptomics

### Protein Clustering, Metatranscriptomic Mapping, and Clade-Level Gene Expression

OrthoMCL [@Li2003] was used to identify clusters of orthologous genes (COGs) in the set of acI genomes. OrthoMCL was run using default options. Then, metatranscriptomic reads were mapped to a single fasta file containing all acI genomes using BBMap (https://sourceforge.net/projects/bbmap/) with the `ambig=all` and `minid=0.85` options. With this command, reads are allowed to map equally well to multiple sites. An 85% percent identity cutoff was chosen as it minimizes the percent of reads which map to multiple sites.

Next, a custom implementation of HTSeq-Count [@Anders2014] was used to count the total number of reads which map to each (genome, gene) pairing. Any read which mapped to multiple sites within the collection of acI genomes was discarded. Using the COGs identified by OrthoMCL, the total number of reads which map to each (clade, COG) pairing was then counted. This gives a measure of gene expression for the clade-level composite genome.

For each clade, read counts were then normalized by sequencing depth and ORF length and expressed on a reads per kilobase million (RPKM) basis [@Mortazavi2008], while accounting for different ORF lengths within a COG. RPKM counts were averaged across the four metatranscriptomes to give the average expression across a 24-hour period. Finally, within each clade, the percentile rank expression for each COG was calculated.

### Re-annotation of Transporter Genes

Many highly-expressed genes were annotated as transport proteins, and these proteins were re-annotated to assign function to COGs identified by OrthoMCL. Protein sequences were BLASTed against the TCDB reference database [@Saier2014] using `blastp` and given the annotation of the best-BLAST hit. If a COG contained genes with multiple annotations, the majority annotation was selected for that COG. Examples are given in the Results section.

## Availability of Data and Materials

All genomic and metatranscriptomic sequences are available through IMG and NCBI, respectively. Reverse ecology calculations were performed using the Python package reverseEcology, written expressly for this purpose and available on the Python Package Index. A reproducible version of this manuscript is available at https://github.com/joshamilton/reverseEcologyMS.

# Results

## Phylogenetics (Figure 1)
* assemble Actino genomes from our extensive FW genome collection and other sources
* identify acI as a monophyletic group w/in phylum Actinobacteria, acI-C as a group within acI - Figure 1 shows just these genomes

## Completeness Estimates (Figure 2)
* reverse ecology tested on complete genomes, SAGs and MAGs are incomplete
* construct composite genomes at tribe- and clade-level
* tribe-level remain incomplete, choose clade-level for Analysis

## Protein Clustering and Metatranscriptomics
* total number of protein clusters
* core and accessory genome for lineage based on presence/absence w/in clades
* read recruiting - % of MG/MT reads which map to acI reference genome collections (still need to map the MGs)
* what does this say about...
  * are acI active or dormant?
  * quality of our reference genomes wrt conditions that day?

## Workflow (Figure 3)
* developed a computational pipeline for reverse ecology analysis on incomplete genomes
* annotate using KBase, no gapfilling (explain justification)
* size of models, fraction of total genes, comparison to other organisms
* metabolic models converted to network graphs and merged
* check size of SCC and reduce to largest component (compare to other studies re: size, % of nodes)
* seed compound calculations - # of seed sets, % which contain one compound, max size

## Evaluation of Potential Seed Compounds (Figure 4a)
* anticipate "noise" in results and evaluate individual compounds
* screen seed compounds to identify a subset for further investigation
  * examples of reasons to reject
    * network pruning (carbamoyl phosphate)
    * missing annotations (fatty acids)
  * examples of reasons to retain
    * homoserine auxotrophy
    * peptide degradation
* manually curate selected compounds
  * false negative - genes missing from the model - lysine auxotrophy
  * false negative - alternative pathways - threonine auxotrophy
  * true result - gene missing from model and genome - homoserine auxotrophy

## Re-annotation of Peptidases and Glycoside Hydrolases (Figure 4b)
* peptidases
  * families and sub-families
  * consistency with MAR-FISH and genomic studies
  * similarities / differences btw clades - presence/absence and expression
* glycoside hydrolases
  * families and sub-families
  * consistency with MAR-FISH and genomic studies
  * similarities / differences btw clades - presence/absence and expression

## Transporters
* summary of transporters
* diverse array of N-rich compounds consistent w/ literature
* transport of simple and complex sugars
* nucleotides and vitamins
  * discuss Sarahi's vitamin auxotrophy idea?
* other transporters not shown in figure
* ecological niche
  * scavenger of cell lysate, exudates/partially degraded plant materials
  * ecological successful b/c ready to eat whatever comes its way
* actinorhodopsin also highly expressed - hypothesis: establishes a protein gradient to drive ATP transport

## Actinorhodopsin (possible Figure 5 or new analysis)
* among most highly-expressed genes
* coupled with high expression of sugar transporters confirms acI are photoheterotrophs
* retinal biosynthesis pathway expressed (and functional in acI-A and acI-B)
* other reasons for such high expression (MMBR review, need to investigate in our genomes)
  * promoter strength?
  * co-localized w/ other highly-expressed genes?

# Discussion

# Acknowledgements

# Conflict of Interest

The authors declare no conflict of interest.

# References
