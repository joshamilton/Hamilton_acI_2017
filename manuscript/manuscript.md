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

### What we know from MAG data

## Reverse Ecology
* genome-gazing is a bottleneck, FISH studies are laborious - need an alternative
* reverse-ecology is an emerging paradigm to infer ecology from genomic information w/ a focus on systems-level signatures
* metabolism one aspect of ecology, long tradition of metabolic network analysis
* example of reverse ecology - seed compounds and environment

## Overview
In this work, we expand existing genome-based analyzes of the acI lineage by analyzing six additional single-cell genomes (SAGs) and 15 metagenome-assembled genomes (MAGs), including for the first time genomes from the acI-C clade. In lieu of manual reconstruction of these genomes, we use a high-throughput reverse ecological analysis to predict seed compounds for each clade, using metabolic network reconstructions generated from KBase (http://kbase.us). Predicted seed compounds confirm the ability of the acI to metabolize N-rich organic compounds and an array of carbohydrates, while also revealing clade-specific differences in auxotrophies and degradation capabilities. We also use metatranscriptomics to profile gene expression across the three acI clades. These data show that the acI express a diverse array of transporters at a high level, which may be key to the success of their passive and planktonic lifestyle. Finally, we observe actinorhodopsin to be among the most highly expressed genes in all three lineages, strongly suggesting a photoheterotrophic lifestyle for the acI.

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

Genomes were binned from each metagenomic co-assembly using MetaBat [@Kang2015], as described previously [@Bendall2016]. Briefly, contigs were classified into bins using tetranucleotide frequency and coverage patterns across the time-series and then manually curated, as previously described [@Bendall2016]. Genome sequences are available through IMG (https://img.jgi.doe.gov/cgi-bin/mer/main.cgi). Genomes were classified using taxonomic assignments from a set of 37 highly-conserved single-copy marker genes using Phylosift [@Darling2014], as previously described [@Bendall2016].

## Metatranscriptome Sampling and Sequencing

Four samples were collected from the top of the water column (depth <1m) from Lake Mendota (Madison, WI, USA) over a twenty-four hour period on August 20 and 21, 2015. For each sample, between 200 and 400 mL lake water was filtered onto a 0.2 μm polyethersulfone filter (Supor, Pall Corp) and stored at -80°C until extraction.

Prior to extraction, three samples were spiked with an internal standard to enable quantification of total transcript abundance, following an established protocol [@Satinsky2014a]. Briefly, a 970-nucleotide-long mRNA standard was synthesized using a a T7 RNA polymerase and the Riboprobe In Vitro Transcription System (Promega, Madison, WI), according to the manufacturer’s protocol. A fixed quantity of each standard (1.172 x 10^10 copies) was added independently to each lysis tube immediately prior to the addition of the sample filter.

Samples were subject to TRIzol-based RNA extraction (Thermo Fisher Scientific, Waltham, MA) followed by on-column DNAse digestion and RNA purification using an RNeasy Mini Kit (Qiagen, Venlo, Netherlands). RNA was then sent to the University of Wisconsin-Madison Biotechnology Center (https://www.biotech.wisc.edu) for sequencing. There, samples were prepared for sequencing using the TruSeq RNA Library Prep Kit v2 (Illumina, San Diego, CA), with the addition of a step for selective ribosomal RNA depletion using the Ribo-Zero rRNA Removal Kit (Bacteria) (Illumina). The resulting cDNA libraries were pooled in an equimolar ratio, and sequenced on an Illumina HiSeq2500.

Raw paired-end reads were then merged using FLASH [@Magoc2011] using default parameters. Finally, additional rRNA and ncRNA sequences were removed using SortMeRNA [@Kopylova2012] using default parameters. SortMeRNA was run using eight built-in databases for bacterial, archaeal, and eukaryotic small and large ribosomal subunits and ncRNAs, derived from the SILVA 119 [@Quast2013] and RFAM [@Nawrocki2015] databases.

Additional information, including all protocols and scripts for RNA analysis, can be found on Github (https://github.com/McMahonLab/OMD-TOILv2). Raw RNA sequences can be found on the National Center for Biotechnology Information (NCBI) website under BioProject PRJNA######.

## Identification of acI SAGs and MAGs

Novel acI SAGs were identified using partial 16S rRNA genes and a reference taxonomy for freshwater bacteria, as described above. To identify acI MAGs, a phylogenetic tree containing all acI SAGs and Actinobacterial MAGs was constructed as described below. Actinobacterial MAGs forming a monophyletic group with the acI SAGs were deemed acI MAGs.

## Genome Completeness and Phylogenetic Relationships

CheckM [@Parks2015] was used to estimate genome completeness based on 204 single-copy marker genes conserved across the phylum Actinobacteria. Phylogenetic analysis of Actinobacterial SAGs and MAGs was performed using a concatenated alignment of single-copy marker genes obtained via Phylosift [@Darling2014]. Maximum likelihood trees were generated using RAxML [@Stamatakis2014] using the automatic protein model assignment option (PROTGAMMAAUTO) and 100 bootstraps.

## Metabolic Network Reconstruction and Reverse Ecology

### Genome Annotation and Model Processing

Genome annotations and metabolic network reconstructions were performed using KBase (http://kbase.us/). Unannotated contigs for each genome were pushed to KBase and annotated using the "Annotate Microbial Contigs" method using default options, which uses components of the RAST toolkit [@Brettin2015, @Overbeek2014] for genome annotation. Genome-scale metabolic network reconstructions were performed using the "Build Metabolic Model" app using default parameters, which relies on the Model SEED framework [@Henry2010a] to build a draft metabolic model without gapfilling. As acI genomes are known to be small with numerous auxotrophies [@Garcia2015], we opted to forgo the gapfilling step and instead leverage multiple genomes to fill gaps, as described below.

Metabolic models were then downloaded via KBase, pruned, and converted to metabolic network graphs. In particular, biomass, exchange, transport, spontaneous, and DNA/RNA biosynthesis reactions were removed from the model, and reactions were mass- and charge-balanced. Next, currency metabolites (compounds used to carry electrons and functional groups) were removed following an established procedure  [@Ma2003]. Finally, the network model was converted to a metabolic network graph, in which nodes denote compounds and edges denote reactions. A directed edge from compound _a_ to compound _b_ means that _a_ is a substrate in a reaction which produces _b_. Supplemental Figure 2 illustrates process by which a genome gets converted to a pruned metabolic network graph, for a genome containing only glycolysis.

Finally, all genome-level metabolic network graphs for a single acI clade were combined to generate a composite clade-level metabolic network graph for that clade. Beginning with two genomes, nodes and edges unique to the second genome are identified and appended to the network graph for the first genome, giving a composite metabolic network graph. The process is repeated for each genome belonging to the clade, until all of the network graphs have been incorporated into the composite. Figure 3A shows metabolic network graphs for three acI-C genomes, and Figure 3B shows the composite metabolic network graph for clade acI-C.

### Calculation and Evaluation of Seed Compounds

Seed compounds for each composite clade-level metabolic network graph were calculated using established methods [@Borenstein2008]. Briefly, the metabolic network is decomposed into its strongly connected components (SCCs), sets of nodes where each node in the set is reachable from every other node. Seed compounds can then be found by identifying source components (components with no incoming edges) on the condensation of the original graph: each source component represents a collection of seed compounds.
Supplemental Figure 3 illustrates this process for a simple network containing only glycolysis, and Figure 3B shows seed compounds for clade acI-C. Finally, all predicted seed compounds were manually evaluated to identify those which may be biologically meaningful. Examples are given in the Results section.

### Re-annotation of Peptidases and Glycoside Hydrolases

Many seed compounds were associated with reactions catalyzed by peptidases or glycoside hydrolases, and genes associated with these reactions were re-annotated. Peptidase sequences were annotated using the MEROPS batch BLAST interface using default parameters [@Rawlings2015]. Glycoside hydrolases were first annotated using dbCAN [@Yin2012] to assign these genes to glycoside hydrolase families, as defined in the Carbohydrate-Active enZYmes Database CAZY [@Lombard2014]. Hidden Markov Models for these sub-families were then downloaded from dbCAN, and HMMER3 [@Eddy2011] was used to assign these genes to individual sub-families using default parameters. Re-annotated peptidases and glycoside hydrolases are given in Supplementary Tables S2 and S3, respectively.

## Integrating Reverse Ecology with Metatranscriptomics

### Protein Clustering, Metatranscriptomic Mapping, and Clade-Level Gene Expression

OrthoMCL [@Li2003] was used to identify clusters of orthologous genes (COGs) in the set of acI genomes. OrthoMCL was run using default options. Annotations were assigned to protein clusters by choosing the most common annotation among all genes assigned to that cluster. Then, metatranscriptomic reads were mapped to a single fasta file containing all acI genomes using BBMap (https://sourceforge.net/projects/bbmap/) with the `ambig=all` and `minid=0.85` options. With this command, reads are allowed to map equally well to multiple sites. An 85% percent identity cutoff was chosen as it minimizes the percent of reads which map to multiple sites.

Next, a custom implementation of HTSeq-Count [@Anders2014] was used to count the total number of reads which map to each (genome, gene) pairing. Any read which mapped to multiple sites within the collection of acI genomes was discarded. Using the COGs identified by OrthoMCL, the total number of reads which map to each (clade, COG) pairing was then counted. This gives a measure of gene expression for the clade-level composite genome. Because ambiguous reads were discarded after mapping, this measure of gene expression provides an underestimate of the true expression level.

For each clade, read counts were then normalized by sequencing depth and ORF length and expressed on a reads per kilobase million (RPKM) basis [@Mortazavi2008], while accounting for different ORF lengths within a COG. RPKM counts were averaged across the four metatranscriptomes to give the average expression across a 24-hour period. Finally, within each clade, the percentile rank expression for each COG was calculated. Figure 3C shows the calculated RPKM values for clade acI-C, along with the presence/absence of each COG in three acI-C genomes.

### Re-annotation of Transporter Genes

Many highly-expressed genes were annotated as transport proteins, and these proteins were re-annotated to assign function to COGs identified by OrthoMCL. Protein sequences were BLASTed against the TCDB reference database [@Saier2014] using `blastp` and given the annotation of the best-BLAST hit. If a COG contained genes with multiple annotations, the majority annotation was selected for that COG. Examples are given in the Results section.

## Availability of Data and Materials

All genomic and metatranscriptomic sequences are available through IMG and NCBI, respectively. Reverse ecology calculations were performed using the Python package reverseEcology, written expressly for this purpose and available on the Python Package Index. A reproducible version of this manuscript is available at https://github.com/joshamilton/reverseEcologyMS.

# Results

## Genome Statistics and Phylogenetic Affiliation

We have assembled a reference genome collection containing 17 SAGs and 19 MAGs from members of the acI lineage. The SAGs, 11 of which have been previously described [@Garcia2013, @Ghylin2014], come from four temperate lakes in the United States and Europe, while the MAGs come from two temperate lakes in the United States (15 MAGs, nine of which have been previously-described [@Bendall2016]), Spanish and American reservoirs (three MAGs [@Ghai2014, @Tsementzi2014]), and a mixed culture from a European temperate lake [@Garcia2015]. The full list of genomes is given in Table 1.

A phylogenetic tree of these genomes is shown in Figure 1. The acI lineage has previously been shown to contain three distinct clades [@Newton2011a], and our concatenated gene tree recapitulates this topology. Of note, three MAGs were classified as belonging to the acI-C clade, and represent the first genomes from this group. Additionally, five MAGs fell into one of the seven tribes defined by our SAGs.

Genome completeness estimates for the new genomes range from 51 to 87% (Table 1), with estimated genome sizes between 1 and 2 MB. The GC content of these genomes was also low (40 to 50%), and both estimated genome size and GC content are consistent with other acI genomes. Estimated genome size and GC content of clade acI-C were not statistically different from clades acI-A and acI-B.

## Estimated Completeness of Tribe- and Clade-Level Composite Genomes

Metabolic network reconstructions created from these genomes will necessarily be missing reactions, as the underlying genomes are incomplete. Previous studies have shown that the percentage of correctly identified seed compounds (true positives) is approximately equal to the completeness of the reaction network [@Borenstein2008], and the number of false positives is approximately equal to the incompleteness of the network [@Borenstein2008].

Using conserved single-copy marker genes [@Parks2015], We estimated the completeness of tribe- and clade-level composite genomes to determine the finest level of taxonomic resolution at which we could confidently compute seed sets (Figure 2). At the tribe level, with the exception of tribe acI-B1, tribe-level composite genomes are estimated to be incomplete (Figure 2A). At the clade level, clades acI-A and B are estimated to be complete, while acI-C remains incomplete (Figure 2B). As a result, seed sets were calculated for composite clade-level genomes, with the understanding that some true seed compounds for the acI-C clade will not be predicted.

## Protein Clustering and Metatranscriptomics

OrthoMCL identified a total of 5013 protein clusters across the three clades (Table S8). Of these, 1078 represent core genes, defined as present in at least one genome belonging to that clade.

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
* diverse array of N-rich compounds
* transport of simple and complex sugars
* nucleotides and vitamins
* other transporters not shown in figure
* actinorhodopsin also highly expressed

## Actinorhodopsin (possible Figure 5 or new analysis)
* among most highly-expressed genes
* coupled with high expression of sugar transporters confirms acI are photoheterotrophs
* retinal biosynthesis pathway expressed (and functional in acI-A and acI-B)
* other reasons for such high expression (MMBR review, need to investigate in our genomes)
  * promoter strength?
  * co-localized w/ other highly-expressed genes?

# Discussion

## (Dis)Agreement w/ previous findings
* diverse array of N-rich compounds & transport of sugars consistent w/ literature
* cyanophycinase?

## Ecology of acI
* ecological niche
  * scavenger of cell lysate, exudates/partially degraded plant materials
  * ecological successful b/c ready to eat whatever comes its way
* actR hypothesis: establishes a protein gradient to drive ATP transport
* discuss Sarahi's vitamin auxotrophy idea?

## Potential and Limitations of Reverse Ecology
* potential applications
* limitations
  * because analysis performed at clade level, cannot capture interactions btw tribes
  * only as good as database, known functions of interest are missing
  * lacks pathway information, enabling alternative metabolic routes - may miss true auxotrophies

# Acknowledgements

# Conflict of Interest

The authors declare no conflict of interest.

# References
