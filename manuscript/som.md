---
title: Supplemental Material for 'High-Throughput Metabolic Network Analysis and Metatranscriptomics of a Cosmopolitan and Streamlined Freshwater Lineage'
author: Joshua J. Hamilton^1,\*^, TBD, and Katherine D. McMahon^1,2,\*^
abstract: \* Correspondence&#58; Joshua J. Hamilton, jjhamilton2@wisc.edu and Katherine D. McMahon, trina.mcmahon@wisc.edu
date: ^1^ Department of Bacteriology, University of Wisconsin-Madison, Madison, WI, USA; ^2^ Department of Civil and Environmental Engineering, University of Wisconsin-Madison, Madison, WI, USA
---

# Supplementary Methods

## A Freshwater Reference Genome Collection

This study relies on an extensive collection of freshwater bacterial genomes, containing metagenome-assembled genomes (MAGs) obtained from two metagenomic time-series from two Wisconsin lakes, as well as single-cell genomes (SAGs) from three lakes in the United States. Additional information about this genome collection can be found below.

### Single-Cell Genome Generation, Classification, and Sequencing

Water samples were collected from the top of the water column (depth <1m) from each of three lakes, Lake Mendota (Madison, WI, USA), Sparkling Lake (Vilas County, WI, USA), and Damariscotta Lake (Lincoln County, ME, USA), in 2009. Samples were cryopreserved and sent to the Single Cell Genomics Center at the Bigelow Laboratory for Ocean Sciences for sorting, as previously described [@Martinez-Garcia2012; @Garcia2013]. Partial 16S rRNA genes amplified previously [@Martinez-Garcia2012] were phylogenetically classified using a controlled nomenclature for freshwater bacteria [@Newton2011a] by insertion into references trees created in the ARB software package [@Ludwig2004].

Actinobacterial SAGs used in this study were then sent to the JGI for sequencing and assembly, also as previously described [@Ghylin2014]. Briefly, shotgun libraries were constructed for each of the SAGs from re-amplified MDA products and sequenced on an Illumina HiSeq2000. All general aspects of and detailed protocols for library construction and sequencing can be found on the JGI website (http://www.jgi.doe.gov/).

For assembly, raw sequence data was first passed through a filtering program developed at JGI to eliminate known sequencing and library preparation artifacts. Assembly was then performed using Velvet [@Zerbino2008] and ALLPATHS-LG [@Gnerre2011]. Additional details of the assembly process have been previously described [@Ghylin2014] and are available through the JGI Genome Portal (http://genome.jgi.doe.gov).f Genome sequences are available through IMG (https://img.jgi.doe.gov/cgi-bin/mer/main.cgi). Genome-specific information can be accessed in both databases by searching for the IMG Taxon OIDs given in Table 1.

### Metagenome Sampling, Sequencing, Assembly, and Binning

Water samples were collected from two lakes, Lake Mendota (Madison, WI, USA) and Trout Bog (Vilas County, WI, USA). For Lake Mendota, depth-integrated water samples were collected from the top 12 meters at 96 time points during ice-free periods from 2008 to 2011. For Trout Bog, depth-integrated water samples were collected from the epilimnion (44 samples) and hypolimnion (45 samples) layers during ice-free periods from 2007 to 2009.

Sample collection, DNA sequencing, and metagenomic assembly, have been described previously [@Bendall2016; @Garcia2016a], as well as genomic binning for the Trout Bog samples [@Bendall2016]. Similar genomic binning procedures were followed for the Lake Mendota samples. A summary is provided here.

All samples were filtered on 0.2 μm polyethersulfone filters (Supor, Pall Corp) prior to storage at -80°C, as described previously [@Bendall2016; @Garcia2016a]. DNA was extracted from these filters using the FastDNA kit (MP Biomedicals) and sent to the JGI for sequencing, as described previously [@Bendall2016; @Garcia2016a].

Shotgun libraries were constructed for each of the samples and sequenced on an Illumina HiSeq2000, following a 2x150 indexed run recipe as previously described [@Bendall2016; @Garcia2016a]. All general aspects of and detailed protocols for library construction and sequencing can be found on the JGI website (http://www.jgi.doe.gov/). Metagenomic sequence reads are publicly available on the JGI Genome Portal (http://genome.jgi.doe.gov/) under Proposal ID 394. The full list of metagenomes is given in Table S12.

Raw sequence data was passed through a filtering program developed at JGI to eliminate known sequencing and library preparation artifacts. Prior to assembly, reads were merged with FLASH [@Magoc2011], as previously described [@Bendall2016; @Garcia2016a]. Merged reads were pooled by lake and layer into three co-assemblies using SOAPdenovo [@Luo2012a], and contigs from the resulting assemblies were assembled into a final assembly using Minimus [@Sommer2007], as previously described [@Bendall2016; @Garcia2016a]. Additional details of the assembly process are available through the JGI Genome Portal (http://genome.jgi.doe.gov) under Proposal ID 394.  The full list of assemblies is given in Table S13.

Genomes were binned from each metagenomic co-assembly using MetaBat [@Kang2015], as described previously [@Bendall2016]. Briefly, contigs were classified into bins using tetranucleotide frequency and coverage patterns across the time-series and then manually curated, as previously described [@Bendall2016]. Genome sequences are available through IMG (https://img.jgi.doe.gov/cgi-bin/mer/main.cgi), and can be accessed by searching for the IMG Taxon OIDs given in Table 1.

Genomes were classified using taxonomic assignments from a set of 37 highly-conserved single-copy marker genes using Phylosift [@Darling2014], as previously described [@Bendall2016]. CheckM [@Parks2015] was used to estimate genome completeness based on 204 single-copy marker genes conserved across the phylum Actinobacteria.

## Metatranscriptome Sampling and Sequencing

This study used four metatranscriptomes obtained as part of a 24-hour sampling experiment designed to identify diel trends in freshwater microbial communities. Samples were collected from the top of the water column (depth <1m) from Lake Mendota (Madison, WI, USA) on August 20 and 21, 2015. For each sample, between 200 and 400 mL lake water was filtered onto a 0.2 μm polyethersulfone filter (Supor, Pall Corp), flash frozen in liquid nitrogen, and stored at -80°C until extraction.

Samples were subject to TRIzol-based RNA extraction (Thermo Fisher Scientific, Waltham, MA) followed by phenol-chloroform separation and RNA precipitation. RNA was purified following an on-column DNAse digestion using the RNase-Free DNase Set (Qiagen, Venlo, Netherlands) and cleaned up with the RNeasy Mini Kit (Qiagen, Venlo, Netherlands). RNA was then sent to the University of Wisconsin-Madison Biotechnology Center (https://www.biotech.wisc.edu) for sequencing. There, samples were prepared for sequencing using the TruSeq RNA Library Prep Kit v2 (Illumina, San Diego, CA), with a ribosomal RNA (rRNA) depletion step using the Ribo-Zero rRNA Removal Kit (Bacteria) (Illumina). The resulting cDNA libraries were pooled in an equimolar ratio, and sequenced on an Illumina HiSeq2500 platform.

Raw paired-end reads were then trimmed using Sickle [@Joshi2011] and merged using FLASH [@Magoc2011]. Sickle was run using default parameters, and FLASH was run with a maximum overlap of 100 nucleotides (`M = 100`). Finally, additional rRNA and ncRNA sequences were removed using SortMeRNA [@Kopylova2012] using default parameters. SortMeRNA was run using eight built-in databases for bacterial, archaeal, and eukaryotic small and large ribosomal subunits and ncRNAs, derived from the SILVA 119 [@Quast2013] and RFAM [@Nawrocki2015] databases.

Metadata about the four samples used in this study can be found in Table S1, and the raw RNA sequences can be found on the National Center for Biotechnology Information (NCBI) website under BioProject PRJNA362825. Additional information, including all protocols and scripts for sample collection, RNA extraction, sequencing, and bioinformatic analysis can be found on Github (https://github.com/McMahonLab/OMD-TOILv2, DOI:######).

## Genome Annotation, Metabolic Network Reconstruction, and Computation and Evaluation of Seed Compounds

Genome annotations were performed and metabolic network reconstructions were built using KBase. Contigs for each genome were uploaded to KBase and annotated using the "Annotate Microbial Contigs" method with default options, which uses components of the RAST toolkit for genome annotation [@Brettin2015; @Overbeek2014]. Metabolic network reconstructions were obtained using the "Build Metabolic Model" app with default parameters, which relies on the Model SEED framework [@Henry2010a] to build a draft reconstruction.

Reconstructions were then pruned and converted to metabolic network graphs (Figure S1). During this process, exchange and transport reactions were removed from the reconstruction, to prevent extracellular metabolites from being identified as seed compounds. The biomass reaction was also removed, as KBase generates generalized biomass equations that may not reflect acI-specific biomass requirements. Finally, DNA/RNA replication reactions were removed, as these reactions do not represent metabolic processes. Reactions in the reconstructions were then mass- and charge-balanced. Next, currency metabolites (compounds used to carry electrons and functional groups) and highly-connected compounds (those that participate in many reactions, such as CO2 and O2) were removed to ensure paths in the resulting metabolic network graph would be biologically meaningful [@Ma2003]. Finally, the metabolic network graph was extracted from the reconstruction, to enable graph-theoretical identification of the network's seed set.

Many of the individual acI genomes are incomplete (see Results). Therefore, composite metabolic network graphs were constructed for each clade, to increase the accuracy of seed identification (Figure S2). To do so, all genome-level metabolic network graphs for all genomes within each acI clade were combined to generate a composite clade-level metabolic network graph. Beginning with two genomes, nodes and edges unique to the second genome are identified and appended to the network graph for the first genome, giving a composite metabolic network graph. The process is repeated for each genome, until all of the network graphs have been incorporated into the composite.

Formally, the seed set of the network is defined as the minimal set of compounds that cannot be synthesized from other compounds in the network, and whose presence enables the synthesis of all other compounds in the network [@Borenstein2008]. Seed compounds for each composite clade-level metabolic network graph were calculated using the seed set framework [@Borenstein2008] (Figure S3). Briefly, the graph is decomposed into its strongly connected components (SCCs), sets of nodes such that each node in the set is reachable from every other node. Seed compounds can then be found by identifying source components (components with no incoming edges) on the condensation of the original graph, a representation in which each SCC is represented as a single vertex. Here, each source component represents a seed set, and the nodes within that vertex represent seed compounds. If a seed set contains multiple seed compounds, each compound is assigned a weight equal to the size of the seed set. Because seed compounds are computed from a metabolic network, it is important to manually evaluate all predicted seed compounds to identify those that may be biologically meaningful, and do not arise from errors in the metabolic network reconstruction. Examples of this process are given below.

All steps were implemented using Python scripts, freely available as part of the reverseEcology Python package (https://pypi.python.org/pypi/reverseEcology/, DOI:######).

# Supplementary Results and Discussion

## Computation of Potential Seed Compounds

Metabolic network reconstructions for individual genomes contained between 110 and 339 genes, encoding between 241 and 587 reactions which interconvert between 374 and 699 metabolites (Table S14). On average, these genes account for 25% of the genes in the genome, a value consistent with metabolic network reconstructions for other organisms [see Supplementary Table 2 in @Feist2009a]. Clade-level composite metabolic network graphs were larger, with between 602 and 811 metabolites (Table S15).

These composite metabolic network graphs contained a large number of disconnected components (groups of nodes that are not connected to the bulk of the network, Figure S4). For simplicity, these components were dropped from the graph, and seed compounds were computed for the single largest component. In all cases, the single largest component contained at least 80% of the nodes in the original graph (Table S15).

Decomposition of the three metabolic network graphs into their strongly connected components (SCCs) resulted in a bow-tie structure, in which a single giant component contains a substantial fraction of the compounds (Figure S4). Across the three clades, the giant component contained 61% of the metabolites, a larger fraction than reported for other organisms [@Ma2003a], which may be a consequence of the small and streamlined genomes of acI bacteria.

The total number of predicted seed sets (source components in the SCC condensation) ranged from 63 to 95, and the number of seed compounds ranged from 70 to 102 (Table S15). This discrepancy arises because some seed sets contain multiple compounds (Table S16). However, such seed sets were rare (4% of all seed sets), and contained at most six compounds (Table S16). The majority of seed compounds (96%) belonged to seed sets containing only a single compound (Table S16).

## Evaluation of Potential Seed Compounds

Here, we present a series of brief vignettes explaining why particular compounds were retained or discarded based on their biological (im)plausibility. For biologically plausible compounds, we also provide examples of manual curation efforts. These vignettes are not intended to provide a comprehensive explanation for all compounds, but to provide a flavor of the types of decisions that went into accepting or rejecting particular compounds.

__Carbamoyl phosphate__. Carbamoyl phosphate was predicted as a seed compound for all three clades. Carbamoyl phosphate synthase is the first step in arginine and pyrimidine biosynthesis, and catalyzes the reaction:

2 ATP + L-glutamine + hydrogen carbonate + H2O → L-glutamate + carbamoyl phosphate + 2 ADP + phosphate + 2 H+

This reaction contains a number of currency metabolites (ATP, ADP glutamine, glutamate), as well as the highly-connected metabolites carbonate, water, phosphate and protons. All of these metabolites were removed from the network during pruning. Thus, the reaction responsible for carbamoyl phosphate synthesis was effectively removed from the network, rendering carbamoyl phosphate a seed compound. Manual inspection of individual genomes revealed the gene for carbamoyl phosphate synthase, confirming carbamoyl phosphate is not an auxotrophy.

__R-enoyl-ACP__. A number of R-enoyl-ACP compounds were predicted to be seed compounds in clades acI-A and acI-B. These compounds were associated with a single COG annotated as an "Enoyl-[acyl-carrier-protein] reductase," the enzyme which catalyzes the final step in fatty acid elongation. Many other seed compounds were predicted to participate in fatty acid and phospholipid biosynthesis, including a number of saturated fatty acids (associated with a COG annotated as a "long-chain-fatty-acid--CoA ligase") and 1-acyl-sn-glycerol 3-phosphate compounds (associated with a COG annotated as an "1-acyl-sn-glycerol-3-phosphate acyltransferase"). Given the broad substrate specificity of these enzymes, KBase automatically assigns these enzymes to the catalysis of a number of reactions. Fatty acid and phospholipid biosynthesis pathways are often organism-specific and unlikely to be properly annotated by automatic metabolic reconstruction pipelines. Thus, these compounds were excluded from further consideration.

__L-Aspartate-4-semialdehyde, L-homoserine, and O-Phospho-L-homoserine.__ Clade acI-C was predicted to have a seed set containing these three compounds. These three compounds can be interconverted via the following reactions:

homoserine dehydrogenase: L-aspartate-4-semialdehyde <--> L-homoserine

homoserine kinase: L-homoserine <--> O-phospho-L-homoserine

making the three compounds members of a SCC. The first reaction, homoserine dehydrogenase, is the final step in homoserine biosynthesis, so these compounds suggest an auxotrophy for homoserine. Homoserine biosynthesis proceeds via the following reactions:

aspartate kinase: aspartate --> L-aspartyl-4-phosphate

aspartate semialdehyde dehydrogenase: L-aspartyl-4-phosphate --> L-aspartate-4-semialdehyde

homoserine dehydrogenase: L-aspartate-4-semialdehyde --> homoserine

The presence of L-aspartate-4-semialdehyde as a seed compound suggests the reaction "aspartate semialdehyde dehydrogenase" is missing, and were unable to manually identify a candidate gene for this reaction.  Furthermore, the acI-C composite genome contains the other two reactions in the pathway. Thus, on the evidence available, we conclude acI-C is auxotrophic for homoserine.

__L-arogenate.__ This compound was predicted as a seed compound for clade acI-C, suggesting an auxotrophy for tyrosine. Tyrosine can be synthesized via the following reactions:

chorismate mutase: chorismate --> prephenate

prephenate aminotransferase: prephanate --> L-arogenate

arogenate dehydrogenase: L-arogenate --> L-tyrosine

L-arogenate was predicted to be a seed compound based on the presence of "arogenate dehydrogenase", the final step in the pathway. The reaction "chorismate mutase" is also present, but we were unable to find a candidate gene for the reaction "prephenate aminotransferase," suggesting an auxotrophy for tyrosine. However, L-tyrosine can be synthesized from chorismate via an alternative pathway:

chorismate mutase: chorismate --> prephenate

prephenate dehydrogenase: prephenate --> 4-hydroxyphenylpyruvate

tyrosine aminotransferase: 4-hydroxyphenylpyruvate --> L-tyrosine

All three genes in this pathway are present in the genome, indicating acI-C is not auxotrophic for tyrosine.

__Ala-Leu and gly-pro-L__. These di-peptides were predicted to be seed compounds for all three clades. The compounds are associated with the following reactions:

H2O + Ala-Leu --> L-Leucine + L-Alanine

H2O + Gly-Pro --> Glycine + L-Proline

These reactions were associated with four COGs, annotated as aminopeptidases. These seed compounds suggest the ability for the acI to degrade peptides, but the broad specificity of aminopeptidases means these particular di-peptides are unlikely to be the only substrates. Similarly, a number of seed compounds were associated with COGs annotated as gluco- and galactosidases, which also have broad substrate specificity. As a result, these genes were further investigated as described in the primary manuscript.

## Compounds Transported by the acI Lineage

The presence of multiple branched-chain amino acid and oligopeptide transporters attests to the importance of these compounds to acI's lifestyle. These ABC transporters are composed of four subunits, including two membrane-associated ATPases and two transmembrane proteins that generally determine the substrate specificity of the transporter [@Higgins1992]. We identified a total of ten distinct oligopeptide transporters within our 36 freshwater acI genomes (Table S10), each with a unique transmembrane (oligopeptide-binding) protein. Six of these transporters are found in all three clades, while the remaining four are present in just one or two clades (Table S10). Similarly, we identified a total of six distinct branched-chain amino acid transporters. In these transporters, an amino acid-binding protein, rather than the transmembrane proteins, determines the substrate specificity [@Adams1990]. Five of the six transporters in our genomes contain the same four transmembrane and ATPase subunits, differing only in the amino acid binding subunit (Table S10). Of these five distinct amino acid binding proteins, only one is found in all three clades, with the others being found in just one or two clades (Table S10). The diversity of these transporters both within and between clades suggests the acI are adapted to a variety of amino acids and oligopeptides, with further specialization within each clade.

# References
