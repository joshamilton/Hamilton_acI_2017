---
title: Metabolic Network Analysis and Metatranscriptomics of a Cosmopolitan and Streamlined Freshwater Lineage
author: Joshua J. Hamilton^1,\*^, TBD, and Katherine D. McMahon^1,2,\*^
abstract: \* Correspondence&#58; Joshua J. Hamilton, jjhamilton2@wisc.edu and Katherine D. McMahon, trina.mcmahon@wisc.edu
date: ^1^ Department of Bacteriology, University of Wisconsin-Madison, Madison, WI, USA; ^2^ Department of Civil and Environmental Engineering, University of Wisconsin-Madison, Madison, WI, USA
---

# Abstract

# Introduction

Microbes are critical players in freshwater systems, where they support essential ecosystem functions such as nutrient cycling. Of the freshwater bacteria, Actinobacteria of the lineage acI are among the most abundant, constituting upwards of 50% of the total bacteria in a variety of aquatic systems [@Zwart1998, @Glockner2000, @Zwart2002]. Despite their abundance, no isolates of the acI lineage have been stably propagated in pure culture.

Nevertheless, the acI lineage has been extensively studied in a community context using both DNA sequencing and single-cell targeted experiments. Most fundamentally, the acI have been phylogenetically divided into three clades (A, B, and C) and thirteen tribes [@Newton2011a] on the basis of their 16S rRNA gene sequences. Several studies have used fluorescent in situ hybridization (FISH) and catalyzed reporter deposition (CARD) or microautoradiography (MAR) to identify substrate uptake capabilities of the acI. These studies reveal that the acI are capable of consuming amino acids generally [@Salcher2010, @Salcher2013]; the individual amino acids arginine, glutamate, glutamine, and leucine [@Buck2009, @Perez2010, @Salcher2010, @Eckert2012, @Salcher2013]); the saccharides glucose [@Buck2009, @Salcher2013], and N-acetylglucosamine (NAG) [@Beier2011, @Eckert2012, @Eckert2013; the deoxynucleoside thymidine [@Perez2010, @Salcher2013], and acetate [@Buck2009]. However, due to limited phylogenetic resolution of some FISH probes, the studies cannot always link the uptake of these substrates to clades or tribes within the lineage.

Instead, metabolic reconstructions of single-cell genomes (SAGs) and metagenome-assembled genomes (MAGs) have been used to propose substrate uptake capabilities of clades aI-A and acI-B. Notably, to date, no genomes from clade acI-C have been included in these studies. These studies indicate both tribes are capable of consuming a wide array of N-containing compounds, including ammonium, branched-chain amino acids, polyamines, di-peptides, and cyanophycin [@Ghylin2014, @Garcia2015], with clade acI-A also capable of consuming oligopeptides [@Ghylin2014]. The lineage is also capable of consuming numerous saccharides, including the five-carbon sugars xylose, ribose, arabinose [@Garcia2013, @Ghylin2014, @Garcia2015] as well as poly- and oligo-saccharides [@Ghylin2014, @Garcia2015]. Notably, transporters for glucose and NAG have not yet been identified [@Garcia2013, @Ghylin2014], despite FISH studies showing uptake of those substrates. Clade acI-B is also predicted to consume sucrose and maltose [@Garcia2015]; it also contains a chitinase for the breakdown of NAG [@Garcia2013, @Garcia2015]. Finally, the acI are predicted to contain the actinobacterial opsin protein actinorhodopsin [@Garcia2013, @Garcia2014, @Ghylin2014 @Garcia2015], a light-harvesting transmembrane protein [@Sharma2008, @Sharma2009], as well as the complete pathway for the biosynthesis of its cofactor retinal [@Ghylin2014].

Finally, a recent study has predicted a number of auxotrophies in a metagenome-assembled genome, including for the amino acids isoleucine, leucine, valine, tyrosine, tryptophan, phenylalanine, asparagine; and the cofactors biotin (Vitamin B7), cobalamin (Vitamin B12), folate (Vitamin B9), niacin (Vitamin B3), pantothenate (Vitamin B5), and riboflavin (Vitamin B2) [@Garcia2015].

These metabolic reconstructions all assume that an organism's genome content says something about its ecology, an assumption underlying the emerging paradigm of reverse ecology [@Levy2012]. While metabolic reconstructions represent a common entry point to reverse ecological analyses, other approaches take cues from systems biology, focusing not just on the "parts" (genes) encoded in the genome, but on the way those parts come together and interact [@Levy2012]. In particular, this approach to reverse ecology analyzes genomes in terms of their metabolic networks with a focus on their topological properties. One such property is an organism's "seed set," the set of compounds that an organism cannot synthesize and must exogenously acquire from its environment [@Borenstein2008]. As such, these compounds may represent both auxotrophies, essential metabolites for which biosynthetic routes are missing, and nutrients, for which routes for degradation (not synthesis) are present in the genome.

In this work, we expand existing genome-based analyzes of the acI lineage trough such a reverse ecological approach. We re-analyze previously-described acI genomes, as well as six additional single-cell genomes (SAGs) and 15 metagenome-assembled genomes (MAGs), including for the first time genomes from the acI-C clade. For this analysis, we have developed a high-throughput reverse ecological analysis to predict seed compounds for each clade, using metabolic network reconstructions generated from KBase (http://kbase.us). Predicted seed compounds confirm the ability of the acI to metabolize N-rich organic compounds and an array of carbohydrates, while also revealing clade-specific differences in auxotrophies and degradation capabilities. We also present the first metatranscriptomic analysis of gene expression across the three acI clades. These data show that the acI express a diverse array of transporters at a high level, which may be key to the success of their passive and planktonic lifestyle. Finally, we observe actinorhodopsin to be among the most highly expressed genes in all three lineages, strongly suggesting a photoheterotrophic lifestyle for the acI.

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

Raw paired-end reads were then trimmed using Sickle [@Joshi2011] and merged using FLASH [@Magoc2011]. Sickle was run using default parameters, and FLASH was run with a maximum overlap of 100 (`M = 100`). Finally, additional rRNA and ncRNA sequences were removed using SortMeRNA [@Kopylova2012] using default parameters. SortMeRNA was run using eight built-in databases for bacterial, archaeal, and eukaryotic small and large ribosomal subunits and ncRNAs, derived from the SILVA 119 [@Quast2013] and RFAM [@Nawrocki2015] databases.

Additional information, including all protocols and scripts for RNA analysis, can be found on Github (https://github.com/McMahonLab/OMD-TOILv2). Raw RNA sequences can be found on the National Center for Biotechnology Information (NCBI) website under BioProject PRJNA######.

## Identification of acI SAGs and MAGs

Novel acI SAGs were identified using partial 16S rRNA genes and a reference taxonomy for freshwater bacteria, as described above. To identify acI MAGs, a phylogenetic tree containing all acI SAGs and Actinobacterial MAGs was constructed as described below. Actinobacterial MAGs forming a monophyletic group with the acI SAGs were deemed acI MAGs.

## Genome Completeness and Phylogenetic Relationships

CheckM [@Parks2015] was used to estimate genome completeness based on 204 single-copy marker genes conserved across the phylum Actinobacteria. Phylogenetic analysis of Actinobacterial SAGs and MAGs was performed using a concatenated alignment of single-copy marker genes obtained via Phylosift [@Darling2014]. Maximum likelihood trees were generated using RAxML [@Stamatakis2014] using the automatic protein model assignment option (PROTGAMMAAUTO) and 100 bootstraps.

## Metabolic Network Reconstruction and Reverse Ecology

### Genome Annotation and Model Processing

Genome annotations and metabolic network reconstructions were performed using KBase (http://kbase.us/). Unannotated contigs for each genome were pushed to KBase and annotated using the "Annotate Microbial Contigs" method using default options, which uses components of the RAST toolkit [@Brettin2015, @Overbeek2014] for genome annotation. Genome-scale metabolic network reconstructions were performed using the "Build Metabolic Model" app using default parameters, which relies on the Model SEED framework [@Henry2010a] to build a draft metabolic model.

Metabolic models were then downloaded via KBase, pruned, and converted to metabolic network graphs. In particular, biomass, exchange, transport, spontaneous, and DNA/RNA biosynthesis reactions were removed from the model, and reactions were mass- and charge-balanced. Next, currency metabolites (compounds used to carry electrons and functional groups) and highly-connected compounds (those which participate in many reactions, such as CO2 and O2) were removed following an established procedure  [@Ma2003]. Finally, the network model was converted to a metabolic network graph, in which nodes denote compounds and edges denote reactions. A directed edge from compound _a_ to compound _b_ means that _a_ is a substrate in a reaction which produces _b_.

Finally, all genome-level metabolic network graphs for a single acI clade were combined to generate a composite clade-level metabolic network graph for that clade. Beginning with two genomes, nodes and edges unique to the second genome are identified and appended to the network graph for the first genome, giving a composite metabolic network graph. The process is repeated for each genome belonging to the clade, until all of the network graphs have been incorporated into the composite.

### Calculation and Evaluation of Seed Compounds

Seed compounds for each composite clade-level metabolic network graph were calculated using established methods [@Borenstein2008]. Briefly, the metabolic network is decomposed into its strongly connected components (SCCs), sets of nodes where each node in the set is reachable from every other node. Seed compounds can then be found by identifying source components (components with no incoming edges) on the condensation of the original graph: each source component represents a seed set, and the nodes within that component represent seed compounds. Finally, all predicted seed compounds were manually evaluated to identify those which may be biologically meaningful. Examples are given in the Results section.

### Re-annotation of Peptidases and Glycoside Hydrolases

Many seed compounds were associated with reactions catalyzed by peptidases or glycoside hydrolases, and genes associated with these reactions were re-annotated. Peptidase sequences were annotated using the MEROPS batch BLAST interface using default parameters [@Rawlings2015]. Glycoside hydrolases were first annotated using dbCAN [@Yin2012] to assign these genes to glycoside hydrolase families, as defined in the Carbohydrate-Active enZYmes Database CAZY [@Lombard2014]. Hidden Markov Models for these sub-families were then downloaded from dbCAN, and HMMER3 [@Eddy2011] was used to assign these genes to individual sub-families using default parameters.

## Integrating Reverse Ecology with Metatranscriptomics

### Protein Clustering, Metatranscriptomic Mapping, and Clade-Level Gene Expression

OrthoMCL [@Li2003] was used to identify clusters of orthologous genes (COGs) in the set of acI genomes. OrthoMCL was run using default options. Annotations were assigned to protein clusters by choosing the most common annotation among all genes assigned to that cluster. Then, metatranscriptomic reads were mapped to a single fasta file containing all acI genomes using BBMap (https://sourceforge.net/projects/bbmap/) with the `ambig=random` and `minid=0.95` options. A 95% identity cutoff was chosen as it approximates the average nucleotide identity of acI sequence-discrete populations [@Garcia2016a], while competitive mapping ensures that reads map only to a single genome.

Next, HTSeq-Count [@Anders2014] was used to count the total number of reads which map to each (genome, gene) pairing. After mapping, (genome, gene) pairs were filtered to remove any pairs which did not map at least one read in all four samples. Using the COGs identified by OrthoMCL, the (genome, gene) pairs which map to each (clade, COG) pairing were then identified.

For each (clade, COG) pairing, gene expression was computed on a reads per kilobase million (RPKM) basis [@Mortazavi2008], while accounting for different sequencing depths across genomes (based on total numbers of mapped reads) and ORF lengths within a COG. Because low abundance (genome, gene) pairs were discarded after mapping, this measure of gene expression provides an underestimate of the true expression level. RPKM counts were then averaged across the four metatranscriptomes. Finally, within each clade, the percentile rank expression for each COG was calculated.

### Identification of Transporter Genes

Many highly-expressed COGs were annotated as transport proteins. We used the metabolic network reconstructions for the acI genomes to systematically characterize the transport capabilities of the acI lineage. For each genome, we identified all transport reactions present in its metabolic network reconstruction. Gene-protein-reaction associations (GPRs) for these reactions were manually curated to remove unannotated proteins, group genes into operons (if applicable), and to identify missing subunits for multi-subunit transporters. These genes were then mapped to their corresponding COGs, and GPRs were grouped on the basis of their mapped COGs. Finally, consensus annotations within each clade were used to identify likely substrates for each of these groups.

## Availability of Data and Materials

All genomic and metatranscriptomic sequences are available through IMG and NCBI, respectively. Reverse ecology calculations were performed using the Python package reverseEcology, written expressly for this purpose and available on the Python Package Index. A reproducible version of this manuscript is available at https://github.com/joshamilton/reverseEcologyMS.

# Results

## Genome Statistics and Phylogenetic Affiliation

We have assembled a reference genome collection containing 17 SAGs and 19 MAGs from members of the acI lineage. The SAGs, 11 of which have been previously described [@Garcia2013, @Ghylin2014], come from four temperate lakes in the United States and Europe, while the MAGs come from two temperate lakes in the United States (15 MAGs, nine of which have been previously-described [@Bendall2016]), Spanish and American reservoirs (three MAGs [@Ghai2014, @Tsementzi2014]), and a mixed culture from a European temperate lake [@Garcia2015]. The full list of genomes is given in Table 1.

A phylogenetic tree of these genomes is shown in Figure 1. The acI lineage has previously been shown to contain three distinct clades [@Newton2011a], and our concatenated gene tree recapitulates this topology. Of note, three MAGs were classified as belonging to the acI-C clade, and represent the first genomes from this group. Additionally, five MAGs fell into one of the seven tribes defined by our SAGs.

Genome completeness estimates for the new genomes range from 51 to 87% (Table 1), with estimated genome sizes between 1 and 2 MB. The GC content of these genomes was also low (40 to 50%), and both estimated genome size and GC content are consistent with other acI genomes [@Ghai2012, @Garcia2013, @Ghylin2014, @Garcia2015, @Tsementzi2014, @Bendall2016]. Estimated genome size and GC content of clade acI-C were not statistically different from clades acI-A and acI-B.

## Estimated Completeness of Tribe- and Clade-Level Composite Genomes

Metabolic network reconstructions created from these genomes will necessarily be missing reactions, as the underlying genomes are incomplete. Previous studies have shown that the percentage of correctly identified seed compounds (true positives) is approximately equal to the completeness of the reaction network [@Borenstein2008], and the number of false positives is approximately equal to the incompleteness of the network [@Borenstein2008].

Using conserved single-copy marker genes [@Parks2015], We estimated the completeness of tribe- and clade-level composite genomes to determine the finest level of taxonomic resolution at which we could confidently compute seed compounds (Figure 2). At the tribe level, with the exception of tribe acI-B1, tribe-level composite genomes are estimated to be incomplete (Figure 2A). At the clade level, clades acI-A and B are estimated to be complete, while acI-C remains incomplete (Figure 2B). As a result, seed compounds were calculated for composite clade-level genomes, with the understanding that some true seed compounds for the acI-C clade will not be predicted.

## Metatranscriptomics and Protein Clustering
Sequencing of cDNA from all four samples yielded approximately 160 billion paired-end reads. After merging, filtering, and _in-silico_ rRNA removal, approximately 81 billion, or 51% of the reads remained (Table S1). These reads were subsequently mapped against our collection of acI SAGs and MAGs. We used the metatranscriptomic reads that mapped to each clade as a proxies for relative activity (Table S2). Overall, our acI genomes accounted for 1.23% of the total activity.

OrthoMCL identified a total of 5013 protein clusters across the three clades (Table S3). Of these, 1078 represent core genes, defined as present in at least one genome belonging to that clade. The COGs were unequally distributed across the three clades, with clade acI-A genomes containing 3175 COGs (63%), clade acI-B genomes containing 3459 COGs (69%), and clade acI-C genomes containing 1365 COGs (27%). Of these COGs, 650 were expressed in clade acI-A, 785 in clade acI-B, and 849 in clade acI-C (Table S4). These COGs account for 0.15% (acI-A), 0.14% (acI-B), and 0.31% (acI-C), of the total activity. Within the acI, the remaining unaccounted for activity comes from non-protein encoding RNA, as we only identified COGs for protein-encoding RNA.

## A Workflow for High-Throughput Reverse Ecological Analysis of Metabolic Networks

A central contribution of this work is a computational pipeline to compute an organism's seed compounds from a graph-based representation of its metabolic network. To recap, unannotated contigs are converted to metabolic network reconstructions using KBase. The reconstructions are then converted to metabolic network graphs and combined to give composite metabolic network graphs for each clade. Seed compounds are then computed for each clade, using its composite metabolic network graph (Figure 3, and Figures S1 and S2).

Metabolic network reconstructions for the acI genomes contained between 110 and 339 genes, encoding between 241 and 587 reactions which interconvert between 374 and 699 metabolites (Table S5). On average, these genes account for 25% of the genes in the genome, a value consistent with metabolic network reconstructions for other organisms. Clade-level composite metabolic network graphs were considerably larger, with between 602 and 811 metabolites (Table S6).

These composite metabolic network graphs contained a large number of disconnected components (groups of nodes that are not connected to the bulk of the network, Figure S3). For simplicity, these components were dropped from the graph, and seed compounds were computed for the single largest component. In all cases, the single largest component contained at least 80% of the nodes in the original graph.

Decomposition of composite metabolic network graphs into their SCCs resulted in a bow-tie structure, in which a single giant component contains a substantial fraction of the compounds (Figure S3). Across the three clades, the giant component contained 61% of the metabolites, a substantially larger fraction than reported for other organisms [@Ma2003a].

The total number of predicted seed sets (source components in the SCC decomposition) ranged from 63 to 95, and the number of seed compounds ranged from 70 to 102 (Table S6). This discrepancy arises because some seed sets contain multiple compounds (an example is discussed below). However, such seed sets were rare (4% of all seed sets), and contained at most six compounds. The majority of seed compounds (96%) belonged to seed sets containing only a single compound (Table S7). A total of 125 unique seed compounds were identified across the three clades, and a complete list can be found in Table S8.

## Evaluation of Potential Seed Compounds
Seed compounds were predicted using the results of an automated annotation pipeline, and as such are likely to contain inaccuracies [@Richardson2013]. As a result, we screened the set of predicted seed compounds to identify those which represented biologically plausible auxotrophies and degradation capabilities. This subset of seed compounds were then manually curated. Tables S9 and S10 contain the final set of proposed auxotrophies and degradation capabilities, respectively, for clades acI-A, B, and C. Here, we present a series of brief vignettes explaining why compounds were retained or discarded discarded as biologically plausible. For biologically plausible compounds, we also provide examples of manual curation efforts.

__Carbamoyl phosphate__. Carbamoyl phosphate was predicted as a seed compound for all three clades. Carbamoyl phosphate synthase is the first step in arginine and pyrimidine biosynthesis, and catalyzes the reaction:

    2 ATP + L-glutamine + hydrogen carbonate + H2O → L-glutamate + carbamoyl phosphate + 2 ADP + phosphate + 2 H+

This reaction contains a number of currency metabolites (ATP, ADP glutamine, glutamate), as well as the highly-connected metabolites carbonate, water, phosphate and protons. All of these metabolites were removed from the network during pruning. Thus, the reaction responsible for carbamoyl phosphate synthesis was effectively removed from the network, rendering this compound a seed compound. Manual inspection of individual genomes revealed the gene for carbamoylphosphate synthase, confirming carbamoyl phosphate is not an auxotrophy.

__R-enoyl-ACP__. A number of R-enoyl-ACP compounds were predicted to be seed compounds in clades acI-A and B. These compounds were associated with a single COG annotated as an "Enoyl-[acyl-carrier-protein] reductase," the enzyme which catalyzes the final step in fatty acid elongation. Many other predicted seed compounds were predicted to partcipate in fatty acid and phospholipid biosynthesis, including a number of saturated fatty acids (associated with a COG annotated as a "long-chain-fatty-acid--CoA ligase") and 1-acyl-sn-glycerol 3-phosphate compounds (associated with a COG annotated as an "1-acyl-sn-glycerol-3-phosphate acyltransferase"). Given the broad substrate specificity of these enzymes, KBase automatically assigns these enzymes to the catalysis of a number of reactions. It is highly likely that manual curation will reveal complete fatty acid and phospholipid biosynthesis pathways, so and other such compounds were excluded from further consideration.

__L-Aspartate-4-semialdehyde, L-homoserine, and O-Phospho-L-homoserine.__ Clade acI-C was predicted to have a seed set containing these three compounds. These three compounds can be interconverted via the following reactions:

    Homoserine dehydrogenase: L-Aspartate-4-semialdehyde <--> L-homoserine

    Homoserine kinase: L-homoserine <--> O-Phospho-L-homoserine

making the three compounds the vertices of a strongly connected component. The first reaction, homoserine dehydrogenase, is the final step in homoserine biosynthesis, so these compounds suggest an auxotrophy for homoserine. Homoserine biosynthesis proceeds via the following reactions:

    aspartate kinase: aspartate --> L-aspartyl-4-phosphate
    aspartate semialdehyde dehydrogenase: L-aspartyl-4-phosphate --> L-Aspartate-4-semialdehyde
    homoserine dehydrogenase: L-Aspartate-4-semialdehyde --> homoserine

The presence of L-Aspartate-4-semialdehyde as a seed compound suggests the reaction "aspartate semialdehyde dehydrogenase" is missing, and were unable to manually identify a candidate gene for this reaction.  Furthermore, the acI-C composite genome contains the other two reactions in the pathway. Thus, on the evidence available, we conclude acI-C is auxotrophic for homoserine.

__L-arogenate.__ This compound was predicted as a seed compound for clade acI-C, suggesting an auxotrophy for tyrosine. Tyrosine can be synthesized via the following route:

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

These reactions were associated with four COGs, annotated as aminopeptidases. These seed compounds suggest the ability for the acI to degrade peptides, but the broad specificity of aminopeptidases means these particular di-peptides are unlikely to be the only substrates. Similarly, a number of seed compounds were associated with COGs annotated as gluco- and galactosidases, which also have broad substrate specificity. As a result, these genes were further investigated as described below.

## Auxotrophies and Degradation Capabilities of the acI Lineage

Figure 4a summarizes predicted auxotrophies for the acI lineage. In all three clades, beta-alanine was identified as a seed compound, suggesting an auxotrophy for Vitamin B5. (Vitamin B5, also known as pantothenic acid, is a precursor to coenzyme A formed from beta-alanine and pantoate). In bacteria, beta-alanine is typically synthesized via the decarboxylation of aspartate, and we were unable to identify a candidate gene for this enzyme in any acI genome (Table S9). Pyridoxine phosphate and pyridoxamine phosphate (forms of the enzyme cofactor Vitamin B6) were also predicted to be seed compounds, and numerous enzymes in the biosynthesis of these compounds were undetected in the genomes (Table S9).

Clades within the acI lineage also exhibited distinct auxotrophies. Clade acI-A was predicted to be auxotrophic for the cofactor tetrahydrofolate (THF), and numerous enzymes for its biosynthesis were missing (Table S9). In turn, acI-C was predicted to be auxotrophic for UMP and the amino acids lysine and homoserine, and in all cases multiple enzymes for the biosynthesis of these compounds went not found in the acI-C genomes. However, because the acI-C composite genome was estimated to be around 80% complete, we cannot rule out the possibility that the missing genes might be found in additional genomes.

Furthermore, both clades acI-A and B were predicted to degrade D-altronate and trans-hydroxy proline, and acI-B was additionally predicted to degrade glycine betaine. __These compounds relate to the ecology of acI...__

Finally, all three clades were predicted to degrade the di-peptides ala-leu and gly-pro-L and the sugar maltose. Clades acI-A and acI-C were also predicted to degrade the polysaccharides stachyose, manninotriose, and cellobiose. In all cases, these compounds were associated with reactions catalyzed by peptidases or glycoside hydrolases, and genes associated with these reactions were re-annotated as described above. In most cases, these annotations were in agreement with annotations given by KBase (Tables S11 and S12). The results of this re-annotation are shown in Figure 4b.

All three clades were predicted to contain both cytosolic- and membrane-bound aminopeptidases capable of releasing a variety of residues from both di- and polypeptides. As discussed below, we identified a number of transport proteins capable of transporting these released residues. The genes for these two enzymes were moderately expressed, being near the 50th percentile for gene expression in all three clades, with log2 RPKM values between 9 and 10.

All three clades were predicted to encode an alpha-glucosidase, which was expressed most strongly in clade acI-C with an log2 RPKM of 10. Clades acI-A and C also encode an additional alpha-glucosidase and an alpha-amylase, though only the alpha-amylase was expressed, and only in clade acI-C. Both of these enzymes release glucose monomers, which acI is known to consume [@Buck2009, @Salcher2013]. Furthermore, these two clades encode an alpha-galactosidase and an enzyme which could be a beta-glucosidase, beta-galactosidase, or a beta-D-fucosidase, though only the alpha-galactosidase was expressed, and only in clade acI-C.

In the aggregate, these results suggest the acI lineage is capable of degrading a diverse array of peptides and polysaccharides, such as putrescine, spermidine, and __some polysaccharides that are abundant in lakes__, which are all known to be abundant in freshwater ecosystems [REFs].

## Compounds Transported by the acI Lineage

All acI clades encode for and express a diverse array of transporters (Figure 5 and Tables S13 and S14). Consistent with the presence of intra- and extra-cellular peptidases, all clades contain numerous genes for the transport of peptides and amino acids, including multiple oligopeptide and  branched-chain amino acid transporters, as well as two distinct transporters for the polyamines spermidine and putrescine. All clades also contain a transporter for ammonium. Of these, the ammonium, branched-chain amino acid, and oligopeptide transporters are among the most highly expressed in these genomes, often above the 75th percentile of all expressed genes. In contrast, while all clades express some genes from the polyamine transporters, only clade acI-B expressed the spermidime/putrescine binding protein. Additionally, clade acI-A contains a third distinct branched-chain amino acid transporter, composed of COGs not found in clades acI-B and C. This transporter is not as highly-expressed as the shared transporters. Finally, clades acI-A and B also contain a transporter for glycine betaine, which is only expressed in clade acI-A.

All clades also strongly express transporters consistent with the presence of glycoside hydrolases, including transporters for the sugars maltose (a dimer of glucose) and xylose (an aldopentose). Clades acI-A and B also contain four distinct transporters for ribose (another aldopentose), although the substrate-binding subunit is not expressed.

The acI lineage also encodes for and expresses a number of transporters which do not have corresponding seed compounds, including a uracil permease, and a xanthine/uracil/thiamine/ascorbate family permease, both of which are highly expressed. Clades acI-A and B also contain a a cytosine/purine/uracil/thiamine/allantoin family permease, though only clade acI-B expresses it. All three clades both contain and strongly express the high-affinity phosphate specific transport system (Pst). In addition, clade acI-A contains but does not express a transporter for Vitamin B12 (cobalamin), and both clades acI-A and B contain but do not express transporters for Vitamins B1 (thiamin) and B7 (biotin). Oddly, despite predicted auxotrophies for Vitamins B5 and B6, we were unable to find transporters for these two compounds.

Finally, all three clades express actinorhodopsin, a light-sensitive opsin protein which functions as an outward proton pump [@Sharma2008]. In all clades, actinorhodopsin is among the top seven most highly-expressed genes in that clade (Table S4). Given that many of the transport proteins are of the ABC type, we speculate that actinorhodopsin may facilitate maintenance of the proton gradient necessary for ATP synthesis. Coupled with high expression levels of the diverse diverse transporters expressed by acI, this result strongly suggests that acI are photoheterotrophs.

The presence of multiple branched-chain amino acid and oligopeptide transporters attests to the importance of these compounds to acI's lifestyle. We identified a total of six distinct branched-chain amino acid transporters within our 36 freshwater acI genomes (Table S13). Five of these contain the same four COGs (LivF, LivG, LivH, and LivM), differing only in the fifth, the amino acid binding subunit. Of these five distinct amino acid binding proteins, only one is found in all three clades, with the others being found in just one or two clades. Similarly, we identified a total of ten distinct oligopeptide transporters (Table S13), each with a unique oligopeptide-binding protein (OppA). Six are found in all three clades, while the remaining four are present in just one or two clades. The diversity of these transporters both within and between clades suggests the acI can take up an array of branched-chain amino acids and oligopeptides under diverse conditions, with different clades capable of uptake under distinct conditions.

# Discussion

Our predictions of substrate utilization capabilities of the acI lineage are largely congruent with previous studies. We predict that the consumption of N-rich compounds is a universal feature of the acI lineage, with all three clades predicted to consume ammonium, branched-chain amino acids (leucine, isoleucine, and valine), the polyamines spermidine and putrescine, and oligopeptides. However, despite experimental observations of arginine, glutamate, and glutamine uptake, we failed to identify a transporter for these compounds.

Furthermore, we confirm the ability of all three clades to consume the five-carbon sugar xylose, as well as the six-carbon sugar maltose (previously thought to be restricted to clade acI-B). The acI-C genomes examined in this study did not contain transporters for ribose, suggesting that the utilization of this five-carbon sugar may be restricted to clades acI-A and B. However, we failed to identify transporters for the saccharides glucose, and N-acetylglucosamine, both of which have been experimentally shown to be consumed by acI bacteria. Furthermore, we identify additional hydrolases capable of acting on beta-glucosides, and alpha- and beta-galactosides in clades acI-A and acI-C, enzymes which had previously been detected only in clade acI-A.

We also identified transporters for the nucleobase uracil, as well as two permeases with broad specificity, capable of acting on both purine and pyrimidine nucleobases (cysotine, uracil, and xanthine), suggesting the acI may obtain these compounds from the environment instead of synthesizing them _de novo_. Finally, all clades within the acI contain actinorhodopsin and the complete retinal biosynthesis pathway. The exception seems to be clade acI-C, which is missing the beta-carotene cleavage enzyme which produces retinal.

Our analysis also suggests that auxotrophies for some vitamins may be universal features of the lineage, as we predict all clades to be auxotrophic for Vitamins B5 and B6 (previously, Vitamin B6 had been identified as an auxotrophy for clade acI-B). In addition, we predict transporters for Vitamins B1, B7 and B12, but do not predict auxotrophies for these vitamins. In addition, our analysis does not identify Vitamins B2, B3, B9, or B12 as auxotrophies for clade acI-B, a result which had been previously suggested.

Transport proteins for many of these metabolites were among the most highly expressed in the genomes, suggesting that the success of acI's passive lifestyle may be due to its ability to consume any substrate in its vicinity, without the need to activate expression of the necessary transport genes. The actinorhodopsin protein is highly expressed, and may facilitate synthesis of the ATP needed to drive acI's many ABC-type transporters. In the aggregate, these results indicate the acI are photoheterotrophs, making a living on a diverse array of N-rich compounds, sugars, and oligo- and poly-saccharides.

Our analysis also provides new insights into auxotrophies within the acI-C lineage, identifying tetrahydrofolate (THF) as an auxotrophy for clade acI-A, and lysine, homoserine, and UMP as auxotrophies for acI-C. THF is a derivative of folic acid (Vitamin B9), which was previously identified as an auxotrophy for clade acI-B. Additionally, clade acI-B was previously identified as auxotrophic for a number of amino acids, though lysine and homoserine were not among them. In the aggregate, these results provide additional support to the hypothesis that distributed metabolic pathways and metabolic complmentarity may be common features of freshwater bacterial communities [@Garcia2015].

While the reverse ecology pipeline presented here identified aspects of acI metabolism not observed previously, our approach does have some significant limitations. First, because our analysis relied on constructing clade-level composite genomes, we of necessity cannot say anything about the metabolism of individual tribes or interactions between them. And second, metabolic network analysis does account for biological pathway organization, and predicts organisms may synthesize compounds via non-canonical routes. An example is presented in the Supplemental Material. As a result, the pipeline may under-predict the number of auxotrophies for a genome, a limitation which may explain why we failed to predict many suspected auxotrophies for clade acI-B.

__Need to discuss why our acI genomes recruit so few reads__

# Acknowledgements

# Conflict of Interest

The authors declare no conflict of interest.

# References
