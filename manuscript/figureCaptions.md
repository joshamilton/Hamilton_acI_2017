---
title: Figure Captions for 'Metabolic Network Analysis and Metatranscriptomics Reveals Auxotrophies and Nutrient Sources of the Cosmopolitan Freshwater Microbial Lineage acI'
author: Joshua J. Hamilton^1\*^, Sarahi L. Garcia^2, Brittany S. Brown^1, Ben O. Oyserman^3, Francisco Moya-Flores^3, Stefan Bertilsson^2,4, Rex R. Malmstrom^5, Katrina T. Forest^1, Katherine D. McMahon^1,3
abstract: \* Correspondence&#58; Joshua J. Hamilton, jjhamilton2@wisc.edu
date: ^1^ Department of Bacteriology, University of Wisconsin-Madison, Madison, WI, USA; ^2^ Department of Ecology and Genetics, Uppsala University, Uppsala, Sweden; ^3^ Department of Civil and Environmental Engineering, University of Wisconsin-Madison, Madison, WI, USA; ^4^ Science for Life Laboratory, Uppsala University, Uppsala, Sweden; ^5^ United States Department of Energy Joint Genome Institute, Walnut Creek, CA, USA
---

&nbsp;

## Figure 1 (figures/fig1-tree/Figure1.pdf)

&nbsp;

Phylogenetic placement of the genomes used in this study within the acI lineage. The tree was built using RAxML [@Stamatakis2014] from a concatenated alignment of protein sequences from 37 single-copy marker genes [@Darling2014]. The order Actinomycetales forms the outgroup. Vertical black bars indicate groups of genomes belonging to defined tribes/clades within the acI lineage, as determined using 16S rRNA gene sequences (for SAGs and bin FNEF8-2 bin_7 acI-B only) and a defined taxonomy [@Newton2011a]. SAGs are indicated with italic text. Supplemental Figure S6 shows the position of the acI lineage relative to other orders within the class Actinobacteria.

&nbsp;

## Figure 2 (figures/fig3-workflow/Figure3.pdf)

&nbsp;

Overview of the seed set framework and metatranscriptomic mapping, using three genomes from the acI-C clade as an example. __(A)__ Metabolic network graphs are created for each genome belonging to clade acI-C. In these graphs, metabolites are represented as nodes (circles) and reactions by arcs (arrows). Grey nodes and edges indicate components of the composite graph missing from that genome graph. Additional information on this step of the workflow is available in Figure S1. __(B)__ A composite network graph is created for each clade by joining graphs for all genomes from that clade, and seed compounds (red) are computed for the composite graph. Additional information on this step of the workflow is available in Figures S2, S3, and S4. __(Inset)__ Three seed compounds which indicate an auxotrophy for L-homoserine, a methionine precursor. __(C)__ Metatranscriptomic reads are mapped to each individual genome using [BBMap](https://sourceforge.net/projects/bbmap/). Orthologous gene clusters are identified using OrthoMCL [@Li2003]. For each cluster, unique reads which map to any gene within that cluster are counted using HTSeq [@Anders2014]. The relative gene expression is computed using RPKM [@Mortazavi2008].

&nbsp;

## Figure 3 (figures/fig4-seeds/Figure4.pdf)

&nbsp;

Seed compounds of members of the acI lineage. __(A)__ Auxotrophies and nutrient sources, not including peptides and glycosides. __(B)__ Peptides and glycosides. These compounds represent those inferred from genome annotations, rather than the seed compounds themselves. In panel (B), the intensity of the color indicates the log2 RPKM of the encoding gene cluster. For compounds acted upon by multiple gene clusters, the percentile of the most highly-expressed cluster was chosen.

&nbsp;

## Figure 4 (figures/fig5-transporters/Figure5.pdf)

&nbsp;

Transporters that are actively expressed by members of the acI lineage, as inferred from consensus annotations of genes associated with transport reactions present in metabolic network reconstructions. The intensity of the color indicates the log2 RPKM of the encoding gene cluster. For multi-subunit transporters, the RPKM of the substrate-binding subunit was chosen (Table S12). For some transporters, consensus annotations have been replaced with broad metabolite classes. Such metabolite classes are indicated with superscripts, and the original annotations are as follows: 1 - spermidine, putrescine; 2 - maltose ; 3 - xylose; 4 - ribose; 5 - uracil; 6 - cytosine / purine / uracil / thiamine / allantoin; 7 - xanthine / uracil / thiamine / ascorbate.

&nbsp;

## Supplementary Figure 1 (figures/fig3-workflow/FigureS1.pdf)

&nbsp;

Converting an unannotated genome to a metabolic network graph, for a simplified genome containing only glycolysis. __(A)__ Microbial contigs are annotated using KBase, and a metabolic network reconstruction is built from the annotations. The reconstruction provides links between protein-encoding genes in the genome and the enzymatic reactions catalyzed by those proteins. __(B)__ The metabolic network reconstruction represents metabolism as a hypergraph, in which metabolites are represented as nodes and reactions as hyperedges. In this representation, an edge can connect more than two nodes. For example, a single hyperedge (denoted by a heavy black line) connects the metabolites glucose and ATP to glucose-6P, ADP, and Pi. For clarity, protons are not shown. __(C)__ However, the algorithm used by the seed set framework requires metabolism to be represented as a metabolic network graph, in which an edge can connect only two nodes. In this representation, a reaction is represented by a set of edges connecting all substrates to all products. For example, the heavy hyperedge in (B) is now denoted by six separate edges connecting glucose to ADP, glucose to Pi, glucose to glucose-6P, ATP to ADP, ATP to Pi, and ATP to glucose-6P (again denoted by heavy black lines). Of these, only one (glucose to glucose-6P) is biologically meaningful. The dotted line surrounds the currency metabolites. __(D)__ The metabolic network graph is then pruned, a process which removes all currency metabolites and any edges in which those metabolites participate. Of the six heavy edges in (C), only the biologically meaningful one is retained, connecting glucose to glucose-6P (again denoted by a heavy black line). The images in (B) and (C) are modified from [@Ma2003].

&nbsp;

## Supplementary Figure 2 (figures/fig3-workflow/FigureS2.pdf)

&nbsp;

Construction of composite metabolic network graph for clade acI-C. Beginning with metabolic network graphs for genomes Actinobacterium_10 and ME00885, nodes and edges unique to ME00885 are identified (in blue). These nodes and edges are added to the Actinobacterium_10 graph, giving the composite metabolic network graph for these two genomes (Actinobacterium_10 + ME00885). Then, this graph is compared to the graph for ME03864, and nodes and edges unique to ME03864 are identified (in blue). These nodes and edges are added to the Actinobacterium_10 + ME00885 metabolic network graph, giving the composite metabolic network graph for clade acI-C.

&nbsp;

## Supplementary Figure 3 (figures/fig3-workflow/FigureS3.pdf)

&nbsp;

Identifying seed compounds in metabolic networks, using the same metabolic network as in Supplemental Figure S1. __(A)__ To identify seed compounds, the metabolic network graph is first decomposed into its strongly connected components (SCCs), sets of nodes such that each node in the set is reachable from every other node. Here, each set of circled nodes corresponds to a unique SCC. __(B)__ SCC decomposition enables seed sets to be identified from source components (components with no incoming edges) on the condensation of the original graph. In the condensation of the original graph shown here, each node corresponds to a unique SCC. This network has a single seed set, SCC_1, enclosed in a dotted circle. __(C)__ Seed compounds can be found from the mapping between SCCs and their constituent metabolites. In this example, glucose is the sole seed compound. While this particular result is probably intuitive, real metabolic networks are considerably more complex.  Note: The visual representations shown here are intended to illustrate the metabolic network reconstruction process, and are not indicative of the data structures used by our pipeline.

&nbsp;

## Supplementary Figure 4 (figures/fig3-workflow/FigureS4.pdf)

&nbsp;

Complete composite metabolic network graph for clade acI-C, showing disconnected components (gray) and the single largest component (green and black). Disconnected components are dropped prior to computing the network's seed sets because these groups of nodes are not connected to the bulk of the network. Within the single largest component, the giant strong component contains a substantial fraction of the compounds (green nodes), giving rise to a bow-tie structure in the metabolic network graph.

&nbsp;

## Supplementary Figure 5 (figures/fig2-sampling/Figure2.pdf)

&nbsp;

Mean estimated completeness of tribe-level (clade-level) population genomes as a function of the number of sampled genomes. For each tribe (clade), genomes were randomly sampled (with replacement) from the set of all genomes belonging to that tribe (clade). Completeness was estimated using 204 single-copy marker genes from the phylum Actinobacteria [@Parks2015]. Error bars represent the 95% confidence interval estimated from 1000 iterations.

&nbsp;

## Supplementary Figure 6 (figures/fig1-tree/FigureS5.pdf)

&nbsp;

Phylogenetic placement of the genomes used in this study within the acI lineage, relative to other sequenced actinobacterial genomes in the class Actinobacteria [@Gao2012] (Table S18). The tree was built using RAxML [@Stamatakis2014] from a concatenated alignment of protein sequences from 37 single-copy marker genes [@Darling2014]. The class Acidimicrobiia forms the outgroup. Vertical black bars indicate groups of genomes belonging to defined tribes/clades within the acI lineage, as determined using 16S rRNA gene sequences (for SAGs and bin FNEF8-2 bin_7 acI-B only) and a defined taxonomy [@Newton2011a]. SAGs are indicated with italic text.

&nbsp;

## References

&nbsp;
