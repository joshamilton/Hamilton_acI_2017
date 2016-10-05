---
title: Figure Captions
author: Joshua J. Hamilton
date: September 16, 2016
---

&nbsp;

## Figure 1 (figures/fig1-tree/tree-abbrev.pdf)

&nbsp;

Phylogenetic placement of the SAGs and MAGs within the acI lineage, relative to other sequenced actinobacterial genomes in the class Actinobacteria [@Gao2012]. The tree was built using RAxML [@Stamatakis2014] from a concatenated alignment of protein sequences from 37 single-copy marker genes [@Darling2014]. The order Actinomycetales forms the outgroup. In addition, Supplemental Figure S1 shows the position of the acI lineage relative to other orders within the class Actinobacteria.

&nbsp;

## Figure 2a(b) (figures/fig2-sampling/completeness-2.svg and figures/fig2-sampling/completeness-1.svg)

&nbsp;

Mean estimated completeness of tribe(clade)-level population genomes. For each tribe(clade), genomes were randomly sampled (with replacement) from the set of all genomes belonging to that tribe(clade). Completeness was estimated using 204 single-copy marker genes from the phylum Actinobacteria [@Parks2015]. Error bars represent the 95% confidence estimated from 1000 iterations.

&nbsp;

## Figure 3 (figures/fig3-workflow/Figure3.pdf)

&nbsp;

Overview of reverse ecology pipeline, using three genomes from the acI-C clade as an example. __(A)__ Annotate microbial contigs using [KBase](http://kbase.us/), and build a metabolic network reconstruction from the annotations. For each genome, convert the metabolic network reconstruction to a metabolic network graph using custom Python scripts. In these graphs, metabolites are represented as nodes (circles) and reactions by arcs. __(B)__ Create a composite network graph for each clade by joining graphs for all genomes from that clade, and compute seed compounds for the composite graph. In individual genome graphs, grey nodes and edges indicate components of the composite graph missing from that genome graph. In the composite graph, seed compounds are shown in red. __(Inset)__ Three seed compounds which indicate an auxotrophy for L-homoserine, a methionine precursor. __(C)__ Map metatranscriptomic reads to each individual genome using [BBMap](https://sourceforge.net/projects/bbmap/). For each clade, identify orthologous gene clusters using OrthoMCL [@Li2003]. For each cluster, count all unique reads which map to any gene within that cluster using HTSeq [@Anders2014] and compute the relative gene expression using RPKM [@Mortazavi2008].

&nbsp;

## Figure 4 (figures/fig4-seeds/acI-aux-GH.pdf)

&nbsp;

Seed compounds of members of the acI lineage. Compounds in the lower panel are degraded by peptidases and glycoside hydrolases. For these compounds, the intensity of the color indicates the percentile average log2 RPKM of the encoding gene cluster. For compounds acted upon by multiple gene clusters, the percentile of the most highly-expressed cluster was chosen.

&nbsp;

## Figure 5 (figures/fig5-transporters/acI-transporters.pdf)

&nbsp;

Compounds actively transported by members of the acI lineage. The intensity of the color indicates the percentile average log2 RPKM of the encoding gene cluster. For multi-subunit transporters, the percentile of the most highly-expressed subunit was chosen.

&nbsp;

## Figure 6 (figures/fig6-actR/actR-expression.pdf)

&nbsp;

acI actinobacteria contain a retinal-based photosystem, comprising the opsin protein actinorhodopsin and the chromophore retinal. __(A)__ Pathway for the biosynthesis of retinal from _trans,trans_-farnesyl diphosphate. __(B)__ Relative expression of the genes for retinal and actinorhodopsin biosynthesis. The intensity of the color indicates the percentile average log2 RPKM of the encoding gene cluster. Note: The gene _blh_ was not found in the acI-C metagenome-assembled genomes examined in this study.

&nbsp;

## Supplementary Figure 1 (figures/fig1-tree/tree-full.pdf)

&nbsp;

Phylogenetic placement of the acI lineage, relative to other sequenced actinobacterial genomes in the class Actinobacteria [@Gao2012]. The tree was built using RAxML [@Stamatakis2014] from a concatenated alignment of single-copy marker genes [@Darling2014]. The class Acidimicrobiia forms the outgroup.

&nbsp;

## Supplementary Figure 2 (figures/fig3-workflow/FigureS2.pdf)

&nbsp;

Converting an unannotated genome to a metabolic network graph, for a simplified genome containing only glycolysis. __(A)__ Microbial contigs are annotated using KBase, and a metabolic network reconstruction is built from the annotations. The reconstruction provides links between protein-encoding genes in the genome and the enzymatic reactions catalyzed by those proteins. The reconstruction can be exported in a variety of formats, including a tabular format similar to the one shown here. __(B)__ The metabolic network reconstruction is converted to a hypergraph, in which metabolites are represented as nodes and reactions as hyperedges. In this representation, an edge can connect more than two nodes. For clarity, protons are not shown. __(C)__ The hypergraph is converted to a metabolic network graph, in which an edge can connect only two nodes. In this representation, a reaction is represented by a set of edges connecting all substrates to all products. The dotted line surrounds the currency metabolites. __(D)__ The metabolic network graph is then pruned, a process which removes all currency metabolites and any edges in which those metabolites participate. Representation of glycolysis after pruning. The images in (B) and (C) are modified from [@Ma2003].

&nbsp;

## Supplementary Figure 3 (figures/fig3-workflow/FigureS3.pdf)

&nbsp;

Identifying seed compounds in metabolic networks, using the same metabolic network as in Supplemental Figure S3. __(A)__ To identify seed compounds, the metabolic network graph is first decomposed into its strongly connected components (SCCs), sets of nodes such that each node in the set is reachable from every other node. Here, each set of circles nodes corresponds to a unique SCC. __(B)__ SCC decomposition enables seed sets to be identified from source components (components with no incoming edges) on the condensation of the original graph. In the condensation of the original graph shown here, each node corresponds to a unique SCC. This network has a single seed set, SCC_1, enclosed in a dotted circle. __(C)__ Seed compounds can be found from the mapping between SCCs and their constituent metabolites. In this example, glucose is the sole seed compound. While this particular result is probably intuitive, real metabolic networks are considerably more complex.

&nbsp;

## Supplementary Figure 4 (figures/fig3-workflow/FigureS4.pdf)

&nbsp;

Complete composite metabolic network graph for clade acI-C, showing disconnected components and the giant strongly connected components. Gray nodes and edges represent disconnected components which are dropped prior to computing the network's seed sets. Red nodes represent those present in the giant strongly connected component which contains the majority of the metabolites in the network.

&nbsp;

## References

&nbsp;
