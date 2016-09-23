---
title: Figure Captions
author: Joshua J. Hamilton
date: September 16, 2016
---

&nbsp;

## Figure 1 (figures/fig1-tree/tree-abbrev.pdf)

&nbsp;

Phylogenetic placement of the SAGs and MAGs within the acI lineage, relative to other sequenced actinobacterial genomes in the class Actinobacteria [@Gao2012]. The tree was built using RAxML [@Stamatakis2014] from a concatenated alignment of single-copy marker genes [@Darling2014]. The class Acidimicrobiia forms the outgroup.

&nbsp;

## Figure 2a(b) (figures/fig2-sampling/completeness-2.svg and figures/fig2-sampling/completeness-1.svg)

&nbsp;

Mean estimated completeness of tribe(clade)-level population genomes. For each tribe(clade), genomes were randomly sampled (with replacement) from the set of all genomes belonging to that tribe(clade). Completeness was estimated using single-copy marker genes [@Parks2015]. Error bars represent the 95% confidence estimated from 1000 iterations.

&nbsp;

## Figure 3 (figures/fig3-workflow/Figure3.pdf)

&nbsp;

Overview of reverse ecology pipeline. __(A)__ Annotate microbial contigs using [KBase](http://kbase.us/), and build a metabolic network reconstruction from the annotations. For each genome, convert the metabolic network reconstruction to a metabolic network graph using custom Python scripts. In these graphs, metabolites are represented as nodes (circles) and reactions by arcs. __(B)__ Create a composite network graph for each clade by joining graphs for all genomes from that clade, and compute seed compounds for the composite graph. In individual genome graphs, grey nodes and edges indicate components of the composite graph missing from that genome graph. In the composite graph, seed compounds are shown in red. __(Inset)__ Three seed compounds which indicate an auxotrophy for L-homoserine, a methionine precursor. __(C)__ Map metatranscriptomic reads to each individual genome using [BBMap](https://sourceforge.net/projects/bbmap/). For each clade, identify orthologous gene clusters using OrthoMCL [@Li2003]. For each cluster, count all unique reads which map to any gene within that cluster using HTSeq [@Anders2014] and compute the relative gene expression using RPKM [@Mortazavi2008].

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

Phylogenetic placement of the SAGs and MAGs within the acI lineage, relative to other sequenced actinobacterial genomes in the class Actinobacteria [@Gao2012]. The tree was built using RAxML [@Stamatakis2014] from a concatenated alignment of single-copy marker genes [@Darling2014]. The class Acidimicrobiia forms the outgroup.

&nbsp;

## Supplementary Figure 2 (figures/fig3-workflow/FigureS2.pdf)

&nbsp;

Pruning a metabolic network graph, using glycolysis as an example. __(A)__ Metabolic reconstruction for glycolysis. This panel also shows the KBase reaction IDs corresponding to each step of glycolysis, as well as the gene(s) encoding the enzyme for that step. Genes are shown only for the TE02754 genome. __(B)__ Glycolysis represented as a metabolic network graph. Currency metabolites are circumscribed by a dotted rectangle. For clarity, protons are not shown. __(C)__ Representation of glycolysis after pruning. The images in (A) and (B) are modified from [@Ma2003].

&nbsp;

## Supplementary Figure 3 (figures/fig3-workflow/FigureS3.pdf)

&nbsp;

Identifying seed compounds in metabolic networks. __(A)__ Metabolic network graph for an imaginary network. Grey circles represent metabolites and arcs represent reactions or parts of reactions. __(B)__ To identify seed compounds, the metabolic network graph is first decomposed into its strongly connected components (SCCs), sets of nodes such that each node in the set is reachable from every other node. Here, each colored set of nodes corresponds to a unique SCC. __(C)__ SCC decomposition enables seed sets to be identified from source components (components with no incoming edges) on the condensation of the original graph. In the condensation of the original graph shown here, each colored node corresponds to a unique SCC (and a potential seed se). Seed sets are enclosed in dotted circles. The pink seed set contains three seed compounds, the orange contains two, and the blue and yellow each contain one seed compound. The images in (A) and (B) are modified from [@Borenstein2008].

&nbsp;

## References

&nbsp;
