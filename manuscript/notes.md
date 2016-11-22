# Stuff We Could Talk About But Probably Shouldn't
* pairwise ANI w/in each clade/tribe - will show presence of multiple species
* core and accessory genome of each clade - could say something based on tribes, but maybe opening a can of worms
* core and accessory genome of each tribe - unable to say anything b/c of completeness issues
* accumulation curves for acIs - haven't yet fully sampled the pangenome

# Removed from Intro in Primary MS

Microbes are critical players in freshwater systems, where they support essential ecosystem functions such as nutrient cycling. Of the freshwater bacteria, Actinobacteria of the acI lineage are among the most abundant, constituting upwards of 50% of the total bacteria in a variety of aquatic systems [@Zwart1998; @Glockner2000; @Zwart2002]. Despite their abundance, no isolates of the acI lineage have been stably propagated in pure culture.

Nevertheless, the acI lineage has been extensively studied in a community context using both DNA sequencing and single-cell targeted experiments. The acI have been phylogenetically divided into three clades (A, B, and C) and thirteen tribes on the basis of their 16S rRNA gene sequences [@Newton2011a]. Several studies have used fluorescent _in situ_ hybridization (FISH) and catalyzed reporter deposition (CARD) or microautoradiography (MAR) to identify substrate uptake capabilities of the acI. These studies reveal that the acI are capable of consuming amino acids [@Salcher2010; @Salcher2013], including the individual amino acids arginine, glutamate, glutamine, and leucine [@Buck2009; @Perez2010; @Salcher2010; @Eckert2012; @Salcher2013]); the saccharides glucose [@Buck2009; @Salcher2013], and N-acetylglucosamine (NAG) [@Beier2011; @Eckert2012; @Eckert2013]; the deoxynucleoside thymidine [@Perez2010; @Salcher2013], and acetate [@Buck2009]. However, due to limited phylogenetic resolution of some FISH probes, the studies cannot always link the uptake of these substrates to clades or tribes within the lineage.

Instead, metabolic reconstructions of single-cell genomes (SAGs) and metagenome-assembled genomes (MAGs) have been used to propose substrate uptake capabilities of clades aI-A and acI-B. To date, no genomes from clade acI-C have been included in these studies. These studies indicate both clades acI-A and acI-B are capable of consuming a wide array of N-containing compounds, including ammonium, branched-chain amino acids, polyamines, di-peptides, and cyanophycin [@Ghylin2014; @Garcia2015], with clade acI-A also capable of consuming oligopeptides [@Ghylin2014]. The lineage is also capable of consuming numerous saccharides, including the five-carbon sugars xylose, ribose, arabinose [@Garcia2013; @Ghylin2014; @Garcia2015] as well as poly- and oligo-saccharides [@Ghylin2014; @Garcia2015]. Notably, transporters for glucose and NAG have not yet been identified [@Garcia2013; @Ghylin2014], despite FISH studies showing uptake of those substrates [@Buck2009; @Beier2011; @Eckert2012; @Eckert2013; @Salcher2013]. Clade acI-B is also predicted to consume sucrose and maltose [@Garcia2015]; it also contains a chitinase for the breakdown of NAG [@Garcia2013; @Garcia2015]. Finally, the acI are predicted to contain the actinobacterial opsin protein actinorhodopsin [@Garcia2013; @Garcia2014; @Ghylin2014; @Garcia2015], a light-harvesting transmembrane protein [@Sharma2008; @Sharma2009], as well as the complete pathway for the biosynthesis of its cofactor retinal [@Ghylin2014]. Finally, a recent study of a metagenome-assembled genome from clade acI-B predicted the clade is unable to synthesize a number of essential compounds, including for the amino acids isoleucine, leucine, valine, tyrosine, tryptophan, phenylalanine, asparagine; and the cofactors biotin (Vitamin B7), cobalamin (Vitamin B12), folate (Vitamin B9), niacin (Vitamin B3), pantothenate (Vitamin B5), and riboflavin (Vitamin B2) [@Garcia2015]. In the aggregate, these results indicate the acI are photoheterotrophs, making a living on a diverse array of N-rich compounds, sugars, and oligo- and poly-saccharides. The acI do not appear to be metabolically self-sufficient, relying on other organisms for the production of essential nutrients.

These metabolic reconstructions all attempt to infer an organism's ecology from its genome content, a concept referred to as "reverse ecology" [@Levy2012]. While metabolic reconstructions represent a common entry point to reverse ecological analyses, other approaches take cues from systems biology, focusing instead on an organism's metabolic network. Here, the chemical reactions of metabolism are represented as connections between substrates and products, and analyzed from the perspective of the entire network [@Levy2012]. One such approach is the "seed set framework", which computes an organism's "seed set," the set of compounds that the organism cannot synthesize on its own and must exogenously acquire from its environment [@Borenstein2008]. As such, these compounds may represent both auxotrophies, essential metabolites for which biosynthetic routes are missing, and nutrients, for which routes for degradation (not synthesis) are present in the genome.

In this work, we expand existing genome-based analyses of the acI lineage by applying the seed set framework to a reference genome collection of 36 freshwater acI genomes, containing all previously-described acI genomes, as well as six additional SAGs and 15 MAGs, including for the first time genomes from the acI-C clade. We have developed a Python package to predict seed compounds for each clade, using metabolic network reconstructions generated from KBase (http://kbase.us). We also present the first metatranscriptomic analysis of gene expression across the three acI clades. Our use of the seed set framework enables the analysis of dozens of genomes without the need for a lengthy manual reconstruction step, and facilitates a focused analysis by identifying those compounds which an organism must obtain from its environment. The seed compounds predicted by our analysis are in agreement with experimental observations, confirming the ability of our method to predict an organism's metabolic requirements. Finally, our metatranscriptomic analysis shows that the acI express a diverse array of transporters, which we hypothesize may contribute to their observed dominance in a wide variety of aquatic systems.

# Removed from Methods in Primary MS

## Single-Cell Genome Generation, Selection, Sequencing, and Assembly

SAGs were collected, sequenced, and assembled as described previously [@Martinez-Garcia2012; @Garcia2013; @Ghylin2014]. SAGs were phylogenetically classified using partial 16S rRNA genes [@Martinez-Garcia2012] and a controlled nomenclature for freshwater bacteria [@Newton2011a] by insertion into references trees created in the ARB software package [@Ludwig2004]. Genome sequences are available through IMG (https://img.jgi.doe.gov/cgi-bin/mer/main.cgi), and can be accessed by searching for the IMG Taxon OIDs given in Table 1. Additional information is available in the Supplemental Online Material.

## Metagenome Sampling, Sequencing, Assembly, and Binning

Sample collection, processing, DNA sequencing, metagenomic assembly, genomic binning, and phylogenetic classification for the Trout Bog samples have been described previously [@Bendall2016]. With the exception of sample collection, identical procedures were followed for the Lake Mendota samples, for which depth-integrated water samples were collected from the top 12 meters at 96 time points during ice-free periods from 2008 to 2011. Metagenomic sequence reads are publicly available on the JGI Genome Portal (http://genome.jgi.doe.gov/) under Proposal ID 394. Metagenome sequences are available through IMG (https://img.jgi.doe.gov/cgi-bin/mer/main.cgi), and can be accessed by searching for the IMG Taxon OIDs given in Table 1. Additional information is available in the Supplemental Online Material.

## Removed b/c no longer used

## Metabolic Network Reconstruction and Reverse Ecology
### Re-annotation of Peptidases and Glycoside Hydrolases

Many seed compounds were associated with reactions catalyzed by peptidases or glycoside hydrolases, and genes associated with these reactions were re-annotated. Peptidase sequences were annotated using the MEROPS batch BLAST interface using default parameters [@Rawlings2015]. Glycoside hydrolases were first annotated using dbCAN [@Yin2012] to assign these genes to glycoside hydrolase families, as defined in the Carbohydrate-Active enZYmes Database CAZY [@Lombard2014]. Hidden Markov Models for these sub-families were then downloaded from dbCAN, and HMMER3 [@Eddy2011] was used to assign these genes to individual sub-families using default parameters.

# Removed from Results in Primary MS

## Non-canonical Biosynthetic Routes

While the reverse ecology pipeline presented here identified aspects of acI metabolism not observed previously, our approach does have some significant limitations. First, because our analysis relied on constructing clade-level composite genomes, we cannot say anything about the metabolism of individual tribes or interactions between them. And second, metabolic network analysis does not account for biological pathway organization, and predicts organisms may synthesize compounds via non-canonical routes (see the Supplemental Online Material for an example). As a result, the pipeline may under-predict the number of auxotrophies for a genome.

I will illustrate this using the biosynthesis of iso-leucine, which has previously been predicted to be an auxotroph for clade acI-B [@Garcia2015]. Iso-leucine is most commonly synthesized from threonine, via the following pathway.

__Need to check in IMG__

The enzyme for the first step is missing.
Step 	EC # 	Reaction
1 	4.3.1.19 	L-threonine --> (2Z)-2-aminobut-2-enoate --> 2-iminobutanoate --> 2-oxobutanoate
2 	2.2.1.6 	2-oxobutanoate --> (S)-2-aceto-2-hydroxybutanoate
3 	1.1.1.86 	(S)-2-aceto-2-hydroxybutanoate --> (R)-2,3-dihydroxy-3-methylpentanoate
4 	4.2.1.9 	(R)-2,3-dihydroxy-3-methylpentanoate --> (S)-3-methyl-2-oxopentanoate
5 	2.6.1.42 	(S)-3-methyl-2-oxopentanoate --> L-isoleucine

Thus, I would anticipate 2-oxobutanoate being predicted as a seed compound.

Iso-leucine can also be synthesized from pyruvate via the following pathway. Figure 2 from Trevor's paper suggests acI uses this pathway. The enzyme for the second step is missing.
Step 	EC # 	Reaction
1 	2.3.1.182 	pyruvate --> R-citramalate
2 	4.2.1.35 	R-citramalate --> (2R,3S)-3-Methylmalate
3 	1.1.1.85 	(2R,3S)-3-Methylmalate --> 2-Oxobutanoate
4 	2.2.1.6 	2-oxobutanoate --> (S)-2-aceto-2-hydroxybutanoate
5 	1.1.1.86 	(S)-2-aceto-2-hydroxybutanoate --> (R)-2,3-dihydroxy-3-methylpentanoate
6 	4.2.1.9 	(R)-2,3-dihydroxy-3-methylpentanoate --> (S)-3-methyl-2-oxopentanoate
7 	2.6.1.42 	(S)-3-methyl-2-oxopentanoate --> L-isoleucine

Thus, I would anticipate R-citramalate to be predicted as a seed compound.

Neither compound is predicted as a seed. So what is the route for iso-leucine biosynthesis?

Using the shortest_path algorithm included as part of the networkX package, I computed the shortest path for L-isoleucine synthesis from both pyruvate and threonine. The path from pyruvate follows:
Step 	Reaction
A 	Pyruvate --> 2-Hydroxyethyl-ThPP
B 	2-Hydroxyethyl-ThPP --> (S)-2-aceto-2-hydroxybutanoate
C 	(S)-2-aceto-2-hydroxybutanoate --> (R)-2,3-dihydroxy-3-methylpentanoate
D 	(R)-2,3-dihydroxy-3-methylpentanoate --> (S)-3-methyl-2-oxopentanoate
E 	(S)-3-methyl-2-oxopentanoate --> L-isoleucine

The path from synthesis from threonine includes two additional steps to convert threonine to pyruvate. In this route, steps C, D, and E are the same as steps 5, 6, 7 in the canonical pathway. With this route, the missing step is bypassed by directly converting pyruvate to (S)-2-aceto-2-hydroxybutanoate via 2-Hydroxyethyl-ThPP instead of via the intermediates citramalate, (2R,3S)-3-Methylmalate and 2-Oxobutanoate.

The presence of alternative pathways reveals why iso-leucine was not identified as a seed compound.

# Removed from Discussion in Primary MS

Our predictions of substrate utilization capabilities of the acI lineage are largely congruent with previous studies. We predict that the consumption of N-rich compounds is a universal feature of the acI lineage, with all three clades predicted to consume ammonium, branched-chain amino acids (leucine, isoleucine, and valine), the polyamines spermidine and putrescine, and oligopeptides. Further specialization may occur within each clade, as evidenced by each clade expressing unique substate binding proteins for some of their amino acid and peptide transporters. However, despite experimental observations of arginine, glutamate, and glutamine uptake, we failed to identify a transporter for these compounds.

Furthermore, we confirm the ability of all three clades to consume the five-carbon sugar xylose, as well as the six-carbon sugar maltose (previously thought to be restricted to clade acI-B). The acI-C genomes examined in this study did not contain transporters for ribose, suggesting that the utilization of this five-carbon sugar may be restricted to clades acI-A and acI-B. However, the possibility that acI-C consumes ribose cannot yet be ruled out, because our acI-C composite metabolic network graph remains incomplete. However, we failed to identify transporters for the saccharides glucose, and N-acetylglucosamine, both of which have been experimentally shown to be consumed by acI bacteria. Furthermore, in clades acI-A and acI-C, we identified additional hydrolases capable of acting on beta-glucosides, as well as alpha- and beta-galactosides, enzymes which had previously been detected only in clade acI-A.

We also identified transporters for the nucleobase uracil, as well as two permeases with broad specificity. These permeases are capable of acting on both purine and pyrimidine nucleobases (cysotine, uracil, and xanthine), suggesting the acI may obtain these compounds from the environment.

Our analysis also suggests that auxotrophies for some vitamins may be universal features of the lineage, as we predict all clades to be auxotrophic for Vitamins B5 and B6, in agreement with previous predictions for clade acI-B [@Garcia2015]. In addition, we predict transporters for Vitamins B1, B7 and B12, but do not predict auxotrophies for these vitamins. Finally, our analysis does not identify Vitamins B2, B3, B9, or B12 as auxotrophies for clade acI-B, a result which had been previously suggested [@Garcia2015]. This discrepancy may arise because we are analyzing the metabolism of the entire clade, while previous predictions were made on the basis of a single genome [@Garcia2015].

Transport proteins for many of these metabolites were among the most highly expressed in the genomes, suggesting that the success of acI's passive lifestyle may be due to constitutive expression of its diverse transport genes, which enable it to consume any substrate in its vicinity without waiting for a regulatory response. We also observe differences in the relative expression of these transporters, which may point to differences in the importance of these substrates to acI. For example, the transporters for oligopeptides and branched-chain amino acids are generally more highly expressed than those for sugars, suggesting a preference for compounds that can supply both nitrogen and carbon. The actinorhodopsin protein is highly expressed, and may facilitate synthesis of the ATP needed to drive acI's many ABC-type transporters.

Our analysis also provides new insights into auxotrophies within the acI lineage, identifying THF as an auxotrophy for clade acI-A, and lysine, homoserine, and UMP as auxotrophies for acI-C. THF is a derivative of folic acid (Vitamin B9), which was previously identified as an auxotrophy for clade acI-B [@Garcia2015]. Additionally, clade acI-B was previously identified as auxotrophic for a number of amino acids, though lysine and homoserine were not among them [@Garcia2015]. These results provide additional support to the hypothesis that distributed metabolic pathways and metabolic complementarity may be common features of freshwater bacterial communities [@Garcia2015].
