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

Metabolic network analysis does not account for biological pathway organization, and metabolic network graphs may suggest that microbes synthesize compounds via non-canonical routes. As a result, the seed set framework may under-predict the number of auxotrophies for a genome.

One example is the biosynthesis of iso-leucine, which has previously been predicted to be an auxotroph for clade acI-B [@Garcia2015]. Iso-leucine is most commonly synthesized from threonine, via the following pathway.

| Step  | EC # | Reaction | Reaction ID |
|-------|------|----------|-------------|
| 1 | 4.3.1.19 | L-threonine --> (2Z)-2-aminobut-2-enoate --> 2-iminobutanoate --> 2-oxobutanoate | rxn00737 |
| 2 | 2.2.1.6 | 2-oxobutanoate --> (S)-2-aceto-2-hydroxybutanoate | rxn03194 |
| 3 | 1.1.1.86 | (S)-2-aceto-2-hydroxybutanoate --> (R)-3-Hydroxy-3-methyl-2-oxopentanoate | rxn03436 |
| 4 | 4.2.1.9 | R)-3-Hydroxy-3-methyl-2-oxopentanoate --> 2,3-Dihydroxy-3-methylvalerate  | rxn03435 |
| 5 | 4.2.1.9 | (R)-2,3-dihydroxy-3-methylpentanoate --> (S)-3-methyl-2-oxopentanoate | rxn03437 |
| 6 | 2.6.1.42 | (S)-3-methyl-2-oxopentanoate --> L-isoleucine | rxn01575 |

The reaction for the first enzyme is absent from the metabolic reconstruction, and the seed set framework should predict 2-oxobutanoate as a seed compound. __Correction, the reaction is present in one genome, so acI-B is not auxotrophic for isoleucine.__

Iso-leucine can also be synthesized from pyruvate via the following pathway.

| Step  | EC # | Reaction | Reaction ID |
|-------|------|----------|-------------|
| 1 | 2.3.1.182 | pyruvate --> S-citramalate | rxn05109 |
| 2 | 4.2.1.35 | D-citramalate --> Citraconate | rxn02749 |
| 3 | No EC assigned | Citraconate --> (2R,3S)-3-Methylmalate | rxn02751 |
| 4 | 1.1.1.85 | D-erythro-3-Methylmalate --> 2-Oxobutanoate | rxn00735 |
| 5 | 2.2.1.6  | 2-oxobutanoate --> (S)-2-aceto-2-hydroxybutanoate | rxn03194 |
| 6 | 1.1.1.86 | (S)-2-aceto-2-hydroxybutanoate --> (R)-3-Hydroxy-3-methyl-2-oxopentanoate | rxn03436 |
| 7 | 4.2.1.9 | R)-3-Hydroxy-3-methyl-2-oxopentanoate --> 2,3-Dihydroxy-3-methylvalerate | rxn03435 |
| 8 | 4.2.1.9 | (R)-2,3-dihydroxy-3-methylpentanoate --> (S)-3-methyl-2-oxopentanoate | rxn03437 |
| 9 | 2.6.1.42 | (S)-3-methyl-2-oxopentanoate --> L-isoleucine  | rxn01575 |

The reaction for enzymes one through four are absent from the metabolic reconstruction, and the seed set framework should again predict 2-oxobutanoate as a seed compound.

However, 2-oxobutanoate is not predicted as a seed compound. Instead, the network indicates 2-oxobutanoate can be synthesized from pyruvate via pyruvate:thiamin diphosphate acetaldehydetransferase.

| Step  | EC # | Reaction | Reaction ID |
|-------|------|----------|-------------|
| 1 | No EC assigned | pyruvate + TPP --> 2-Hydroxyethyl-ThPP + CO2 | rxn00011 |
| 2 | 2.2.1.6  | 2-oxobutanoate + 2-Hydroxyethyl-ThPP--> (S)-2-aceto-2-hydroxybutanoate | rxn03194 |
| 3 | 1.1.1.86 | (S)-2-aceto-2-hydroxybutanoate --> (R)-3-Hydroxy-3-methyl-2-oxopentanoate | rxn03436 |
| 4 | 4.2.1.9 | R)-3-Hydroxy-3-methyl-2-oxopentanoate --> 2,3-Dihydroxy-3-methylvalerate | rxn03435 |
| 5 | 4.2.1.9 | (R)-2,3-dihydroxy-3-methylpentanoate --> (S)-3-methyl-2-oxopentanoate | rxn03437 |
| 6 | 2.6.1.42 | (S)-3-methyl-2-oxopentanoate --> L-isoleucine  | rxn01575 |

This pathway is for isoleucine biosynthesis is present in both Model SEED and KEGG, so it is probably a valid pathway in addition to those in MetaCyc. Thus, I conclude that acI-A is NOT auxotrophic for isoleucine.

After reviewing other proposed auxotrophies, I am unable to come up with a good example, so I will not discuss the issue of alternative pathways in the MS.

# Removed from Discussion in Primary MS

Our predictions of substrate utilization capabilities of the acI lineage are largely congruent with previous studies. We predict that the consumption of N-rich compounds is a universal feature of the acI lineage, with all three clades predicted to consume ammonium, branched-chain amino acids (leucine, isoleucine, and valine), the polyamines spermidine and putrescine, and oligopeptides. Further specialization may occur within each clade, as evidenced by each clade expressing unique substate binding proteins for some of their amino acid and peptide transporters. However, despite experimental observations of arginine, glutamate, and glutamine uptake, we failed to identify a transporter for these compounds.

Furthermore, we confirm the ability of all three clades to consume the five-carbon sugar xylose, as well as the six-carbon sugar maltose (previously thought to be restricted to clade acI-B). The acI-C genomes examined in this study did not contain transporters for ribose, suggesting that the utilization of this five-carbon sugar may be restricted to clades acI-A and acI-B. However, the possibility that acI-C consumes ribose cannot yet be ruled out, because our acI-C composite metabolic network graph remains incomplete. However, we failed to identify transporters for the saccharides glucose, and N-acetylglucosamine, both of which have been experimentally shown to be consumed by acI bacteria. Furthermore, in clades acI-A and acI-C, we identified additional hydrolases capable of acting on beta-glucosides, as well as alpha- and beta-galactosides, enzymes which had previously been detected only in clade acI-A.

We also identified transporters for the nucleobase uracil, as well as two permeases with broad specificity. These permeases are capable of acting on both purine and pyrimidine nucleobases (cysotine, uracil, and xanthine), suggesting the acI may obtain these compounds from the environment.

Our analysis also suggests that auxotrophies for some vitamins may be universal features of the lineage, as we predict all clades to be auxotrophic for Vitamins B5 and B6, in agreement with previous predictions for clade acI-B [@Garcia2015]. In addition, we predict transporters for Vitamins B1, B7 and B12, but do not predict auxotrophies for these vitamins. Finally, our analysis does not identify Vitamins B2, B3, B9, or B12 as auxotrophies for clade acI-B, a result which had been previously suggested [@Garcia2015]. This discrepancy may arise because we are analyzing the metabolism of the entire clade, while previous predictions were made on the basis of a single genome [@Garcia2015].

Transport proteins for many of these metabolites were among the most highly expressed in the genomes, suggesting that the success of acI's passive lifestyle may be due to constitutive expression of its diverse transport genes, which enable it to consume any substrate in its vicinity without waiting for a regulatory response. We also observe differences in the relative expression of these transporters, which may point to differences in the importance of these substrates to acI. For example, the transporters for oligopeptides and branched-chain amino acids are generally more highly expressed than those for sugars, suggesting a preference for compounds that can supply both nitrogen and carbon. The actinorhodopsin protein is highly expressed, and may facilitate synthesis of the ATP needed to drive acI's many ABC-type transporters.

Our analysis also provides new insights into auxotrophies within the acI lineage, identifying THF as an auxotrophy for clade acI-A, and lysine, homoserine, and UMP as auxotrophies for acI-C. THF is a derivative of folic acid (Vitamin B9), which was previously identified as an auxotrophy for clade acI-B [@Garcia2015]. Additionally, clade acI-B was previously identified as auxotrophic for a number of amino acids, though lysine and homoserine were not among them [@Garcia2015]. These results provide additional support to the hypothesis that distributed metabolic pathways and metabolic complementarity may be common features of freshwater bacterial communities [@Garcia2015].
