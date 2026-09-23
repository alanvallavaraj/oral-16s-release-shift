# Paired full-length versus V3–V4 oral 16S reference classification under a HOMD release shift

**Author manuscript, 23 September 2026 — exploratory preprint draft for author review**  
**Alan Immanuel Benjamin Vallavaraj**  
Affiliation and correspondence: to be confirmed by the author.

## Abstract

**Background.** Full-length 16S has previously improved taxonomic resolution in measured oral samples, but versioned reference collections allow a distinct test: how do full-length and V3–V4 representations perform on *the same new-release sequences* against an older reference collection?

**Methods.** We used Human Oral Microbiome Database (HOMD) v15.23 as the historical training reference and v16.03 as the temporal evaluation release. We removed every query identical to an old full-length sequence. For an exploratory paired comparison, degenerate 341F/805R primer extraction with up to three substitutions selected V3–V4 fragments from both releases. We also removed queries whose fragments duplicated a historical fragment. Each surviving v16.03 record was classified by binary 7-mer cosine nearest reference using either its full-length nucleotide sequence or its V3–V4 fragment. Both representations used the same 2,742 query IDs, the same HMT labels and corresponding older-release reference collections. The primary outcome was exact human microbial taxon (HMT) top-hit accuracy among 2,236 HMTs represented historically. The secondary outcome was AUROC separating known from 506 HMT-new query sequences. We estimated paired accuracy and AUROC differences with a 1,500-replicate stratified HMT-cluster bootstrap.

**Results.** Full-length classification correctly assigned **90.03%** of known-HMT queries, versus **80.23%** for V3–V4: an absolute paired gain of **9.79 percentage points** (95% HMT-cluster-bootstrap percentile interval **5.84 to 14.67 points**). Of 2,236 known queries, 1,763 were correct under both representations, 250 under full length alone, 31 under V3–V4 alone and 192 under neither. Paired novelty-detection AUROC was **0.8867** for full length and **0.8477** for V3–V4; the difference interval **−0.0128 to 0.0930** included zero. On the larger unfiltered temporal full-length benchmark, 7-mer cosine HMT assignment accuracy was **89.84%** and novelty AUROC **0.8353**.

**Conclusions.** Full-length reference information improves HMT assignment in this paired, reference-derived temporal test. The novelty-discrimination difference is uncertain. The paired analysis was devised after initial release inspection; independent isolate or specimen data are required to confirm clinical or prospective utility.

**Keywords:** oral microbiome; 16S rRNA; full-length sequencing; V3–V4; reference release; paired comparison; taxonomic assignment.

## 1. Introduction

HOMD establishes stable human oral taxon identifiers and curates reference sequences for named and provisionally named oral organisms [1,2]. Habitat-specific training sets and other curated oral references show why both classifier and database choice matter [3,4,6,16,17,30,31,40]. A database HMT identifier does not always denote a formally valid species, and taxonomy may change between reference versions [1,24]. Temporal reference evaluation therefore answers a different question from splitting randomly within a single release.

Oral microbiome studies increasingly compare full-length and short-region 16S measurements. Platform and primer choices, in addition to sequence length, affect taxonomic resolution and compositional profiles [5,9,11,23,36,38,39]. The 2024 PacBio comparison used saliva and oral biofilm samples, while a 2024 oral primer study varied trim length, primer and reference database [38,39]. **The general benefit of longer 16S reads is prior art**. Our contribution is a reproducible, same-query, versioned-reference comparison of HMT-level assignment and known-versus-new-HMT detection on curated oral references. The experiment uses in-silico sequences and should not be equated with a direct platform trial.

Reference-informed classifiers include RDP naive Bayes, IDTAXA, specialized k-mer methods, SpeciateIT and deep architectures [7,8,12–15,29]. SSUnique addresses novelty from a different phylogenetic perspective [32]. Benchmarks of profilers and region-specific assignment illustrate why labels, release timing and query construction must be explicit [10,18]. We compare one transparent 7-mer nearest-reference implementation under paired sequence representations and report auxiliary classifiers only as descriptive context; we do not claim to beat the published systems.

## 2. Materials and methods

### 2.1 Reference sources and temporal set

We obtained full-length 16S FASTA and matching taxonomy for HOMD v15.23 and v16.03 from the public HOMD reference archive. The historical release has 1,015 sequences and 781 HMT identifiers; the later release has 6,600 sequences. We marked a later-release record *known* if its HMT identifier was represented in v15.23 and *new* otherwise. We excluded every later sequence whose exact full-length string appeared anywhere in v15.23, even if its HMT label differed. The resulting complete temporal set has 5,158 known sequences from 530 known HMTs and 677 new sequences from 74 newly represented HMTs. Input URLs and four SHA-256 hashes are in the accompanying code and source_hashes.json. Exact-overlap exclusion does not remove near matches or guarantee independent isolate origins.

### 2.2 Paired full-length and V3–V4 comparison

We extracted V3–V4 from historical and newer full-length records using the degenerate primer strings 341F CCTACGGGNGGCWGCAG and 805R GACTACHVGGGTATCTAATCC. Edlib semi-global primer matching allowed up to three substitutions per primer; position and product-length bounds in the script restricted plausible regions [35,36]. These are reference-derived fragments, not PCR products, error-processed amplicon reads or an estimate of amplification bias [19,34].

We excluded each newer fragment exactly duplicated by *any* historical fragment. We then linked each surviving fragment by record ID to its full-length original. The paired set contains precisely the same 2,742 v16.03 query records under both representations: 2,236 known queries spanning 333 known HMTs and 506 new queries spanning 67 new HMTs. All paired full-length query strings were also checked against old full-length sequences for exact overlap. The paired full-length model used full historical reference sequences, and the paired fragment model used their extracted fragments. Both used the same HMT ground-truth labels, so the difference is not caused by comparing different query IDs. The primer/exact-fragment eligibility filter still changes which HMTs enter the experiment; the paired results generalize only to this selected subset.

### 2.3 Classification and statistical outcomes

We represented each reference and query as a binary character 7-mer vector and assigned the HMT of the reference with maximum cosine similarity. Reference order resolved ties. The main paired outcome was correct top-hit HMT assignment among known queries; newly added HMTs lack a matching historically available label and were excluded from this accuracy denominator. We separately treated maximum cosine as a novelty confidence score and computed known-versus-new HMT AUROC. AUROC summarizes ranking across thresholds rather than performance at one deployable cutoff [27,28]. The identical full-length and V3–V4 algorithm isolates the effect of the represented nucleotide span under this specific baseline.

For uncertainty we resampled *known HMT groups* with replacement and carried all their paired query IDs into the accuracy difference. We separately resampled known and new HMT groups for the AUROC difference. Both methods were evaluated on the same sampled indices in each of 1,500 draws (seed 482), giving percentile intervals for differences conditional on these data and references. Sequence counts vary across HMTs; cluster resampling accounts for some within-HMT dependence but does not sample patients or historical releases. We report the 2×2 correctness table instead of treating thousands of sequences as independent clinical replicates.

### 2.4 Context analyses

The complete temporal full-length set provides a non-paired reference outcome. It also includes a five-hit edit-similarity shortlist, a single untuned multinomial naive-Bayes baseline and cosine margins; those implementations are not exhaustive BLAST/VSEARCH or matched, optimised comparisons with RDP [20,29,35]. Five same-release HMT-withholding experiments quantify the optimism from testing on a current-release random split. A genome-derived 16SGOSeq reference collection [6] supplies a limited genus-label sanity check, but its sources partly overlap HOMD curation and do not provide independent species-level truth. DADA2, mothur, QIIME and QIIME 2 perform important specimen/read processing absent from our reference-only setup [19,22,33,37].

## 3. Results

### 3.1 Paired taxon assignment

| Paired result on 2,236 known-HMT queries | Full-length 16S | V3–V4 fragment | Full minus V3–V4 |
|:--|--:|--:|--:|
| Correct HMT top hit | 2,013 (90.03%) | 1,794 (80.23%) | +219 queries (+9.79 percentage points) |
| 95% HMT-cluster-bootstrap interval for accuracy difference | | | +5.84 to +14.67 percentage points |

Both representations assigned 1,763 queries correctly. Full length alone correctly assigned 250, the V3–V4 fragment alone assigned 31, and both missed 192. Thus the positive paired difference holds for an identical set of query IDs and a common HMT ground truth. Figure 1 reports these accuracy rates and the taxon-cluster interval.

### 3.2 Novelty ranking and release shift

On the same 2,742 query IDs, full-length cosine novelty AUROC was **0.8867** and V3–V4 cosine AUROC was **0.8477**, a difference of **0.0390**. The paired HMT-cluster-bootstrap 95% interval (**−0.0128 to 0.0930**) spans zero; we do not infer an established novelty-detection improvement. These AUROCs should not be directly compared with the larger complete temporal set because fragment eligibility selects the query HMT mix.

On all 5,835 full-length temporal queries after full-string overlap exclusion, raw cosine achieved known-HMT top-hit accuracy **0.8984** and known-versus-new AUROC **0.8353**. Five within-v16.03 taxon-withholding splits gave AUROCs **0.9289–0.9604**, a materially more favorable setting than the historical-to-newer-release comparison (Figure 2). In a retrospective test-defined 90%-known-coverage summary, the complete-set raw cosine score accepted 41.65% of new-HMT sequences; that operating threshold uses test labels and must not be interpreted as deployable. Shortlisted edit confidence yielded temporal AUROC 0.8116, an untuned naive-Bayes posterior 0.5472, and cosine margin 0.7412. Their limitations and detailed inventories appear in homd_deep_results.json.

### 3.3 Sensitivity and data sourcing

One-sequence-per-HMT resampling of the complete temporal set (300 sequence draws, with the same historical reference) gave raw-cosine AUROC 5th/median/95th percentiles **0.8212/0.8360/0.8503**; this represents sequence-choice variation within this dataset, not an independent-population confidence interval. Genome-derived 16SGOSeq genus concordance was 99.70% on 3,947 filtered sequences from 1,744 genomes; because source curation partly overlaps HOMD, that result is a sanity check, not external validation. Alternative catalogues and cohort studies show broader sources of variation beyond our design [21,25,26,31].

## 4. Discussion

The paired HMT accuracy gain is the principal positive result: full-length representation recovered roughly one in ten known-HMT assignments missed by V3–V4 on an identical set of later-release references. The 250-versus-31 discordance and taxon-bootstrap interval support the direction of that conditional association. This aligns with prior measured-platform reports [5,9,38], so the manuscript's value rests on its temporal release design, explicit identity checks, pairing and open reproduction, not on a claim that full-length 16S resolution was previously unknown.

A longer sequence did not show a similarly decisive improvement in novelty AUROC: its paired interval crossed zero. An HMT can be hard to distinguish even with more nucleotides when taxonomic boundaries are revised, references are sparse, or closely related organisms share sequence. Primer yield and exact-fragment exclusions also select a subset. The single later release is a proxy for future discovery rather than an independent biological cohort. In particular, a correct HMT top hit against versioned references does not establish diagnostic species identification, disease association, or specimen-level sensitivity [1,24,32].

This study is **exploratory**. The paired region design followed inspection of a previous temporal benchmark in the same data, so these results cannot be presented as preregistered confirmation. Some queries may share originating isolates or near-identical sequence variants with references after exact-string removal. We did not simulate sequencing errors, PCR efficiency, relative abundance, extraction, contamination, patient dependence or computational costs. Random test splits within a release performed better, reinforcing the need to use held-out releases and independent isolates or mock specimens before making prospective claims. Full-length versus V3–V4 wet-lab comparisons with frozen models, uniform taxonomy and subject-level separation are the decisive next experiment [5,38–40]. Power and sample-size work should be specified for that external study rather than inferred from the 2,236 correlated reference sequences [28].

## 5. Reproducibility and declarations

The repository contains the two evaluation scripts, archived JSON results, source hashes, dependency pins and figures. Run homd_data.py to download HOMD; fetch_external_data.py downloads the published 16SGOSeq reference and verifies its SHA-256. Run homd_deep_benchmark.py for contextual results and paired_region_experiment.py for the primary matched analysis; make_figures.py regenerates the plots. Public FASTA files are excluded from version control. No human participants, specimens or patient-level records were analysed. Author affiliation, funding, competing interests and target-journal formatting remain for author verification. This is an unreviewed manuscript draft, not a publication or clinical claim.

## References


1. Chen, T.; Yu, W.-H.; Izard, J.; Baranova, O. V.; et al. (2010). The Human Oral Microbiome Database: a web accessible resource for investigating oral microbe taxonomic and genomic information. *Database*. https://doi.org/10.1093/database/baq013

2. Escapa, Isabel F.; Chen, Tsute; Huang, Yanmei; Gajare, Prasad; et al. (2018). New Insights into Human Nostril Microbiome from the Expanded Human Oral Microbiome Database (eHOMD): a Resource for the Microbiome of the Human Aerodigestive Tract. *mSystems*. https://doi.org/10.1128/msystems.00187-18

3. F. Escapa, Isabel; Huang, Yanmei; Chen, Tsute; Lin, Maoxuan; et al. (2020). Construction of habitat-specific training sets to achieve species-level assignment in 16S rRNA gene datasets. *Microbiome*. https://doi.org/10.1186/s40168-020-00841-w

4. Griffen, Ann L.; Beall, Clifford J.; Firestone, Noah D.; Gross, Erin L.; et al. (2011). CORE: A Phylogenetically-Curated 16S rDNA Database of the Core Oral Microbiome. *PLoS ONE*. https://doi.org/10.1371/journal.pone.0019051

5. Esberg, Anders; Fries, Niklas; Haworth, Simon; Johansson, Ingegerd (2024). Saliva microbiome profiling by full-gene 16S rRNA Oxford Nanopore Technology versus Illumina MiSeq sequencing. *npj Biofilms and Microbiomes*. https://doi.org/10.1038/s41522-024-00634-1

6. Vázquez-González, Lara; Regueira-Iglesias, Alba; Balsa-Castro, Carlos; Tomás, Inmaculada; et al. (2025). A curated bacterial and archaeal 16S rRNA Gene Oral Sequences dataset. *Scientific Data*. https://doi.org/10.1038/s41597-025-05050-4

7. Salah, Rana; AbdElaal, Khlood R; Ghonaim, Lobna; Awe, Olaitan I; et al. (2026). DeepTaxa: a hybrid CNN-BERT framework for 16S rRNA taxonomic classification. *Bioinformatics Advances*. https://doi.org/10.1093/bioadv/vbag166

8. Murali, Adithya; Bhargava, Aniruddha; Wright, Erik S. (2018). IDTAXA: a novel approach for accurate taxonomic classification of microbiome sequences. *Microbiome*. https://doi.org/10.1186/s40168-018-0521-5

9. Wagner, Josef; Coupland, Paul; Browne, Hilary P.; Lawley, Trevor D.; et al. (2016). Evaluation of PacBio sequencing for full-length bacterial 16S rRNA gene classification. *BMC Microbiology*. https://doi.org/10.1186/s12866-016-0891-4

10. Odom, Aubrey R.; Faits, Tyler; Castro-Nallar, Eduardo; Crandall, Keith A.; et al. (2023). Metagenomic profiling pipelines improve taxonomic classification for 16S amplicon sequencing data. *Scientific Reports*. https://doi.org/10.1038/s41598-023-40799-x

11. Martínez-Porchas, Marcel; Villalpando-Canchola, Enrique; Vargas-Albores, Francisco (2016). Significant loss of sensitivity and specificity in the taxonomic classification occurs when short 16S rRNA gene sequences are used. *Heliyon*. https://doi.org/10.1016/j.heliyon.2016.e00170

12. Chaudhary, Nikhil; Sharma, Ashok K.; Agarwal, Piyush; Gupta, Ankit; et al. (2015). 16S Classifier: A Tool for Fast and Accurate Taxonomic Classification of 16S rRNA Hypervariable Regions in Metagenomic Datasets. *PLOS ONE*. https://doi.org/10.1371/journal.pone.0116106

13. Vinje, Hilde; Liland, Kristian Hovde; Almøy, Trygve; Snipen, Lars (2015). Comparing K-mer based methods for improved classification of 16S sequences. *BMC Bioinformatics*. https://doi.org/10.1186/s12859-015-0647-4

14. Ziemski, Michal; Wisanwanichthan, Treepop; Bokulich, Nicholas A.; Kaehler, Benjamin D. (2021). Beating Naive Bayes at Taxonomic Classification of 16S rRNA Gene Sequences. *Frontiers in Microbiology*. https://doi.org/10.3389/fmicb.2021.644487

15. Holm, Johanna B.; Gajer, Pawel; Ravel, Jacques (2024). SpeciateIT and vSpeciateDB: novel, fast, and accurate per sequence 16S rRNA gene taxonomic classification of vaginal microbiota. *BMC Bioinformatics*. https://doi.org/10.1186/s12859-024-05930-3

16. Myer, Phillip R.; McDaneld, Tara G.; Kuehn, Larry A.; Dedonder, Keith D.; et al. (2020). Classification of 16S rRNA reads is improved using a niche-specific database constructed by near-full length sequencing. *PLOS ONE*. https://doi.org/10.1371/journal.pone.0235498

17. Hsieh, Yu-Peng; Hung, Yuan-Mao; Tsai, Mong-Hsun; Lai, Liang-Chuan; et al. (2022). 16S-ITGDB: An Integrated Database for Improving Species Classification of Prokaryotic 16S Ribosomal RNA Sequences. *Frontiers in Bioinformatics*. https://doi.org/10.3389/fbinf.2022.905489

18. Seppey, Mathieu; Benavides, Andres; Berkeley, Matthew R.; Manni, Mosè; et al. (2026). LEMMIv2: benchmarking framework for metagenomic and 16S amplicon profilers with a catalogue of evaluated tools. *Genome Biology*. https://doi.org/10.1186/s13059-026-04089-9

19. Callahan, Benjamin J; McMurdie, Paul J; Rosen, Michael J; Han, Andrew W; et al. (2016). DADA2: High-resolution sample inference from Illumina amplicon data. *Nature Methods*. https://doi.org/10.1038/nmeth.3869

20. Edgar, Robert C. (2010). Search and clustering orders of magnitude faster than BLAST. *Bioinformatics*. https://doi.org/10.1093/bioinformatics/btq461

21. Wood, Derrick E.; Lu, Jennifer; Langmead, Ben (2019). Improved metagenomic analysis with Kraken 2. *Genome Biology*. https://doi.org/10.1186/s13059-019-1891-0

22. Caporaso, J Gregory; Kuczynski, Justin; Stombaugh, Jesse; Bittinger, Kyle; et al. (2010). QIIME allows analysis of high-throughput community sequencing data. *Nature Methods*. https://doi.org/10.1038/nmeth.f.303

23. Mukherjee, Chiranjit; Beall, Clifford J.; Griffen, Ann L.; Leys, Eugene J. (2018). High-resolution ISR amplicon sequencing reveals personalized oral microbiome. *Microbiome*. https://doi.org/10.1186/s40168-018-0535-z

24. Hackmann, Timothy J. (2025). Setting new boundaries of 16S rRNA gene identity for prokaryotic taxonomy. *International Journal of Systematic and Evolutionary Microbiology*. https://doi.org/10.1099/ijsem.0.006747

25. Rashidi, Armin; Gem, Hakan; McLean, Jeffrey S.; Kerns, Kristopher; et al. (2024). Multi-cohort shotgun metagenomic analysis of oral and gut microbiota overlap in healthy adults. *Scientific Data*. https://doi.org/10.1038/s41597-024-02916-x

26. Cha, Jun Hyung; Kim, Nayeon; Ma, Junyeong; Lee, Sungho; et al. (2025). A high-quality genomic catalog of the human oral microbiome broadens its phylogeny and clinical insights. *Cell Host & Microbe*. https://doi.org/10.1016/j.chom.2025.10.001

27. DeLong, Elizabeth R.; DeLong, David M.; Clarke-Pearson, Daniel L. (1988). Comparing the Areas under Two or More Correlated Receiver Operating Characteristic Curves: A Nonparametric Approach. *Biometrics*. https://doi.org/10.2307/2531595

28. Mattiello, Federico; Verbist, Bie; Faust, Karoline; Raes, Jeroen; et al. (2016). A web application for sample size and power calculation in case-control microbiome studies. *Bioinformatics*. https://doi.org/10.1093/bioinformatics/btw099

29. Wang, Qiong; Garrity, George M.; Tiedje, James M.; Cole, James R. (2007). Naïve Bayesian Classifier for Rapid Assignment of rRNA Sequences into the New Bacterial Taxonomy. *Applied and Environmental Microbiology*. https://doi.org/10.1128/aem.00062-07

30. Quast, Christian; Pruesse, Elmar; Yilmaz, Pelin; Gerken, Jan; et al. (2012). The SILVA ribosomal RNA gene database project: improved data processing and web-based tools. *Nucleic Acids Research*. https://doi.org/10.1093/nar/gks1219

31. McDonald, Daniel; Jiang, Yueyu; Balaban, Metin; Cantrell, Kalen; et al. (2023). Greengenes2 unifies microbial data in a single reference tree. *Nature Biotechnology*. https://doi.org/10.1038/s41587-023-01845-1

32. Lynch, Michael D. J.; Neufeld, Josh D. (2016). SSUnique: Detecting Sequence Novelty in Microbiome Surveys. *mSystems*. https://doi.org/10.1128/msystems.00133-16

33. Schloss, Patrick D.; Westcott, Sarah L.; Ryabin, Thomas; Hall, Justine R.; et al. (2009). Introducing mothur: Open-Source, Platform-Independent, Community-Supported Software for Describing and Comparing Microbial Communities. *Applied and Environmental Microbiology*. https://doi.org/10.1128/aem.01541-09

34. Bokulich, Nicholas A; Subramanian, Sathish; Faith, Jeremiah J; Gevers, Dirk; et al. (2013). Quality-filtering vastly improves diversity estimates from Illumina amplicon sequencing. *Nature Methods*. https://doi.org/10.1038/nmeth.2276

35. Šošić, Martin; Šikić, Mile (2017). Edlib: a C/C ++ library for fast, exact sequence alignment using edit distance. *Bioinformatics*. https://doi.org/10.1093/bioinformatics/btw753

36. Klindworth, Anna; Pruesse, Elmar; Schweer, Timmy; Peplies, Jörg; et al. (2012). Evaluation of general 16S ribosomal RNA gene PCR primers for classical and next-generation sequencing-based diversity studies. *Nucleic Acids Research*. https://doi.org/10.1093/nar/gks808

37. Bolyen, Evan; Rideout, Jai Ram; Dillon, Matthew R.; Bokulich, Nicholas A.; et al. (2019). Reproducible, interactive, scalable and extensible microbiome data science using QIIME 2. *Nature Biotechnology*. https://doi.org/10.1038/s41587-019-0209-9

38. Nagai, Takahiko; Shiba, Takahiko; Komatsu, Keiji; Watanabe, Takayasu; et al. (2024). Optimal 16S rRNA gene amplicon sequencing analysis for oral microbiota to avoid the potential bias introduced by trimming length, primer, and database. *Microbiology Spectrum*. https://doi.org/10.1128/spectrum.03512-23

39. Buetas, Elena; Jordán-López, Marta; López-Roldán, Andrés; D’Auria, Giuseppe; et al. (2024). Full-length 16S rRNA gene sequencing by PacBio improves taxonomic resolution in human microbiome samples. *BMC Genomics*. https://doi.org/10.1186/s12864-024-10213-5

40. Sierra, Maria A.; Li, Qianhao; Pushalkar, Smruti; Paul, Bidisha; et al. (2020). The Influences of Bioinformatics Tools and Reference Databases in Analyzing the Human Oral Microbial Community. *Genes*. https://doi.org/10.3390/genes11080878
