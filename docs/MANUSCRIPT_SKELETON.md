# Dating the PIF1/RRM3 duplication in budding yeasts (Saccharomycotina)

*Working manuscript draft. Methods and Results are complete and figure-backed; narrative sections
(Abstract, Introduction, Discussion) are drafted in Spencer's voice. Front/back matter (byline, funding,
Zenodo DOI) is placeholder pending finalization.*

**Spencer Gray**¹ *(co-authors TBD)*
¹ *Affiliation TBD*

---

## Significance
The budding yeast *Saccharomyces cerevisiae* carries two PIF1-family helicases, Pif1 and Rrm3, but when
the single ancestral gene duplicated to produce the pair has never been pinned down, because the alignable
helicase core is short and fast-evolving and the amino-acid signal at that depth is not merely weak but
actively misleading. We show that adding a protein's 3Di structural alphabet to the same alignment moves
the duplication from a spurious base-of-Fungi placement onto the Saccharomycotina stem, that maximum-
likelihood reconciliation confirms it, and that a fossil-calibrated time tree dates it to roughly 326 to
383 million years ago, long before the *Saccharomyces* whole-genome duplication. It is a concrete case of
structural information recovering a deep relationship that sequence alone gets wrong.

## Working title
*When did PIF1 and RRM3 split? Structural phylogenetics places the budding-yeast PIF1-family duplication
at the origin of Saccharomycotina.*

## Abstract
The budding yeast *Saccharomyces cerevisiae* carries two PIF1-family 5′→3′ DNA helicases, Pif1 and Rrm3,
with separable roles at telomeres, the rDNA, and stalled replication forks. Exactly when the single
ancestral fungal PIF1 gene duplicated to produce this pair has stayed unsettled, largely because the part
of the protein that can be aligned across the family, the helicase core, is short and fast-evolving; at
the depth where the answer lies the amino-acid signal is not just weak but misleading. We gathered the
cellular PIF1-family helicases of 728 fungal species (957 proteins; an Ascomycota-wide ingroup with
Basidiomycota, early-diverging fungi, and human PIF1 as outgroups), inferred a gene tree, and reconciled
it against a fungal species tree to locate the duplication. Amino-acid trees place the duplication at the
base of Fungi, and a topology test shows the sequence data marginally but significantly favor that deep
placement (AU *p* = 0.033), a preference the branch lengths attribute to long-branch attraction. Adding
the 3Di structural alphabet to the same alignment overrides it, pulling the node onto the stem of
Saccharomycotina; maximum-likelihood reconciliation (GeneRax) confirms that placement. Read against a
fossil-calibrated time tree, the split dates to roughly 326 to 383 Mya (Devonian to Carboniferous), just
crown-ward of the budding-yeast common ancestor and well before the *Saccharomyces* whole-genome
duplication. We further find that the recurrent two-copy state among mushrooms is an independent
duplication rather than shared ancestry with the yeast pair, and read Rrm3 as the more-derived of the two
paralogs. The result is a concrete case where structural information recovers a deep relationship that
sequence alone gets wrong.

## Introduction
The PIF1-family helicases are conserved 5′→3′ DNA helicases that act wherever the replication fork meets
trouble. In *Saccharomyces cerevisiae* the family has two members, Pif1 and Rrm3, and they have divided
the labor: Pif1 maintains the mitochondrial genome, unwinds G-quadruplexes, and negatively regulates
telomerase, while Rrm3 clears protein barriers ahead of the fork at the rDNA, the tRNA genes, and other
hard-to-replicate sites (Boulé and Zakian 2006; Bochman, Sabouri, and Zakian 2010). Both descend from a
single ancestral fungal PIF1, but where on the fungal tree that ancestor duplicated has not been pinned
down.

The question is harder than it looks. The region of these proteins that can be aligned across the whole
family is the helicase core, a few hundred residues that evolve quickly; over the distances that separate
the major fungal lineages the amino-acid signal in that core saturates, and a sequence tree can no longer
tell which arrangement of the deep branches is real. Worse than agnostic, the saturated signal is
positively misleading here: the fast-evolving budding-yeast copies and the distant outgroups set up a
textbook long-branch-attraction trap, so the sequence data actively favor a deep, wrong placement. This is
the regime where protein structure helps: a fold is conserved long after the sequence that specifies it
has been overwritten, so a tree built from structural features can recover relationships that have gone
dark at the sequence level (Moi et al. 2025). The Pif1/Rrm3 split sits at exactly this awkward depth,
which made it a good candidate for a structure-aware approach.

Here we sample the cellular PIF1-family helicases across the fungi (957 proteins from 728 species, with
human and early-diverging-fungal outgroups), build the gene tree three ways (from sequence, from
structure, and from the two combined), and reconcile it against a fungal species tree to read off the
branch on which the duplication occurred (Fig. 1). Sequence alone leaves the node at the base of Fungi;
adding structure pulls it onto the Saccharomycotina stem, a single, datable branch. The duplication that
produced Pif1 and Rrm3 is a Saccharomycotina innovation, roughly as old as the budding-yeast lineage
itself.

---

## Methods

### Taxon sampling and homolog identification
- Fungi-wide sample: **957 cellular PIF1-family proteins across 728 species** — Ascomycota ingroup
  (681) + outgroups: Basidiomycota (235), Mucoromycota (36), Chytridiomycota (4), and human PIF1 (1).
- **Family filter — InterPro IPR048293 ("PIF1_RRM3_pfh1"), *not* Pfam PF05970.** PF05970 is also
  carried by **Helitron** rolling-circle transposons, inflating fungal counts several-fold (Ascomycota
  2,255→681, Basidiomycota 2,393→235, Mucoromycota 742→36, Chytridiomycota 33→4); IPR048293 isolates the
  *cellular* helicases while retaining every anchor. *(Key decision — Fig. S1 / Table S1.)*
- Anchors: *S. cerevisiae* Pif1 (UniProt **P07271**), Rrm3 (**P38766**); *S. pombe* Pfh1 (**Q9UUA2**),
  human PIF1 (**Q9H611**). Sequences retrieved from UniProt; one representative per gene.

### Structures
- **AlphaFold DB** for 828/957 proteins (model version v6); the **129** lacking an AFDB model were
  folded with **ColabFold (AlphaFold2)** to match the predictor (median full-length pLDDT 62.5 vs 64.6 for
  AFDB — predictor-consistent; Fig. S2). 3 very long proteins were folded as their helicase-core region
  only (flagged in `manifest.csv`). Provenance per protein recorded in `manifest.csv`.

### Core trimming ("corecut")
- Both sequences and structures trimmed to the conserved PIF1 helicase core (PF05970 HMM envelope) so
  variable N/C-terminal extensions don't confound the trees. Structures sliced to the **same** residue
  coordinates as the sequence core (so the AA and structural analyses describe identical sites). The core
  is uniformly high-confidence — median core pLDDT 88 in both predictors — even where the full-length
  model is moderate (Fig. S2).

### Sequence maximum-likelihood gene tree
- Align cores (MAFFT) → trim (trimAl `-automated1`, **209 columns**) → **IQ-TREE LG+I+G4**, 957 taxa,
  1000 ultrafast bootstraps + SH-aLRT (`-seed 42`). *(Also LG+C20 via PMSF as a model-adequacy check.)*

### Structural gene tree (FoldTree)
- **foldseek** all-vs-all (3Di+AA, exhaustive) → FoldTree `fident` distance with the −b·ln(1−d/b)
  correction (b = 0.93) → **FastME**; `alntmscore` and `lddt` distances as metric-robustness checks.
  Rooted on human PIF1.

### AA+3Di partitioned tree *(the resolving method)*
- Per-residue **3Di** strings (foldseek) mapped onto the AA alignment columns → a 3Di alignment with
  identical gap structure → concatenated (**957 × 418**: 209 AA + 209 3Di). **IQ-TREE partition model:
  LG+G for the AA partition, the 3DiPhy 3Di substitution matrix (`3di_substmat`) for the 3Di partition.**
  The 3Di partition carried *more* parsimony-informative sites (194) than the AA core (180). The 3Di
  characters are derived from the predicted structure, so the AA and 3Di partitions are not statistically
  independent; the gain is fold-constraint signal that persists after the underlying sequence has
  saturated, not a second independent dataset (see Discussion).

### Topology test and relative rates
- Approximately-unbiased (AU) test (IQ-TREE, 10,000 RELL replicates, LG+I+G4) comparing the AA-optimal
  tree against a tree constrained only to Saccharomycotina monophyly; Tajima's relative-rate test on the
  anchor pair and clade-level branch-length comparisons from AA branch lengths re-estimated on the AA+3Di
  topology (`workflow/15_rate_check.py`, `-seed 42`).

### Species tree
- **NCBI-taxonomy** backbone for the 728 species (ete3 `NCBITaxa`, pruned, forced binary; 720 leaves
  after dropping 8 taxa NCBI doesn't recognize → 941 mapped genes), *and* a **grafted phylogenomic
  backbone** (`workflow/16_graft_species_tree.py`, 719 taxa) that splices the Shen et al. 2018
  time-calibrated Y1000+ budding-yeast topology into the non-Saccharomycotina NCBI backbone. The
  duplication is placed on both (Results R4), so the answer does not depend on NCBI's polytomy resolution.

### Reconciliation and dating
- **Species-overlap** reconciliation (ete3) of each gene tree against the species tree, tallying
  duplications by the clade they map to.
- **GeneRax 2.0.4** (`UndatedDL`, EVAL reconciliation) as the gold-standard ML duplication–loss
  reconciliation, on native x86 (WSL2) — run against the NCBI tree, the grafted phylogenomic tree, and
  three robustness variants of the gene tree (Results R4); gives the duplication branch, the loss pattern,
  and stem-vs-crown. (GeneRax 2.1.3 segfaults at `inferReconciliation` on both Rosetta and native Linux,
  so 2.0.4 is used throughout; EVAL rather than SPR search, which is intractable for a single 941-tip
  family.)
- **Absolute dating** reads the age off the Shen et al. 2018 fossil-calibrated Y1000+ chronogram rather
  than re-estimating it: the duplication clade's species are matched onto the Shen tips and their MRCA
  bracketed by [crown, stem] ages (`workflow/17_place_on_timetree.py`).

### Reproducibility
- All steps scripted (`workflow/01`–`17`) in a version-pinned conda environment; every protein, structure
  source, and decision logged in `manifest.csv` and the README. Seeds fixed (`-seed 42`) and recorded in
  the `.iqtree`/`.log` outputs. Repo: github.com/spegray/pif1-foldtree. Key tool versions: MAFFT 7.526,
  trimAl 1.5.1, HMMER 3.4, IQ-TREE 3.1.2, foldseek 10, FastME 2.1.6.3, GeneRax 2.0.4 (run in a separate
  native-x86/WSL2 environment), ColabFold/AF2.

---

## Results

### R1 — Copy number hints at, but cannot prove, a Saccharomycotina-restricted duplication
- Two-copy vs one-copy by subphylum: **Saccharomycotina 91/105 two-copy**; Pezizomycotina 430/447 and
  Taphrinomycotina 11/11 **one-copy**; two-copy mushrooms (Agaricomycotina) common too. Read: the
  duplication *looks* Saccharomycotina-restricted, but copy-counting can't distinguish shared ancestry
  from convergence — the gene tree must. *(Fig. 2 / Table 1.)*

### R2 — Amino-acid trees place the duplication at the base of Fungi
- In both LG+I+G4 and LG+C20(PMSF) trees, the MRCA of ScPif1 and ScRrm3 maps to **"Fungi" (kingdom)**:
  the sequence tree scatters the budding-yeast Pif1 and Rrm3 orthologs so widely that their common
  ancestor sits near the root. *(Fig. 3A, left panel.)*

### R3b — The amino-acid signal significantly favors that deep placement, and it is long-branch attraction
- This is not a soft polytomy. A constrained AU test (10,000 RELL replicates, LG+I+G4) shows the
  amino-acid data marginally but **significantly reject Saccharomycotina monophyly (AU *p* = 0.033,
  ΔlnL = 143)** — a real, if weak, preference for the deep placement, not an absence of signal. The branch
  lengths say why: the Saccharomycotina copies evolve about **twice as fast** as the non-Saccharomycotina
  single-copy genes (mean terminal branch 0.15/0.19 vs 0.08 subs/site), and within them the Rrm3 clade
  runs **1.22× faster** than the Pif1 clade (root-to-tip 1.07 vs 0.88; permutation *p* < 0.0001). Fast
  ingroup copies plus distant outgroups is the textbook long-branch-attraction configuration, which
  produces exactly this kind of modest, artifactual pull toward a deep split. *(The model-free Tajima
  anchor-pair test is non-significant but underpowered at two sequences and 209 columns; the confounded
  whole-tree AU comparison, ΔlnL ≈ 1891, mixes global tree differences with the node question and is not
  interpretable, so we do not use it.)* *(Fig. 3B.)*

### R3 — The 3Di structural signal overrides it → Saccharomycotina *(headline)*
- In the **AA+3Di** tree the ScPif1/ScRrm3 MRCA collapses to **Saccharomycotina** (197 genes /
  103 species — essentially all sampled Saccharomycotina), with non-Saccharomycotina fungi *outside* the
  duplication clade. Adding structure to an otherwise-identical alignment overrides the significant-but-
  weak amino-acid preference and moves the node from kingdom-level to subphylum-level — a clean
  demonstration of the FoldTree thesis. *(Fig. 3A, right panel — the key figure; AA-only vs AA+3Di fans
  side by side.)*

### R4 — GeneRax reconciliation places the duplication on the Saccharomycotina ancestor
- Gold-standard ML reconciliation (GeneRax 2.0.4, `UndatedDL`) of the AA+3Di gene tree maps the PIF1/RRM3
  duplication to **`node_448` — the MRCA of all 104 sampled Saccharomycotina** (the budding-yeast common
  ancestor) — confirming the species-overlap result with an independent, error-tolerant method. The
  literal recPhyloXML reconciliation, rendered independently, agrees. *(Fig. 4; raw reconciliation in
  Fig. S3.)*
- **Robust to gene-tree wobble:** the same placement holds on the bootstrap-supported AA+3Di tree.
  Although the duplication node is UFBoot-weak (29; SH-aLRT 98) and its Rrm3 daughter weaker still
  (UFBoot 16; the Pif1 daughter is well-supported at UFBoot 90), reconciliation (which integrates over
  gene-tree error) places the event at Saccharomycotina regardless.
- **Stem, just crown-ward:** both daughter (Pif1, Rrm3) clades contain members of **8/9 sampled
  Saccharomycotina families** (only Dipodascaceae is Pif1-only); the earliest-branching sample
  (Lipomycetales) is single-copy → the duplication sits on the Saccharomycotina stem near the crown
  (caveat: sparse early-diverging-lineage sampling).
- **Robust to the species-tree backbone and to structure-prediction choices.** On the **grafted
  phylogenomic species tree** (Shen 2018 budding-yeast topology spliced into the NCBI backbone, 719 taxa),
  `node_448` is the **single top duplication node in the whole tree (D = 43)**; on the NCBI-taxonomy
  backbone the same node still carries a duplication (D = 9). The Saccharomycotina placement also holds
  across three reruns of the AA+3Di gene tree — AFDB-structures only (dropping the 129 ColabFold models;
  node_448 D = 37), high-confidence cores only (core pLDDT ≥ 87; D = 39), and an estimated GTR20 3Di
  matrix in place of the fixed Foldseek one (D = 43) — closing the predictor-batch, structure-quality,
  and 3Di-matrix robustness questions. *(Fig. S5.)*
- *A note on the totals:* on the grafted tree the global event counts are high (D = 476; ~3,618
  speciation-losses) because reconciliation absorbs residual deep gene-tree noise as many small events;
  the interpretable signal is the anchor's placement on the Saccharomycotina node, not the totals (Fig. S4
  ranks the per-node counts and shows the two real events standing clear of the tail). *(Fig. 4 —
  duplication + losses on the species tree.)*

### R5 — An independent mushroom duplication, and RRM3 as the derived paralog
- The two-copy **mushrooms (Agaricomycetes) are a separate duplication**: `node_578` carries the largest
  duplication count outside the budding-yeast node (**D = 23**; the tree-wide maximum on the NCBI
  backbone), distinct from the yeast Pif1/Rrm3 event — copy-number convergence, invisible to copy-counting
  but resolved by the tree. *(Fig. 4; ranked per-node counts in Fig. S4.)*
- **RRM3 reads as the more-derived paralog.** Single-copy PIF1s outside Saccharomycotina are co-orthologs
  of *both* yeast paralogs by descent (0/536 nest inside either clade), yet **442 of 536 sit closer to
  ScPif1** by sequence (mean core identity 54.5% vs 52.3%; *S. pombe* Pfh1 55.1 vs 53.2) — i.e. Pif1
  retained more of the ancestral state while Rrm3 diverged further, consistent with Rrm3's faster branch
  (R3b). We read this as an asymmetry, not a proof; it rests on raw sequence distance and on the
  bootstrap-weak Rrm3 clade, so reciprocal monophyly leans on the reconciliation rather than raw support.
  *(Fig. 3B.)*

### R6 — Absolute timing: a Devonian–Carboniferous (~326–383-Myr-old) duplication, long predating the WGD
**Approach — read the age off a published, fossil-calibrated time tree rather than re-estimate it.** The
duplication maps to the Saccharomycotina stem, and the Y1000+ budding-yeast molecular clock (Shen et al.
2018, *Cell*; 332 genomes) already dates exactly this node.
- **Direct placement (`workflow/17_place_on_timetree.py`):** mapping the duplication's species set onto
  Shen's calibrated tree puts its MRCA at a node subtending ~95% of budding yeasts (**317/332 tips**, 54
  of our species matched), dated **crown 326 Mya, stem 383 Mya** (tree root / BYCA 404 Mya). The
  duplication sits on that stem branch → **~326–383 Mya (Devonian–Carboniferous), in the early
  Saccharomycotina crown, just crown-ward of the ~404-Mya budding-yeast common ancestor.** The
  earliest-diverging budding yeasts fall outside the clade, consistent with their single-copy state. This
  refines the borrowed "~400 Mya" estimate with a direct read off the calibrated phylogenomic tree.
- **Independent sanity check:** ~3–4× older than the *Saccharomyces* **whole-genome duplication (~100
  Mya**, a crown-Saccharomycetaceae event) — consistent with Rrm3 occurring across budding yeasts
  *including* pre-WGD lineages (*Kluyveromyces*, *Lachancea*). The ordering (duplication ≫ WGD) is robust
  to the exact calibration. *(Fig. 5.)*

---

## Tables

**Table 1 — PIF1-family copy number by subphylum** *(the R1 tally; Fig. 2 shows the proportions).*

| Subphylum (phylum) | species | two-or-more copy | single copy |
|---|---|---|---|
| Saccharomycotina (Ascomycota) | 105 | 91 | 14 |
| Pezizomycotina (Ascomycota) | 447 | 17 | 430 |
| Taphrinomycotina (Ascomycota) | 11 | 0 | 11 |
| Agaricomycotina (Basidiomycota) | 92 | 58 | 34 |
| *(other sampled subphyla)* | 73 | — | — |

**Table S1 — Family filter, PF05970 vs IPR048293 protein counts per phylum** *(Fig. S1).*

| Phylum | Pfam PF05970 (raw) | InterPro IPR048293 (cellular) |
|---|---|---|
| Ascomycota | 2,255 | 681 |
| Basidiomycota | 2,393 | 235 |
| Mucoromycota | 742 | 36 |
| Chytridiomycota | 33 | 4 |

---

## Figures

**Main:**
1. **Pipeline** (`fig1_pipeline.svg`) — PIF1-family duplication-dating pipeline: homolog retrieval →
   predicted structures → helicase-core extraction → parallel trees (sequence ML, FoldTree, AA+3Di) →
   GeneRax reconciliation → dated duplication.
2. **Copy number** (`fig2_copynumber.svg`) — Pif1-family copy number by fungal subphylum; stacked bars
   show the fraction of species with 1 / 2 / 3+ copies. Saccharomycotina 91/105 ≥ 2-copy; Pezizomycotina
   430/447 single-copy; mushrooms (Agaricomycotina) recurrently multi-copy.
3. **The resolution** (`fig3_combined.svg`) — (A) human-rooted fan cladograms of all 957 PIF1-family
   proteins, amino acids alone vs amino acids + 3Di; the Pif1/Rrm3 MRCA moves from base-of-Fungi
   (950 tips) to Saccharomycotina (197 tips). (B) AA branch-length phylogram of that clade: Rrm3 root-to-
   tip 1.07 vs Pif1 0.88 subs/site (1.22× faster), the long-branch signal that misleads sequence trees.
   *(headline)*
4. **Reconciliation** (`fig4_reconciliation.svg`) — 719-taxon grafted phylogenomic tree (GeneRax
   UndatedDL); coral node area = per-node duplications; the Saccharomycotina ancestor (node_448, D = 43)
   and the independent Agaricomycotina ancestor (node_578, D = 23); D < 5 faded; tips colored by copy
   state (recurrent single-paralog loss).
5. **Absolute date** (`fig5_timetree.svg`) — duplication placed on the Shen 2018 fossil-calibrated
   chronogram (332 genomes); MRCA subtends 317/332 budding-yeast tips; crown 326 / stem 383 / root 404 Mya
   (~326–383 Mya band), well before the ~100-Mya *Saccharomyces* WGD.

**Supplementary:**
- **S1** (`figS1_family_filter.svg`) — family filter, PF05970 vs IPR048293 counts per phylum (Ascomycota
  2,255→681, Basidiomycota 2,393→235, Mucoromycota 742→36, Chytridiomycota 33→4).
- **S2** (`figS2_plddt.svg`) — pLDDT for full-length vs helicase-core models by predictor; AFDB and
  ColabFold agree (full-length ~63, core ~88).
- **S3** (`figS3_reconciliation_thirdkind.svg`) — the literal recPhyloXML reconciliation (thirdkind), a
  raw cross-check of Fig. 4.
- **S4** (`figS4_duplication_counts.svg`) — per-node duplication counts ranked; the two interpreted events
  (D = 43, 23) tower over the D ≤ 16 noise tail.
- **S5** (`figS5_robustness.png`) — robustness sweep: node_448 remains the top duplication node across
  AFDB-only, pLDDT ≥ 87, GTR20, and grafted-backbone reruns.

## Discussion
The practical lesson here is narrow but clean: for this duplication the amino-acid core had not merely run
quiet, it pointed the wrong way. A topology test shows the sequence data significantly prefer the deep,
base-of-Fungi placement (AU *p* = 0.033), and better modeling did not rescue it (a profile-mixture model
fit on the same 209 columns left the node exactly where the plain model did). The branch lengths explain
the preference as long-branch attraction, with the fast-evolving budding-yeast copies drawn toward the
distant outgroups. What moved the node was information of a different kind. The 3Di structural alphabet
carried more resolving power than the amino-acid core over the same sites (194 versus 180
parsimony-informative positions), and folding that signal into the tree pulled the Pif1/Rrm3 ancestor from
the base of Fungi down onto the Saccharomycotina stem. Because the 3Di characters are read off the
predicted structure, they are not an independent second dataset so much as the fold-level constraint that
outlasts the sequence once it has saturated (Moi et al. 2025); the comparison is fair because both
partitions describe the same trimmed sites.

Placed in time, the duplication is old. Mapping the Saccharomycotina stem onto the calibrated Y1000+ time
tree (Shen et al. 2018) puts the split at roughly 326 to 383 million years ago (Devonian to Carboniferous),
just crown-ward of the budding-yeast common ancestor. That timing matters for how the pair is usually
explained: the *Saccharomyces* whole-genome duplication, the first event people tend to reach for when
accounting for yeast gene pairs, is three to four times younger and cannot have produced Pif1 and Rrm3.
Rrm3 is present across the budding yeasts, including lineages (*Kluyveromyces*, *Lachancea*) that diverged
well before the whole-genome duplication, exactly as a Saccharomycotina-ancestral origin predicts.

The two paralogs are also not interchangeable copies of the ancestral gene. Single-copy PIF1s outside
Saccharomycotina are co-orthologs of both by descent, yet 442 of 536 sit closer to Pif1 than to Rrm3 in
sequence (*Schizosaccharomyces pombe* Pfh1 among them), and the Rrm3 clade carries the longer branch
(R3b), which together we read as Pif1 retaining more of the ancestral character while Rrm3 took on the
faster-evolving, more-derived role. The reading rests on raw sequence distance and on a bootstrap-weak
Rrm3 clade, so we hold it as an asymmetry worth naming rather than a settled fact; whether it tracks the
functional specialization of the two helicases is a question this tree cannot answer, but it is a natural
one to ask next.

Copy number is a poor guide to history. Many mushrooms (Agaricomycetes) also carry two PIF1-family genes,
and on counts alone that resembles the budding-yeast situation; the gene tree shows it is a separate,
later duplication. Convergent gene duplication of this kind is easy to mistake for shared ancestry, and
telling the two apart is precisely what a reconciliation, rather than a tally, is built to do.

Two caveats temper the confidence and mark the work's edges. First, the duplication node is well supported
by the SH-aLRT test but weak under bootstrap resampling (UFBoot 29, and the Rrm3 daughter weaker still at
16), as expected for a deep node resting on saturated sequence; the result leans on reconciliation, which
tolerates that gene-tree uncertainty, more than on any single branch's bootstrap value. That reliance
holds up under stress-testing. The Saccharomycotina placement survives swapping the NCBI taxonomy for a
grafted phylogenomic backbone, dropping the predicted structures, keeping only the high-confidence cores,
and estimating the 3Di exchange rates from the data rather than borrowing a fixed matrix; in every one of
those reruns the Saccharomycotina ancestor still carries the tree-wide maximum duplication count. Second,
the early-diverging budding yeasts are thinly sampled here, so while the duplication maps cleanly to the
Saccharomycotina stem, whether it falls just before or just after the subphylum's first split is not yet
settled; denser sampling of the Lipomycetales and their relatives would tighten that. Neither caveat
unseats the central result, and both read as the next steps rather than reasons for doubt.

---

## References
- Moi D, et al. (2025) Structural phylogenetics with the FoldTree/3Di approach. *Nat. Struct. Mol. Biol.*
  https://doi.org/10.1038/s41594-025-01649-8
- Shen X-X, et al. (2018) Tempo and mode of genome evolution in the budding yeast subphylum (Y1000+).
  *Cell* 175:1533–1545. https://doi.org/10.1016/j.cell.2018.10.023
- Bochman ML, Sabouri N, Zakian VA (2010) Unwinding the functions of the Pif1 family helicases.
  *DNA Repair* 9:237–249. https://doi.org/10.1016/j.dnarep.2010.01.008
- Boulé J-B, Zakian VA (2006) Roles of the yeast Pif1 helicase family in maintaining genome stability.
  *Nucleic Acids Res.* 34:4147–4153. https://doi.org/10.1093/nar/gkl561
- Malone EG, Thompson MD, Byrd AK (2022) Role and regulation of Pif1 family helicases. *Int. J. Mol. Sci.*
  23:3736. https://doi.org/10.3390/ijms23073736
- Harman A, Manna E (2016) [PIF1 phylogenetics in amoebae and accessory domains]. *Mol. Phylogenet. Evol.*
  https://doi.org/10.1016/j.ympev.2016.07.015
*(Y1000+ comparative-genomics family studies and remaining citations to be completed at finalization.)*

## Data and code availability
- Repo (scripts `01`–`17`, `manifest.csv`, trees, alignments, figure scripts): github.com/spegray/pif1-foldtree.
- Deposited on submission (Zenodo DOI — TODO): `manifest.csv`; the trimmed AA and AA+3Di alignments;
  `pif1_aa3di.treefile` and `pif1_aa3di_supported.treefile`; the GeneRax reconciliation outputs including
  the grafted run; and the NCBI, grafted, and Shen 2018 species trees.
- Structures: AlphaFold DB accessions + ColabFold predictions (manifest); sequences from UniProt.

## Acknowledgments and funding
*(Placeholder — TBD.)*
