# Dating the PIF1/RRM3 duplication in budding yeasts (Saccharomycotina)

*Working manuscript draft. Methods and Results are complete, figure-backed, and written out in prose;
Abstract, Introduction, and Discussion are in Spencer's voice. Front/back matter (byline, affiliations,
funding, Zenodo DOI) is placeholder pending finalization.*

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

### R1 — Copy number points to a Saccharomycotina-restricted duplication but cannot prove one
Counting PIF1-family copies per genome already suggests where the answer lies. Among the Saccharomycotina
we sampled, 91 of 105 species carry two or more copies, while the Pezizomycotina (430 of 447) and the
Taphrinomycotina (all 11) are single-copy (Fig. 2; Table 1). Read at face value, the duplication looks
restricted to the budding yeasts. Copy number is a treacherous witness, though; many mushrooms
(Agaricomycotina) carry two copies as well, and a tally cannot tell whether a shared two-copy state
reflects one ancestral duplication or several convergent ones. That distinction is exactly what a gene
tree, reconciled against the species tree, is built to make.

### R2 — Amino-acid trees place the duplication at the base of Fungi
Under both a plain LG+I+G4 model and a profile-mixture LG+C20 model (fit by PMSF), the most recent common
ancestor of ScPif1 and ScRrm3 maps to the base of Fungi rather than to any one subphylum (Fig. 3A, left
panel). The reason is visible in the tree itself: the sequence analysis scatters the budding-yeast Pif1
and Rrm3 orthologs into two clusters far apart on the tree, so the smallest clade containing both anchors
takes in almost the whole sample and their common ancestor sits near the root.

### R2b — The amino-acid data significantly favor that deep placement, and the cause is long-branch attraction
The deep placement is not a soft polytomy to be broken either way. A constrained approximately-unbiased
test (10,000 RELL replicates, LG+I+G4), which asks only whether Saccharomycotina monophyly is compatible
with the sequence data, rejects it: the amino-acid alignment marginally but significantly prefers the deep
tree (AU *p* = 0.033, ΔlnL = 143). The branch lengths say why. The Saccharomycotina copies evolve roughly
twice as fast as the non-Saccharomycotina single-copy genes (mean terminal branch 0.15 and 0.19 versus
0.08 substitutions per site), and within them the Rrm3 clade runs 1.22 times faster than the Pif1 clade
(root-to-tip 1.07 versus 0.88; permutation *p* < 0.0001). Fast-evolving ingroup copies together with
distant outgroups is the textbook configuration for long-branch attraction, which pulls exactly this kind
of modest, artifactual preference toward a deep split (Fig. 3B). A model-free Tajima relative-rate test on
the anchor pair is non-significant, but with two sequences and 209 columns it is underpowered; the
confounded whole-tree comparison (ΔlnL ≈ 1891) mixes global tree differences with the node question and is
not interpretable, so we set it aside.

### R3 — Adding the 3Di structural signal overrides the sequence preference and resolves Saccharomycotina *(headline)*
Folding the 3Di structural alphabet into the same alignment moves the node. In the AA+3Di tree the
ScPif1/ScRrm3 ancestor collapses to Saccharomycotina, a clade of 197 genes across 103 species that takes
in essentially every budding yeast we sampled and excludes the non-Saccharomycotina fungi (Fig. 3A, right
panel). Adding structure to an otherwise-identical alignment overrides the significant-but-weak amino-acid
preference and pulls the duplication from kingdom-level down to a single subphylum, the resolution the
FoldTree approach was built to provide.

### R4 — Maximum-likelihood reconciliation places the duplication on the Saccharomycotina ancestor
Reconciling the AA+3Di gene tree against the species tree with GeneRax (UndatedDL) maps the PIF1/RRM3
duplication to node_448, the common ancestor of all 104 sampled Saccharomycotina, confirming the
species-overlap result by an independent and error-tolerant method (Fig. 4); the literal recPhyloXML
reconciliation, rendered on its own, agrees (Fig. S3). The placement survives the tree's own uncertainty.
Although the duplication node is weak under bootstrap resampling (UFBoot 29, SH-aLRT 98) and its Rrm3
daughter weaker still (UFBoot 16, against a well-supported Pif1 daughter at 90), reconciliation integrates
over that gene-tree error and recovers the same event regardless. It also sits where a Saccharomycotina-
ancestral origin predicts: both daughter clades draw members from 8 of the 9 sampled budding-yeast
families (only Dipodascaceae is Pif1-only), and the earliest-branching sample, the Lipomycetales, is
single-copy, placing the duplication on the Saccharomycotina stem near the crown (with the caveat that
these early-diverging lineages are thinly sampled).

The placement holds when we change the species tree and when we perturb the structures. On a grafted
phylogenomic backbone that splices the Shen 2018 budding-yeast topology into the NCBI tree (719 taxa),
node_448 is the single largest duplication node in the whole tree (D = 43); on the NCBI-taxonomy backbone
the same node still carries the event (D = 9). Three reruns of the AA+3Di gene tree reach the same
Saccharomycotina placement: dropping the 129 ColabFold models and keeping only AFDB structures (D = 37),
keeping only the high-confidence cores (pLDDT ≥ 87; D = 39), and estimating a GTR20 exchange matrix from
the data in place of the fixed 3Di matrix (D = 43), which together close the predictor-batch,
structure-quality, and 3Di-matrix robustness questions (Fig. S5). One honest qualification: on the grafted
tree the global event totals are large (D = 476, with roughly 3,618 speciation-losses), because
reconciliation absorbs residual deep gene-tree noise as many small events. The interpretable signal is the
anchor's placement on the Saccharomycotina node, not the totals, and ranking the per-node counts shows the
two real events standing clear of that tail (Fig. S4).

### R5 — The mushroom two-copy state is a separate duplication, and Rrm3 is the more-derived paralog
The recurrent two-copy state among mushrooms is a second, independent duplication rather than shared
ancestry with the yeast pair. Its node, node_578, carries the largest duplication count outside the
budding-yeast ancestor (D = 23, and the tree-wide maximum on the NCBI backbone), a case of copy-number
convergence that copy-counting cannot see but the reconciliation resolves (Fig. 4; ranked per-node counts
in Fig. S4). The two yeast paralogs are also not symmetric descendants of the ancestral gene. The
single-copy PIF1s outside Saccharomycotina are co-orthologs of both by descent, none of 536 nesting inside
either yeast clade, yet 442 of them sit closer to ScPif1 than to ScRrm3 by sequence (mean core identity
54.5 versus 52.3 per cent; *S. pombe* Pfh1, at 55.1 versus 53.2, among them), and the Rrm3 clade carries
the longer branch (R2b). We read this as Pif1 retaining more of the ancestral character while Rrm3 diverged
further into the more-derived role (Fig. 3B). The reading rests on raw sequence distance and on a
bootstrap-weak Rrm3 clade, so we hold it as an asymmetry worth naming rather than a settled fact.

### R6 — The duplication dates to the Devonian–Carboniferous, long before the whole-genome duplication
Because the duplication maps to the Saccharomycotina stem, we can read its age off a published,
fossil-calibrated time tree rather than re-estimate it. The Y1000+ molecular clock of Shen et al. (2018),
built from 332 budding-yeast genomes, already dates this node. Mapping the duplication's species set onto
that calibrated tree (54 of our species matched) places its ancestor at a node subtending 317 of 332 tips,
roughly 95 per cent of the budding yeasts, with a crown age of 326 million years and a stem age of 383
million years; the tree root, the budding-yeast common ancestor, sits at 404 (Fig. 5). The duplication
occupies that stem branch, which puts it at roughly 326 to 383 million years ago, in the Devonian and
Carboniferous, just crown-ward of the budding-yeast common ancestor. The earliest-diverging budding yeasts
fall outside the clade, consistent with their single-copy state, and the direct read off the calibrated
tree refines the borrowed "~400 Mya" estimate. The split is three to four times older than the
*Saccharomyces* whole-genome duplication (about 100 million years, a crown-Saccharomycetaceae event), and
that ordering fits the biology: Rrm3 is present across the budding yeasts, including pre-whole-genome-
duplication lineages such as *Kluyveromyces* and *Lachancea*, exactly as a Saccharomycotina-ancestral
origin predicts.

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
(R2b), which together we read as Pif1 retaining more of the ancestral character while Rrm3 took on the
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
