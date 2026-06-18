# Dating the PIF1/RRM3 duplication in Ascomycota — manuscript skeleton

*Working scaffold. Methods are near-complete (the pipeline is done); Results carry the findings we
have, with `[PENDING]` for the GeneRax confirmation and absolute dating. Narrative sections (Abstract,
Intro, Discussion) are stubs — draft those in Spencer's voice at finalization.*

---

## Working title
*When did PIF1 and RRM3 split? Structural phylogenetics places the budding-yeast PIF1-family
duplication at the origin of Saccharomycotina.*

## Abstract — [DRAFT, Spencer's voice]
The budding yeast *Saccharomyces cerevisiae* carries two PIF1-family 5′→3′ DNA helicases, Pif1 and
Rrm3, with separable roles at telomeres, the rDNA, and stalled replication forks. When the single
ancestral fungal PIF1 gene duplicated to produce this pair has stayed unsettled, largely because the
part of the protein that can be aligned across the family, the helicase core, is short and fast-evolving,
and at the depth where the answer lies the amino-acid signal has saturated. We gathered the cellular
PIF1-family helicases of 728 fungal species (957 proteins; an Ascomycota-wide ingroup with
Basidiomycota, early-diverging fungi, and human PIF1 as outgroups), inferred a gene tree, and
reconciled it against a fungal species tree to locate the duplication. Amino-acid trees left the
duplication node unresolved, collapsing it to the base of Fungi; adding the 3Di structural alphabet to
the same alignment resolved it onto the stem of Saccharomycotina. Maximum-likelihood reconciliation
(GeneRax) confirms that placement, which dates the Pif1/Rrm3 split to the budding-yeast common ancestor
(~400 Mya, Devonian), well before the *Saccharomyces* whole-genome duplication. We further find that
Rrm3 is the more-derived of the two paralogs, and that the recurrent two-copy state among mushrooms is
an independent duplication rather than shared ancestry with the yeast pair. The result is a concrete
case where structural information recovers a deep relationship that sequence alone cannot.

## Introduction — [DRAFT, Spencer's voice]
The PIF1-family helicases are conserved 5′→3′ DNA helicases that act wherever the replication fork
meets trouble. In *Saccharomyces cerevisiae* the family has two members, Pif1 and Rrm3, and they have
divided the labor: Pif1 maintains the mitochondrial genome, unwinds G-quadruplexes, and negatively
regulates telomerase, while Rrm3 clears protein barriers ahead of the fork at the rDNA, the tRNA genes,
and other hard-to-replicate sites. Both descend from a single ancestral fungal PIF1, but where on the
fungal tree that ancestor duplicated has not been pinned down.

The question is harder than it looks. The region of these proteins that can be aligned across the whole
family is the helicase core, a few hundred residues that evolve quickly; over the distances that
separate the major fungal lineages the amino-acid signal in that core saturates, and a sequence tree
can no longer tell which arrangement of the deep branches is real. This is the regime where protein
structure helps: a fold is conserved long after the sequence that specifies it has been overwritten, so
a tree built from structural features can recover relationships that have gone dark at the sequence
level (Moi et al. 2025). The Pif1/Rrm3 split sits at exactly this awkward depth, which made it a good
candidate for a structure-aware approach.

Here we sample the cellular PIF1-family helicases across the fungi (957 proteins from 728 species, with
human and early-diverging-fungal outgroups), build the gene tree three ways (from sequence, from
structure, and from the two combined), and reconcile it against a fungal species tree to read off the
branch on which the duplication occurred. Sequence and structure disagree, and the disagreement is
itself informative: it points to a single, datable event. The duplication that produced Pif1 and Rrm3
is a Saccharomycotina innovation, roughly as old as the budding-yeast lineage itself.

---

## Methods

### Taxon sampling and homolog identification
- Fungi-wide sample: **957 cellular PIF1-family proteins across 728 species** — Ascomycota ingroup
  (681) + outgroups: Basidiomycota (235), Mucoromycota (36), Chytridiomycota (4), and human PIF1 (1).
- **Family filter — InterPro IPR048293 ("PIF1_RRM3_pfh1"), *not* Pfam PF05970.** PF05970 is also
  carried by **Helitron** rolling-circle transposons, inflating fungal counts 3–10× (e.g. Ascomycota
  2,255→681); IPR048293 isolates the *cellular* helicases while retaining every anchor. *(Key
  decision — Fig. S1 / Table S1.)*
- Anchors: *S. cerevisiae* Pif1 (UniProt **P07271**), Rrm3 (**P38766**); *S. pombe* Pfh1 (**Q9UUA2**),
  human PIF1 (**Q9H611**). Sequences retrieved from UniProt; one representative per gene.

### Structures
- **AlphaFold DB** for 828/957 proteins (model version v6); the **129** lacking an AFDB model were
  folded with **ColabFold (AlphaFold2)** to match the predictor (median pLDDT 62.5 vs 64.6 for AFDB —
  predictor-consistent). 3 very long proteins folded as their helicase-core region only (flagged in
  `manifest.csv`). Provenance per protein recorded in `manifest.csv`.

### Core trimming ("corecut")
- Both sequences and structures trimmed to the conserved PIF1 helicase core (PF05970 HMM envelope) so
  variable N/C-terminal extensions don't confound the trees. Structures sliced to the **same**
  residue coordinates as the sequence core (so the AA and structural analyses describe identical sites).

### Sequence maximum-likelihood gene tree
- Align cores (MAFFT) → trim (trimAl `-automated1`, **209 columns**) → **IQ-TREE LG+I+G4**, 957 taxa,
  1000 ultrafast bootstraps + SH-aLRT. *(Also LG+C20 via PMSF as a model-adequacy check.)*

### Structural gene tree (FoldTree)
- **foldseek** all-vs-all (3Di+AA, exhaustive) → FoldTree `fident` distance with the −b·ln(1−d/b)
  correction (b = 0.93) → **FastME**; `alntmscore` and `lddt` distances as metric-robustness checks.
  Rooted on human PIF1.

### AA+3Di partitioned tree *(the resolving method)*
- Per-residue **3Di** strings (foldseek) mapped onto the AA alignment columns → a 3Di alignment with
  identical gap structure → concatenated (**957 × 418**: 209 AA + 209 3Di). **IQ-TREE partition model:
  LG+G for the AA partition, the 3DiPhy 3Di substitution matrix (`3di_substmat`) for the 3Di partition.**
  The 3Di partition carried *more* parsimony-informative sites (194) than the AA core (180).

### Species tree
- **NCBI-taxonomy** backbone for the 728 species (ete3 `NCBITaxa`, pruned, forced binary; 720 leaves
  after dropping 8 taxa NCBI doesn't recognize → 941 mapped genes). *Rigorous backbone (Y1000+ / Li
  et al. 2021) pruned to our taxa: `[PENDING — confirmation pass]`.*

### Reconciliation = the dating step
- **Species-overlap** reconciliation (ete3) of each gene tree against the species tree, tallying
  duplications by the NCBI clade they map to.
- **GeneRax** (UndatedDL, `--unrooted-gene-tree`, EVAL + SPR) as the gold-standard ML
  duplication–loss reconciliation, on native x86 (WSL2). `[PENDING — Windows run]` → gives the
  duplication branch, the loss pattern, and stem-vs-crown.

### Reproducibility
- All steps scripted (`workflow/01`–`14`) in a version-pinned conda environment; every protein,
  structure source, and decision logged in `manifest.csv` and the README. Repo:
  github.com/spegray/pif1-foldtree. Key tool versions: IQ-TREE 3.1.2, foldseek 10, GeneRax 2.1.3,
  ete3 3.1.3, ColabFold/AF2.

---

## Results

### R1 — Copy number hints at, but cannot prove, a Saccharomycotina-restricted duplication
- Two-copy vs one-copy by subphylum: **Saccharomycotina 91/105 two-copy**; Pezizomycotina 430/447 and
  Taphrinomycotina 11/11 **one-copy**; two-copy mushrooms (Agaricomycotina) common too. Read: the
  duplication *looks* Saccharomycotina-restricted, but copy-counting can't distinguish shared ancestry
  from convergence — the gene tree must. *(Fig. 2 / Table 1.)*

### R2 — Amino-acid trees cannot resolve the deep node
- In both LG+I+G4 and LG+C20(PMSF) trees, the MRCA of ScPif1 and ScRrm3 maps to **"Fungi" (kingdom)** —
  the deep relationship is unresolved (209-aa core saturated). *(Fig. 3a.)*

### R3 — The 3Di structural signal resolves it → Saccharomycotina *(headline)*
- In the **AA+3Di** tree the ScPif1/ScRrm3 anchor collapses to **Saccharomycotina** (197 genes /
  103 species — i.e. essentially all sampled Saccharomycotina), with non-Saccharomycotina fungi
  *outside* the duplication clade. Adding structure to an otherwise-identical alignment moves the
  node from kingdom-level to subphylum-level — a clean demonstration of the FoldTree thesis.
  *(Fig. 3b — the key figure; AA-vs-AA+3Di side by side.)*

### R4 — GeneRax reconciliation places the duplication on the Saccharomycotina ancestor
- Gold-standard ML reconciliation (GeneRax 2.0.4, `UndatedDL`, `--unrooted-gene-tree`) of the AA+3Di
  gene tree maps the PIF1/RRM3 duplication to **`node_448` — the MRCA of all 104 sampled
  Saccharomycotina** (the budding-yeast common ancestor) — confirming the species-overlap result with
  an independent, error-tolerant method.
- **Robust to gene-tree wobble:** the same placement holds on the bootstrap-supported AA+3Di tree;
  although the duplication node is UFBoot-weak (29; SH-aLRT 98 — expected for a deep, saturated split),
  reconciliation (which integrates over gene-tree error) places it at Saccharomycotina regardless.
- **Stem, just crown-ward:** both daughter (Pif1, Rrm3) clades contain members of **8/9 sampled
  Saccharomycotina families** (only Dipodascaceae is Pif1-only); the earliest-branching sample
  (Lipomycetales) is single-copy → the duplication sits on the Saccharomycotina stem near the crown
  (caveat: sparse early-diverging-lineage sampling).
- *Honest note:* the global event totals are high (D ≈ 464; ~3,226 speciation-losses) because
  reconciliation absorbs residual deep gene-tree noise as many small D/L events; the interpretable
  signal is the anchor's placement on the Saccharomycotina node, not the global counts.
  *(Fig. 4 — duplication + losses mapped onto the species tree.)*

### R5 — An independent mushroom duplication, and RRM3 as the derived paralog
- The two-copy **mushrooms (Agaricomycetes) are a separate duplication**: `node_578` carries the single
  largest per-node duplication count (**D = 23**), distinct from the yeast Pif1/Rrm3 event — copy-number
  convergence, invisible to copy-counting but resolved by the tree. *(Fig. 4 inset / Fig. S3.)*
- **RRM3 is the more-derived paralog.** Single-copy PIF1s outside Saccharomycotina are co-orthologs of
  *both* yeast paralogs by descent (0/536 nest inside either clade), but lean toward Pif1 by sequence
  (442/536 closer to ScPif1; mean core identity 54.5% vs 52.3%; *S. pombe* Pfh1 55.1 vs 53.2) — i.e.
  Pif1 retained more of the ancestral state while RRM3 diverged further.

### R6 — Absolute timing: a Devonian (~400-Myr-old) duplication, long predating the WGD
**Approach — read the age off a published, fossil-calibrated time tree rather than re-estimate it.**
The duplication maps to the Saccharomycotina stem, and the Y1000+ budding-yeast molecular clock
(Shen et al. 2018, *Cell*; 332 genomes) already dates exactly this node:
- **Direct placement (this study, `workflow/17_place_on_timetree.py`):** mapping the duplication's
  species set onto Shen's calibrated tree puts its MRCA at a node subtending ~96% of budding yeasts
  (317/330 tips, 54 of our species matched), dated **crown 326 Mya, stem 383 Mya** (tree root / BYCA
  404 Mya). The duplication sits on that stem branch → **~330–383 Mya (Devonian–Carboniferous), in the
  early Saccharomycotina crown, just crown-ward of the ~404-Mya budding-yeast common ancestor.** The
  earliest-diverging budding yeasts fall outside the clade, consistent with their single-copy state.
  This refines (and supersedes) the borrowed "~400 Mya" estimate with a direct read off the calibrated
  phylogenomic tree.
- **Independent sanity check:** ~3–4× older than the *Saccharomyces* **whole-genome duplication
  (~100 Mya**, a crown-Saccharomycetaceae event) — consistent with RRM3 occurring across budding yeasts
  *including* pre-WGD lineages (*Kluyveromyces*, *Lachancea*). The ordering (duplication ≫ WGD) is robust
  to the exact calibration.
- **Optional own relaxed clock** (only if a de-novo estimate is wanted): treePL or MCMCtree on a
  calibrated species tree with secondary calibrations from Shen 2018 — but Y1000+ already did this
  rigorously for this subphylum, so the marginal value is low. *(Fig. 5 — calibrated placement.)*

---

## Figures
1. **Pipeline schematic** — homologs → structures → corecut → {AA tree, FoldTree, AA+3Di} → species tree → reconciliation.
2. **Sampling + copy number** — taxa across the fungal tree; two-copy vs one-copy by subphylum.
3. **The resolution** — (a) AA-only gene tree: duplication node maps to "Fungi"; (b) AA+3Di: maps to Saccharomycotina. *(headline)*
4. **Reconciliation** — duplication + losses mapped onto the species tree (the answer). `[PENDING GeneRax]`
5. **Absolute date** — calibrated placement. `[PENDING]`
- **S1** Family-filter validation (IPR048293 vs PF05970 / Helitron contamination).
- **S2** Structure-source + pLDDT (AFDB vs ColabFold consistency).
- **S3** Per-metric structural trees (fident / alntmscore / lddt) + support values.
- **S4** Sequence-vs-structure concordance.

## Discussion — [DRAFT, Spencer's voice]
The practical lesson here is narrow but clean: for this duplication the amino-acid core had simply run
out of signal, and better modeling did not recover it (a profile-mixture model fit on the same 209
columns left the node exactly where the plain model did). What moved the node was information of a
different kind. The 3Di structural alphabet carried more resolving power than the amino-acid core over
the same sites (194 versus 180 parsimony-informative positions), and folding that signal into the tree
pulled the Pif1/Rrm3 ancestor from the base of Fungi down onto the Saccharomycotina stem.

Placed in time, the duplication is old. Mapping the Saccharomycotina stem onto the calibrated Y1000+
time tree (Shen et al. 2018) puts the split near the budding-yeast common ancestor, on the order of 400
million years ago. That timing matters for how the pair is usually explained: the *Saccharomyces*
whole-genome duplication, the first event people tend to reach for when accounting for yeast gene pairs,
is roughly four times younger and cannot have produced Pif1 and Rrm3. Rrm3 is present across the budding
yeasts, including lineages (*Kluyveromyces*, *Lachancea*) that diverged well before the whole-genome
duplication, exactly as a Saccharomycotina-ancestral origin predicts.

The two paralogs are also not interchangeable copies of the ancestral gene. Single-copy PIF1s outside
Saccharomycotina are co-orthologs of both by descent, yet they sit consistently closer to Pif1 than to
Rrm3 in sequence (442 of 536 such proteins, *Schizosaccharomyces pombe* Pfh1 among them), which we read
as Pif1 retaining more of the ancestral character while Rrm3 took on the faster-evolving, more-derived
role. Whether that asymmetry tracks the functional specialization of the two helicases is a question
this tree cannot answer, but it is a natural one to ask next.

Copy number is a poor guide to history. Many mushrooms (Agaricomycetes) also carry two PIF1-family
genes, and on counts alone that resembles the budding-yeast situation; the gene tree shows it is a
separate, later duplication. Convergent gene duplication of this kind is easy to mistake for shared
ancestry, and telling the two apart is precisely what a reconciliation, rather than a tally, is built
to do.

Two caveats temper the confidence and mark the work's edges. First, the duplication node is well
supported by the SH-aLRT test but weak under bootstrap resampling (UFBoot 29), as expected for a deep
node resting on saturated sequence; the result leans on reconciliation, which tolerates that gene-tree
uncertainty, more than on any single branch's bootstrap value. Second, the early-diverging budding
yeasts are thinly sampled here, so while the duplication maps cleanly to the Saccharomycotina stem,
whether it falls just before or just after the subphylum's first split is not yet settled; denser
sampling of the Lipomycetales and their relatives, with the published phylogenomic backbone in place of
the NCBI topology, would tighten that. Neither caveat unseats the central result, and both read as the
next steps rather than reasons for doubt.

## Data and code availability
- Repo (scripts `01`–`14`, manifest, trees, alignments): github.com/spegray/pif1-foldtree.
- Structures: AlphaFold DB accessions + ColabFold predictions (manifest); sequences from UniProt.
