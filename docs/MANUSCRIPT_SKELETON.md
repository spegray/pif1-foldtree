# Dating the PIF1/RRM3 duplication in Ascomycota — manuscript skeleton

*Working scaffold. Methods are near-complete (the pipeline is done); Results carry the findings we
have, with `[PENDING]` for the GeneRax confirmation and absolute dating. Narrative sections (Abstract,
Intro, Discussion) are stubs — draft those in Spencer's voice at finalization.*

---

## Working title
*When did PIF1 and RRM3 split? Structural phylogenetics places the budding-yeast PIF1-family
duplication at the origin of Saccharomycotina.*

## Abstract — [STUB]
One paragraph: the question; that amino-acid phylogenetics alone cannot resolve the deep node; that
adding the 3Di structural alphabet resolves it; the answer (Saccharomycotina); the independent
mushroom duplication; (absolute date `[PENDING]`).

## Introduction — [STUB]
- PIF1-family 5′→3′ DNA helicases; *S. cerevisiae* has two paralogs, **Pif1** (mito + nuclear) and
  **Rrm3**, with distinct roles at telomeres/rDNA/replication-fork barriers.
- The evolutionary question: one ancestral fungal PIF1 → when/where did the PIF1+RRM3 pair arise?
- Why this is hard: a short, fast-evolving helicase core saturates at the relevant depth → classical
  sequence trees lose the signal. Motivates **structural phylogenetics** (Moi et al. 2025) — fold is
  conserved past sequence saturation.
- State the contribution: a fungi-wide, structure-aware gene tree + species-tree reconciliation.

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
- **Point estimate ≈ 400 Mya** (Devonian) — the budding-yeast common ancestor (BYCA) age (Shen et al.
  2018); the event is just crown-ward of it.
- **Bracket:** lower = BYCA (~400 Mya); upper = the **Saccharomycotina–Pezizomycotina divergence**
  (older; Ascomycota-wide clocks place major Ascomycota splits in the ~450–590 Mya range — read exact
  node ages off Shen 2018 / Mende et al. 2021 *Sci Adv* at finalization).
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

## Discussion — [STUB]
- Why structure succeeded where sequence failed (saturation; the 3Di information gain — 194 vs 180 sites).
- Biological reading: PIF1/RRM3 as a Saccharomycotina innovation; relation to the whole-genome-duplication era `[after dating]`.
- Independent recruitment of a second PIF1 in mushrooms — convergence.
- Limitations: gene-tree error at depth; NCBI vs phylogenomic backbone; distance-tree noise in FoldTree.

## Data and code availability
- Repo (scripts `01`–`14`, manifest, trees, alignments): github.com/spegray/pif1-foldtree.
- Structures: AlphaFold DB accessions + ColabFold predictions (manifest); sequences from UniProt.
