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

### R4 — Reconciliation places the PIF1/RRM3 duplication on the Saccharomycotina branch
- Species-overlap reconciliation maps the duplication to the **Saccharomycotina** ancestor.
- GeneRax confirmation + loss pattern (which lineages retain both paralogs vs one) + stem-vs-crown:
  `[PENDING]`. *(Fig. 4 — duplication mapped onto the species tree, with losses.)*

### R5 — The two-copy mushrooms are an independent duplication
- Across all trees, Agaricomycetes (mushroom) PIF1 duplications form a separate event, not orthologs
  of yeast Pif1/Rrm3 — copy-number convergence, resolved by the tree. *(Fig. 4 inset / Fig. S3.)*

### R6 — Absolute timing — [PENDING]
- Map the Saccharomycotina branch onto a calibrated fungal time tree (bracket) + an independent
  relaxed-clock estimate (treePL/MCMCtree). *(Fig. 5.)*

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
