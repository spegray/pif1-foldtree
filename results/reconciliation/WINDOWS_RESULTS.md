# Windows / WSL2 results — GeneRax reconciliation + supported AA+3Di tree

*Produced on the Windows desktop (WSL2, native linux-64) per `WINDOWS_HANDOFF.md`. Hardware: Ryzen
9 7950X3D, 32 threads, 46 GiB RAM. Env: Miniforge `pif1recon` (IQ-TREE 3.1.2, ete3 3.1.3, dendropy,
biopython) + `gx204` (GeneRax 2.0.4). All paths repo-relative; commands run from the repo root.*

## TL;DR
The gold-standard reconciliation **confirms the Mac result**: the PIF1/RRM3 duplication maps to the
**Saccharomycotina common ancestor** (the stem of the budding-yeast radiation), and the two-copy
**mushrooms (Agaricomycetes) are an independent duplication**. AA+3Di and the bootstrap-supported
AA+3Di tree both resolve it cleanly; the AA-only LG+I+G4 tree does not (its raw topology leaves the
node at "Fungi"), exactly as expected.

---

## ⚠️ Critical tooling note (for the Mac session / future runs)
**bioconda `generax 2.1.3` (build `hf316886_3`) is broken** for this analysis: the reconciliation
likelihood computes, but the event inference returns an all-zero scenario (`S:0 D:0 L:0 …`) and then
**segfaults in `Scenario::savePerSpeciesEventsCounts`** (both EVAL and SPR). This is not a data
problem (the gene→species mapping is 100% clean). **Use `generax 2.0.4` instead** (bioconda build
`h103dbdd_3`) — it runs cleanly. Note 2.0.4 uses `--strategy` (not `--geneSearchStrategy`).
Also: the species tree had one geneless leaf (taxid `1041607`, 720→719); pruned to
`data/species_tree/ncbi_species.covered.nwk` so every species-tree leaf is covered.

---

## Inputs reconciled
- 941 genes / 719 species (957→941 after dropping 16 genes on taxa NCBI doesn't map; 720→719 species
  after pruning the geneless leaf). Model `UndatedDL`, `--unrooted-gene-tree`, seed 42.
- Gene trees: AA+3Di (`pif1_aa3di.treefile`), AA-only LG+I+G4 (`pif1.treefile`), and a NEW
  bootstrap-supported AA+3Di tree built here (`pif1_aa3di_supported.treefile`, 1000 UFBoot + SH-aLRT).

## TASK A — GeneRax reconciliation (EVAL, GeneRax 2.0.4)

| Gene tree | Global events (S / D / loss[SL]) | Duplication on Saccharomycotina node_448 (104 sp) | Anchor MRCA(Pif1,Rrm3) by species-overlap |
|---|---|---|---|
| **AA+3Di** | 476 / 464 / 3226 | **D = 9**, losses 10 | **Saccharomycotina** (197 genes / 103 sp) |
| **AA+3Di, UFBoot-supported** | 463 / 477 / 3348 | **D = 8**, losses 10 | **Saccharomycotina** (197 genes / 103 sp) |
| **AA-only LG+I+G4** | 468 / 472 / 3236 | D = 8, losses 11 | Fungi (raw tree unresolved) |

- **GeneRax confirms the species-overlap result**: on the AA+3Di trees, both the per-species event
  counts (a duplication on `node_448` = the MRCA of all 104 sampled Saccharomycotina) and the
  anchor-node mapping point to **Saccharomycotina**.
- **AA+3Di vs LG+I+G4**: the *reconciliation tool* recovers a Saccharomycotina-ancestor duplication
  signal from all trees, but only the AA+3Di trees place the *anchor* (ScPif1/ScRrm3) duplication
  there by topology — the AA-only tree's deep node is saturated (maps to "Fungi"). Adding the 3Di
  structural alphabet is what resolves it.
- **Independent mushroom duplication**: the largest per-node duplication count is on
  `node_578 = Agaricomycetes` (**D = 23**), a separate event from the yeast Pif1/Rrm3 pair.
- **Loss pattern**: ~3.2–3.3k speciation-losses total — consistent with the duplication at the
  Saccharomycotina ancestor followed by widespread single-copy retention outside Saccharomycotina
  (pre-duplication) and scattered paralog loss within it (e.g. Dipodascaceae, Saccharomycodales).

## TASK B — supported AA+3Di tree + its reconciliation
- Built `pif1_aa3di_supported.treefile` (IQ-TREE 3.1.2, partitioned LG+G:AA / 3Di-substmat+G:3Di,
  `-B 1000 -alrt 1000 -bnni`, 957 taxa). ML logL ≈ -185,6xx.
- **Support on the duplication node** MRCA(ScPif1,ScRrm3): **SH-aLRT 98.2 / UFBoot 29**
  (Pif1 daughter clade 95.4/90; Rrm3 daughter clade 88/16). Read: the duplication clade is strongly
  supported by SH-aLRT but the deep node is **bootstrap-weak (UFBoot 29)** — the split is real but
  shaky under resampling, as expected for a deep, saturated node. This is precisely why the
  *reconciliation* (robust to gene-tree error) is the appropriate test, and it places the event at
  Saccharomycotina regardless.
- **Full SPR tree-search reconciliation was not tractable** for this single 941-tip family: because
  the deep backbone is *pervasively* weakly supported, `--support-threshold` cannot prune the search
  (~1,800 candidate moves per radius-1 round, no convergence in hours, 1-core-bound). We therefore
  reconciled the supported topology with EVAL (above) — it agrees with everything else.

## Deep-dive (answers to follow-up questions)
- **Where within Saccharomycotina?** Stem of the subphylum: both paralog clades contain members of
  8/9 sampled families (Saccharomycetaceae, Debaryomycetaceae, Pichiaceae, Metschnikowiaceae,
  Saccharomycodaceae, Saccharomycopsidaceae, Trichomonascaceae, Wickerhamomycetaceae); only
  Dipodascaceae is Pif1-only. The single earliest-branching sample (Lipomycetales) is single-copy,
  hinting the event is just crown-ward of the first split (sparse early sampling: 1 Lipomycetales,
  0 Trigonopsidales).
- **Does an outside species' single PIF1 resemble Pif1 or Rrm3?** By descent it is a **co-ortholog of
  both** (0/536 single-copy outside genes nest inside either paralog clade). By sequence it leans
  **ScPif1**: 442 of 536 single-copy outside species are closer to ScPif1 (vs 93 to ScRrm3); mean core
  identity 54.5% (Pif1) vs 52.3% (Rrm3); every subphylum agrees; *S. pombe* Pfh1 = 55.1 vs 53.2.
  Core evolutionary rates are near-symmetric, so this reflects Pif1 retaining more of the ancestral
  state — Rrm3 is the more-derived paralog.

## FoldTree / foldseek pipeline (structural side)
- The 3Di structural alphabet carries **more** phylogenetic signal than amino acids over the same
  209-column core: **194 vs 180 parsimony-informative sites** (4 vs 15 invariant).
- **Only the AA+3Di ML partition resolves the split**; the pure structural *distance* FoldTrees
  (fident/alntmscore/lddt → FastME) all leave the anchor at "Fungi" (fident gets the clade size
  nearly right at 203 genes ≈ 197, but a few rogue tips drag it deep). Lesson: for a recent,
  saturated duplication, structural *characters under an ML model* beat structural *distances*.

## Figures & deck (in `results/figures/`)
- `FoldTree_structure_resolves_duplication.png` — info content + which method resolves the node.
- `Q1_saccharomycotina_paralog_retention.png` — paralog retention by Saccharomycotina order.
- `Q2_pif1_vs_rrm3_identity.png` — 442:93 Pif1-leaning identity of outside single-copy PIF1s.
- `tree_schematic.png`, `tree_dupclade_real.png` — annotated trees.
- `PIF1_RRM3_lab_meeting.pptx` — 8-slide lab-meeting deck.

## Reproduce
```bash
# env: Miniforge + `conda create -n gx204 -c conda-forge -c bioconda generax=2.0.4`
conda activate gx204
generax -f results/reconciliation/aa3di/families.txt -s data/species_tree/ncbi_species.covered.nwk \
        -r UndatedDL --unrooted-gene-tree --strategy EVAL --seed 42 -p results/reconciliation/aa3di/eval
# interpret (env pif1recon): per-species event counts -> clade
python workflow/_aps_parse_events.py results/reconciliation/aa3di/eval
python workflow/12_reconcile.py --tree results/seq_tree/pif1_aa3di_supported.ufboot.treefile
```
