# PIF1/RRM3 reconciliation — reviewer robustness (R1, R2, R3)

*2026-07-07. Three requested robustness variants of the AA+3Di gene tree, each reconciled against the
grafted species tree to test whether the PIF1/RRM3 duplication still maps to the Saccharomycotina ancestor.*

## Reviewer goals

| | Variant | Tips | Purpose |
|---|---|---|---|
| **R1** | AA+3Di, AFDB-only | 828 | drop the 129 ColabFold-predicted structures — does including predicted models drive the result? |
| **R2** | AA+3Di, core-pLDDT ≥ 87 | 748 | keep only high-confidence structures — does structure quality drive the result? |
| **R3** | AA+3Di, GTR20 3Di matrix | 957 | replace the fixed Foldseek 3Di substitution matrix with an estimated GTR20 — model dependence? |

None had been run before: every prior AA+3Di tree is the full 957-tip set on the fixed 3Di matrix.

## Method

Each variant: build the AA+3Di ML gene tree (IQ-TREE 3.1.2, partition `AA=LG+G` / `TDi=<matrix>+G`,
`-alrt 1000`, no UFBoot), reconcile with **GeneRax 2.0.4, UndatedDL, EVAL** against the binary grafted
species tree (`grafted_species.binary.nwk`, 719 taxa), and read the placement two ways: the reconciliation's
per-species duplication counts (`_aps_parse_events.py`), and — as a cross-check — the raw gene-tree MRCA of
ScPif1 (P07271)/ScRrm3 (P38766) (`_aps_deepdive.py`).

- **R1** = the alignment subset by manifest `structure_source == AFDB` (828 tips).
- **R2** = the subset by **core-region** mean pLDDT ≥ 87 (748 tips). The manifest `mean_plddt` is *full-length*
  (median 64, only 4 ≥ 87); the reviewer's ≥ 87 is the helicase-core pLDDT (median ≈ 88). Core pLDDT was
  computed per structure over `core_from..core_to`: AFDB structures from the AlphaFold-DB confidence JSON
  (`AF-<acc>-F1-confidence_v6.json`; the AFDB `.cif` files are not local), ColabFold structures from the
  local PDB B-factor.
- **R3** = the same alignment with the 3Di partition model set to `GTR20+G` (`part_gtr20.nex`).

## Result — the placement is robust (see `reviewer_robustness.png`)

The reconciliation puts **Saccharomycotina as the top duplication node in every variant**, and recovers the
independent Agaricomycetes (mushroom) duplication in every variant:

| Variant | mapped genes | Dups (D) | Losses (SL) | **Saccharomycotina dups** (node_448) | Agaricomycetes dups (node_578) |
|---|---|---|---|---|---|
| Main (grafted) | 941 | 476 | 3,618 | **43 — top** | 23 |
| **R3** GTR20 3Di | 941 | 489 | 3,920 | **43 — top** | 31 |
| **R1** AFDB-only | 815 | 402 | 3,245 | **37 — top** | 24 |
| **R2** pLDDT ≥ 87 | 736 | 386 | 3,352 | **39 — top** | 20 |

The PIF1/RRM3 duplication maps to the Saccharomycotina common ancestor under the GTR20 matrix (R3), on the
AFDB-only set (R1), and on the high-pLDDT set (R2). The conclusion does not depend on the predicted
structures, on lower-confidence structures, or on the specific 3Di matrix.

## Honest nuance: reconciliation robust, raw two-anchor MRCA fragile

The raw gene-tree MRCA(ScPif1, ScRrm3) — a shortcut that reports the smallest clade containing both anchors —
is **clean for R2** (subtends 94 species, both daughter clades tidy Saccharomycotina, like the main tree's
103) but **broad for R1 and R3** (634 and 718 species): the Pif1 daughter stays clean Saccharomycotina while
ScRrm3's exact position slips, inflating the two-anchor clade. This is a fragility of the shortcut metric,
**not** evidence the duplication moved — the reconciliation (which accounts for the whole tree and its
losses) still concentrates the duplication on the Saccharomycotina node regardless of where the single
ScRrm3 tip sits, because the many two-paralog budding-yeast species force it there. The tell is R2, the
highest-quality subset, which is clean on both readouts.

The wobble traces to these being **ML-only trees (no UFBoot / `-bnni`)**: the duplication node was already
bootstrap-weak in the main analysis (UFBoot 29), so without bootstrap refinement ScRrm3's placement is
unstable in R1/R3. If the raw-tree reciprocal Pif1/Rrm3 monophyly is needed for R1/R3 (not just the
reconciliation placement), re-run those two with `-B 1000 -bnni` (~4.5 h each) to stabilize it.

## Caveats / method notes

- ML trees are bootstrap-free (`-alrt 1000`, no UFBoot) for speed; the reconciliation placement is robust,
  the raw-tree topology of the weakly-supported deep nodes is not (see above).
- GeneRax needs the mapping restricted to genes present in each subset tree; the standard prep script
  hard-codes the full 941-gene map, so R1/R2 required a per-variant `gene_species.map` subset before GeneRax
  would run (`_aps_fix_maps.py`).
- Reconciliation totals (D/SL) scale with tip count; the Saccharomycotina-node concentration is the
  comparable quantity, and it is stable at 37–43 across all conditions.

## Files (in `results/reviewer/` unless noted)

- `REVIEWER_ROBUSTNESS.md` — this document
- `reviewer_robustness.png` — the comparison figure
- `R3_gtr20.treefile`, `R1_afdb.treefile`, `R2_plddt87.treefile` — the three ML gene trees (+ `.iqtree` reports)
- `part_gtr20.nex` — the GTR20 3Di partition; `aln_afdb.fasta`, `aln_plddt87.fasta` — the R1/R2 alignment subsets
- `core_plddt.tsv` — per-tip core-region pLDDT (the R2 filter)
- `runall.log` — the full build/reconcile log
- `../reconciliation/rev_{R3_gtr20,R1_afdb,R2_plddt87}/run/` — the GeneRax reconciliations (eventCounts,
  speciesEventCounts, reconciled trees)
