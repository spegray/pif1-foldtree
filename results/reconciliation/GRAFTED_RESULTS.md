# GeneRax reconciliation on the grafted (phylogenomic) species tree — R4

> **Provenance note.** Numbers were first transcribed from a results artifact
> (https://claude.ai/code/artifact/8df3c2f7-c498-4a26-820a-ccddd682e8e3) generated remotely on the
> Windows/WSL box on 2026-07-07. The **authoritative GeneRax output has since been recovered into the
> repo** at `results/reconciliation/aa3di_grafted/run/` (2026-07-08), and the headline numbers below were
> **verified directly against it**: `node_448` D=43 (tree-wide maximum), `node_578` D=23, totals D=476 /
> SL=3618 (gene-tree speciation nodes = 940−476 = 464). Only the NCBI-tree comparison run
> (`aa3di_ncbi/run/`) is still on the Windows box; its counts here are corroborated by the
> reviewer-robustness writeup (`results/reviewer/REVIEWER_ROBUSTNESS.md`).

## What was run
- **GeneRax 2.0.4**, `UndatedDL`, **EVAL** reconciliation (not SPR tree-search — intractable
  for this single 941-tip family).
- Gene tree: `results/seq_tree/pif1_aa3di.treefile` (the AA+3Di tree).
- Species tree: **grafted binary species tree**, 719 taxa (Shen 2018 budding-yeast topology
  grafted into the non-Saccharomycotina NCBI backbone; `workflow/16_graft_species_tree.py`,
  `data/species_tree/grafted_species.nwk`).
- 941 genes reconciled across 719 species.
- Only prep needed: resolve the gene tree's single 27-way *Candida*/*Debaryomyces* polytomy
  to make it binary (GeneRax requirement).

## Verdict — PASS
**The PIF1/RRM3 duplication maps to the Saccharomycotina common ancestor** on the grafted
phylogenomic tree — unchanged from the original NCBI-tree result. Robust to the species-tree
backbone.

## Evidence 1 — species-tree side
- Saccharomycotina ancestor = `node_448` (the 104-leaf budding-yeast clade) is the
  **top duplication node in the entire tree: 43 duplications** on the grafted tree
  (vs 9 on the NCBI tree).
- Independent **Agaricomycetes (mushroom) duplication** recovered identically: `node_578`,
  **23 duplications** in both runs.

## Evidence 2 — gene-tree-intrinsic (species-tree-independent; the one to trust)
- MRCA of ScPif1 (`P07271`) and ScRrm3 (`P38766`) subtends **197 genes across 103 species,
  all Saccharomycotina**.
- Splits cleanly into a **Pif1 clade** (98 genes, 9 families) and an **Rrm3 clade**
  (99 genes, 8 families).
- **8 budding-yeast families carry both paralogs** → the duplication predates their radiation:
  Saccharomycetaceae · Debaryomycetaceae · Pichiaceae · Metschnikowiaceae · Saccharomycodaceae
  · Saccharomycopsidaceae · Trichomonascaceae · Wickerhamomycetaceae.

## Grafted vs NCBI tree — same answer, both backbones
| Metric                        | Grafted tree | NCBI tree |
|-------------------------------|-------------:|----------:|
| Duplications (total)          |          476 |       464 |
| Losses — SL (total)           |        3,618 |     3,226 |
| Speciations (total)           |          464 |       476 |
| **Saccharomycotina-node dups**|       **43** |     **9** |
| Agaricomycetes-node dups      |           23 |        23 |

*Same gene tree and method throughout; only the species tree differs.*

## Loss pattern
Losses concentrate on the Saccharomycotina nodes in both runs — after the ancestral
duplication many budding-yeast lineages dropped one of the two paralogs (recurrent
single-paralog loss). The grafted tree carries ~12% more total losses (3,618 vs 3,226),
from its re-grafted yeast topology conflicting slightly more with the gene tree, but the
**spatial pattern is unchanged**.

## Caveat
The grafted tree **inflates** duplication/loss counts *within* the *Candida* clade, because
its 27-way polytomy had to be resolved arbitrarily — so per-node counts *there* are not
literal. The MRCA placement and overall loss pattern are robust, which is why the
gene-tree-intrinsic test (Evidence 2) is the one to trust.
