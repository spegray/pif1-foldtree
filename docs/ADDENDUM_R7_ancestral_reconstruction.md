# Addendum — R7 ancestral reconstruction (drop-in blocks)

*Everything new from the ancestral-reconstruction thread, packaged for merging into the working draft.
Each block is labeled with where it belongs. Underlying data: `results/asr/` (reconstruction,
`ANCESTOR_SUMMARY.md`; fold scores, `ancestor_fold_scores.tsv`).*

---

## New citation this section needs

R7 uses **AlphaFold3**, which is a different tool from the AlphaFold2/ColabFold already cited. Add:

> Abramson J, Adler J, Dunger J, Evans R, Green T, Pritzel A, et al. 2024. Accurate structure prediction of
> biomolecular interactions with AlphaFold 3. *Nature* 630(8016):493–500. https://doi.org/10.1038/s41586-024-07487-w

*(Verify the DOI/pages against the journal when you finalize the reference list.)*

---

## Abstract — insert this sentence

Place it right before the closing "The result is a concrete case…" sentence:

> Finally, reconstructing the pre-duplication ancestor shows that the single Saccharomycotina PIF1 that
> duplicated already carried the G4-engaging basic wedge, with a predicted structure placing that wedge on
> the quadruplex, so the motif predates the split.

---

## Introduction — two inserts

**(a) Roadmap sentence** — extend the "…reconcile it against a fungal species tree to read off the branch on
which the duplication occurred" clause to:

> …reconcile it against a fungal species tree to read off the branch on which the duplication occurred,
> **then reconstruct the ancestral gene that produced the pair** (Fig. 1).

**(b) Framing-paragraph closer** — append to "…a Saccharomycotina innovation, roughly as old as the
budding-yeast lineage itself":

> …roughly as old as the budding-yeast lineage itself, **and the ancestral gene that duplicated already
> carried the G4-engaging basic wedge that both paralogs kept.**

---

## Methods — new subsection

### Ancestral reconstruction and folding
- Marginal (empirical-Bayes) ancestral-sequence reconstruction with IQ-TREE (`-asr`;
  `workflow/19_asr_ancestor.py`), the AA+3Di topology held fixed and branch lengths plus LG+I+G4
  parameters re-optimized; the reconstructed node is the MRCA of ScPif1 and ScRrm3 (the 197-species
  Saccharomycotina clade), giving the single pre-duplication PIF1 core. Per-site posteriors and the
  wedge-column mapping in `results/asr/` (`extract_ancestor.R`).
- The reconstructed core (206 aa, after dropping three columns gapped across the clade) was folded with
  **AlphaFold3** (alphafoldserver.com), alone and on the same parallel-G4 substrate with ATP and two K⁺
  used for the extant family; G4 formation and the wedge-to-5′-tetrad distance were scored exactly as for
  the modern models (`results/asr/score_ancestor_fold.py`).

---

## Results — new section R7

### R7 — The pre-duplication ancestor is reconstructable, and it already carries the G4-engaging wedge
Having placed and dated the duplication, we can ask what the gene that duplicated looked like. Marginal
ancestral-sequence reconstruction on the AA+3Di tree, with the topology held fixed, recovers the helicase
core of the single pre-duplication PIF1 at the base of the Saccharomycotina clade, the node whose two
daughters are the Pif1 and Rrm3 lineages. The reconstruction is well-resolved: across the 209-column core
the mean per-site posterior is 0.94, with 89 per cent of sites at or above 0.8, because the node sits above
197 densely-sampled descendants rather than at the saturated depth that defeats the sequence tree. Only the
accessory N- and C-terminal domains, trimmed away because they cannot be aligned family-wide, are beyond
reach; the core, which carries the catalytic and DNA-engaging machinery, comes back cleanly.

That ancestral core already carries the G4-engaging wedge. The residue that in *S. cerevisiae* Pif1 is the
Arg324 wedge reconstructs as arginine with posterior 0.96, so the basic wedge predates the Pif1/Rrm3 split
rather than arising in the Pif1 lineage after it. Folding the reconstructed core with AlphaFold3, on the
same parallel-G4 substrate used for the extant family, bears this out in three dimensions: the quadruplex
forms in all five models, and the ancestral wedge reaches the 5′ G-tetrad at a median of 5.6 Å (best model
4.0 Å), indistinguishable from the modern PIF1 aggregate (median 5.5 Å); without the substrate the wedge
still folds as a confident, ordered residue (pLDDT 96). The gene that duplicated on the Saccharomycotina
stem was, on this evidence, already a G4-engaging PIF1 helicase with the wedge in place, a capability both
Pif1 and Rrm3 inherited (Fig. 6, Table S3; the wedge's conservation across the extant family is treated
separately).

---

## Table S3

**Table S3 — Ancestral-core fold on the G4 substrate** *(the R7 structural test; AlphaFold3, five models).*
The reconstructed pre-duplication wedge (core position 89, Arg) engages the 5′ G-tetrad as the modern PIF1s do.

| model | G4 formed | wedge→5′-tetrad (Å) | ipTM |
|---|---|---|---|
| 0 | yes | 5.57 | 0.78 |
| 1 | yes | 6.82 | 0.77 |
| 2 | yes | 4.52 | 0.77 |
| 3 | yes | 4.04 | 0.76 |
| 4 | yes | 10.02 | 0.72 |
| **all / median** | **5/5** | **5.57** | **0.77** |

*Apo (no substrate): the wedge folds as an ordered Arg in all five models (pLDDT 96). Modern PIF1 aggregate
for comparison (companion analysis): G4 folded in 94%, wedge→5′-tetrad median 5.5 Å.*

---

## Figure 6

Image: `docs/figures/fig6_ancestor_g4.png` (PyMOL; script `docs/figures/scripts/fig6_ancestor_g4.py`).

**Figure 6 — The reconstructed pre-duplication ancestor engages the G4 through the Arg-wedge.** Top-ranked
AlphaFold3 model of the reconstructed pre-duplication Saccharomycotina PIF1 core (grey cartoon) folded on
the parallel-G4 substrate. The quadruplex forms (orange; the two channel K⁺ ions in purple), and the
ancestral wedge, Arg89 (teal sticks), projects from the 1A subdomain to cap the 5′ G-tetrad (nearest
wedge-to-tetrad atom 5.6 Å in this model; the G4 forms in all five models, median wedge-to-tetrad 5.6 Å,
Table S3). The pose is contact-level, not atomic (AlphaFold3 places the quadruplex ~10 Å off its
crystallographic position), so the figure shows that the ancestral wedge reaches the 5′ face, not its exact
geometry.

---

## Reconstructed ancestral core sequence (for the supplement / data record)

`results/asr/ancestor_dupnode.fasta` — the pre-duplication Saccharomycotina PIF1 core (206 aa; MRCA of
ScPif1/ScRrm3). Arg-wedge at position 89.

```
TLSEEQQHVLDMVVQGKSIFFTGSAGTGKSVLLREIIKRLRKKYGPDSVAVTASTGLAACNIGGTTLHSFAGIGLGNESVEQLVKKIRRN
KKSRQRWRNTKVLIIDEISMIDGELFDKLDQIARKIRKNDKPFGGIQLVITGDFFQLPPVSKDNNQPAKFCFESESWKECIKHTIVLTQVF
RQKDNEFIDMLNEMRLGKLSPETEQ
```
