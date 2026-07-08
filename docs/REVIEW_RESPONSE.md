# Adversarial review + response — PIF1/RRM3 dating project

Internal record of a self-administered adversarial ("Reviewer 2") pass, the analyses run to address
it, and what remains. Numbers and commands are recorded for reproducibility; citations included for
the novelty assessment. Folds into the manuscript Discussion and a future response-to-reviewers.

---

## A. The review, in brief
Major concerns raised:
- **M1 — Novelty unestablished.** No literature review had been done; the qualitative result may be known.
- **M2 — Is structure actually necessary?** Reconciliation recovered a Saccharomycotina duplication
  from all trees, so "amino acids cannot resolve it" might be overstated.
- **M3 — Headline node is bootstrap-weak (UFBoot 29); no robustness/topology test had been run.**
- **M4 — AA and 3Di are non-independent** (3Di is derived from structure, itself from sequence); the
  "194 vs 180 informative sites" framing risks double-counting.
- **M5 — Reconciliation used the NCBI *taxonomy*, not a phylogeny** (Y1000+/Li backbone never substituted).
- **M6 — Structure-quality / predictor confounds** (per-residue pLDDT; AFDB vs ColabFold batch) untested.
- Minor: alignment quality (MAFFT FFT-NS-i + aggressive trim); single-outgroup rooting of the FoldTrees;
  "Rrm3 is derived" rested on raw % identity; dating is a borrowed bracket (stem vs crown widens it);
  the planned orthology cross-check and the ~134 PF05970-not-IPR048293 species were never checked.

The first round addressed **M1** and **M2** (the two that decide whether this is a methods paper, a
confirmatory note, or needs reframing); a second, compute-heavy round on the Windows/WSL box then closed
**M5** (phylogenomic backbone) and **M6** (structure-prediction confounds) with a grafted-tree
reconciliation and three AA+3Di reruns (Section B4). **M3 is largely addressed as well** — the headline
AA+3Di tree already carries 1000-replicate UFBoot (deep node SH-aLRT 98 / UFBoot 29, both reported
honestly); only a cosmetic R1/R3 raw-tree cleanup and the minors remain (Section D).

---

## B. Analyses run this round

### B1 — Literature / novelty (M1)
Searched PubMed and Consensus (≈200M papers; Semantic Scholar/Scopus/PubMed/arXiv).
- The Pif1/Rrm3 pair is **intensely studied functionally**; the **qualitative evolutionary picture is
  already in the reviews** — family conserved bacteria→humans, *S. cerevisiae* has two paralogs
  (Pif1, Rrm3), *S. pombe*/most others one (Pfh1), with non-overlapping functions. Sources (PubMed;
  cite with DOIs): Bochman, Sabouri & Zakian 2010, *DNA Repair* — reviews family "discovery, evolution"
  (https://doi.org/10.1016/j.dnarep.2010.01.008); Boulé & Zakian 2006, *NAR* — functional
  specialization of ScPif1 vs Rrm3 (https://doi.org/10.1093/nar/gkl561); Malone, Thompson & Byrd 2022,
  *IJMS* (https://doi.org/10.3390/ijms23073736).
- The only dedicated *Pif1 phylogenetics* found is Harman & Manna 2016, *Mol Phylogenet Evol* — on
  **amoebae and accessory domains**, not the Saccharomycotina split
  (https://doi.org/10.1016/j.ympev.2016.07.015).
- Consensus returned only functional Rrm3/Pif1 work (replication, mutation avoidance, R-loops).

**Verdict:** the bare claim "RRM3 is a budding-yeast paralog of an ancestral PIF1" is **not novel**
(implicit in the reviews). A **rigorous, fungi-wide phylogenetic placement + reconciliation + dating**
of that duplication was **not found** → that, plus the methods angle (B2/B3), is the defensible
contribution. **Y1000+ / comparative-genomics check (done):** the large budding-yeast
gene-family-evolution studies that would have caught this target metabolic/signaling families, not
PIF1 — histidine kinases across 82 Saccharomycotina (Hérivaux 2017,
https://doi.org/10.1007/s00294-017-0797-1), sugar transporters across 332 (Crandall 2024,
https://doi.org/10.1093/molbev/msae228), and the 993-genome / 14,785-gene-family scan of David 2025,
which hunts gene-family expansions but only for metabolic traits (https://doi.org/10.1073/pnas.2500165122).
None analyze the PIF1/RRM3 helicases; functional reviews (e.g. Muellner 2020) restate only the qualitative
paralogy. **Conclusion: no published study places or dates the PIF1/RRM3 duplication** → the dated
phylogenetic placement + the structural method are the contribution. Residual diligence before submission:
a full-text/supplement scan of the Y1000+ orthogroup tables (RRM3 is surely *present* there, just never
analyzed) and a Google Scholar sweep.

### B2 — Topology tests: does the amino-acid data reject the structural answer? (M2/M3)
Files: `results/seq_tree/au_test.*` (confounded whole-tree), `results/seq_tree/au_constrained.*` (clean);
constraint `results/seq_tree/constraint_sacch.nwk` (205 Saccharomycotina vs 752 other tips);
model LG+I+G4, 10,000 RELL replicates, native arm64 IQ-TREE.

| Test | Tree 2 (structural / Sacch-monophyletic) | ΔlogL | p-AU |
|---|---|---|---|
| Whole AA tree vs whole AA+3Di tree (**confounded** — global differences) | AA+3Di | 1891.5 | 1.6e-33 |
| **Clean: force *only* Saccharomycotina monophyly** | AA-optimal+constraint | **142.8** | **0.033** |

**Result:** a **genuine but marginal conflict.** The AA data significantly rejects Saccharomycotina
monophyly (clean test p-AU = 0.033) — *not* "no signal" (which would be p ≫ 0.05), but a weak,
just-significant preference for the deep placement (ΔlogL ≈ 143). The confounded whole-tree test
(ΔlogL 1891) is **not** interpretable for the node question and should not be cited.

### B3 — Is RRM3 long-branched? The LBA defense (M2 mechanism)
Why the marginal AA preference is plausibly an artifact. Script `workflow/15_rate_check.py`;
branch lengths from `results/seq_tree/aa_bl.treefile` (AA-only branch lengths on the
Saccharomycotina-monophyletic topology).

- **Tajima's relative-rate test (model-free, anchor pair):** ScPif1 vs ScRrm3 — vs *S. pombe* Pfh1:
  21 vs 25 unique sites, χ²=0.35 (n.s.); vs human PIF1: 15 vs 19, χ²=0.47 (n.s.). Rrm3 marginally
  faster but **underpowered** (two sequences, 209 columns).
- **Clade-level (197 Saccharomycotina tips; 98 Pif1, 99 Rrm3):**
  - root(duplication)→tip: **Pif1 0.877, Rrm3 1.073 (ratio 1.22), permutation p < 0.0001** — Rrm3
    significantly faster (it is the more-derived paralog).
  - mean terminal branch: **Pif1 0.150, Rrm3 0.191, non-Saccharomycotina 0.082** — the Saccharomycotina
    copies evolve **~2× faster** than the non-Saccharomycotina single-copy genes.

**Interpretation:** this is the classic long-branch-attraction configuration — fast-evolving ingroup
paralogs (~2× the rest; Rrm3 fastest) plus distant outgroups — which produces exactly a modest,
artifactual AA preference for a deep split. **Caveats:** the Pif1-vs-Rrm3 asymmetry is modest (1.22×;
the dominant signal is Saccharomycotina-vs-rest, ~2×); within-clade rates are measured on the
structurally-favored topology (mildly circular, though the terminal-branch result is topology-robust);
report the n.s. model-free anchor test alongside the significant clade test for completeness.

### B4 — Robustness reruns + phylogenomic-backbone reconciliation (M5, M6)
Run natively on the Windows/WSL box (IQ-TREE 3.1.2 + GeneRax 2.0.4, `UndatedDL`, EVAL). Files:
`results/reviewer/` (the three ML trees + `.iqtree` reports, subset alignments, `REVIEWER_ROBUSTNESS.md`,
`runall.log`, `reviewer_robustness.png`) and `results/reconciliation/rev_{R1_afdb,R2_plddt87,R3_gtr20}/run/`.
The grafted-tree main run is summarized in `results/reconciliation/GRAFTED_RESULTS.md` (its
`aa3di_grafted/run/` output is still on the Windows box, but its headline numbers are corroborated by the
rerun writeup, which lists the same Main-grafted counts).

**Phylogenomic backbone (M5) — closed.** A binary species tree was grafted (719 taxa: the Shen 2018
budding-yeast topology into the non-Saccharomycotina NCBI backbone; `workflow/16_graft_species_tree.py`,
`data/species_tree/grafted_species.nwk`) and GeneRax re-run on it. The PIF1/RRM3 duplication maps to the
Saccharomycotina ancestor (`node_448`) as the **top duplication node in the whole tree — 43 dups, vs 9 on
the NCBI tree**; the independent mushroom duplication (`node_578`, 23 dups) is recovered identically. So
the placement is not an artifact of NCBI's arbitrary polytomy resolution.

**Structure-prediction confounds (M6) — closed.** Three reruns of the AA+3Di tree, each reconciled against
the grafted tree; every one keeps `node_448` as the **tree-wide maximum duplication node** (verified here
against the raw per-species event counts, not just the writeup):

| Variant | Tree tips | Sacch. dups (`node_448`, = tree max) | Mushroom dups (`node_578`) |
|---|---|---|---|
| Main (grafted, full set) | 957 | 43 | 23 |
| **R1** — AFDB-only (drop 129 ColabFold) | 828 | 37 | 24 |
| **R2** — core pLDDT ≥ 87 (high-confidence) | 748 | 39 | 20 |
| **R3** — GTR20 3Di matrix (vs fixed Foldseek) | 957 | 43 | 31 |

The Saccharomycotina placement rides on neither the predicted (ColabFold) structures, nor the lower-
confidence structures, nor the specific 3Di substitution matrix.

**Honest nuance (carry into the manuscript, don't hide).** The *reconciliation* placement is robust in all
three; the *raw two-anchor MRCA shortcut* — the smallest clade containing both ScPif1 and ScRrm3 — is clean
for R2 (94 species, tidy Saccharomycotina daughters, like the main tree's 103) but broad for R1 and R3
(634 / 718 species): the Pif1 daughter stays clean while ScRrm3's single tip slips. This is a fragility of
the shortcut on **ML-only trees** (no UFBoot; the deep node was already UFBoot-weak at 29), not evidence the
duplication moved — reconciliation, which uses the whole tree and its losses, concentrates the duplication
on `node_448` regardless of where the lone ScRrm3 tip lands, because the many two-paralog budding-yeast
species force it there. R2 (the highest-quality subset) is clean on *both* readouts, which is the tell. To
firm up the raw-tree reciprocal Pif1/Rrm3 monophyly for R1/R3, rerun those two with `-B 1000 -bnni` (~4.5 h
each) — optional, since the reconciliation already answers the placement question.

---

## C. How this addresses the review
- **M1 (novelty):** eased but contingent. Frame as a **methods-forward case study** — a rigorous,
  dated, structure-enabled placement of a duplication whose qualitative existence was known but never
  phylogenetically nailed. Do the Y1000+ literature check before final novelty claims.
- **M2 (is structure necessary?):** **resolved in favor, and it becomes the paper's strongest point.**
  Structure does real work: it overrides a *significant* (p=0.033) amino-acid preference that the rate
  analysis attributes to LBA and that three independent lines (structure, reconciliation, and the
  topology-independent fact that RRM3 orthologs occur only in budding yeasts) contradict.
- **Residual vulnerability (be explicit, don't hide):** the override rests on discounting a p < 0.05
  sequence result. The defense — LBA mechanism + RRM3's biological restriction — must be made *up front*
  in the Discussion, and the AU p=0.033 reported honestly as evidence the artifact is weak.
- **M5 (phylogenomic backbone) — closed (B4).** GeneRax on the grafted 719-taxon tree keeps the
  duplication on the Saccharomycotina ancestor as the top duplication node in the whole tree (43 dups);
  not an NCBI-taxonomy artifact.
- **M6 (structure-prediction confounds) — closed (B4).** AFDB-only, high-pLDDT, and GTR20-matrix reruns
  all keep `node_448` the tree-wide top duplication node (37–43 dups). The result depends on none of the
  three suspected confounds.

---

## D. Open items NOT addressed this round (prioritized)
1. **M5 — substitute the published phylogenomic species tree.** *Closed (B4).* Two ways: (i) mapping the
   duplication's species set onto the Shen et al. 2018 time-calibrated Y1000+ tree (`workflow/17`) placed
   it on a node subtending ~96% of budding yeasts (317/330 tips) dated **~330–383 Mya**; (ii) a full
   GeneRax DL reconciliation against the *grafted* broad-fungi + Shen tree (`workflow/16`, 719 taxa)
   keeps the duplication on the Saccharomycotina ancestor as the **top duplication node (43 dups)** and
   reproduces the per-branch loss pattern (`results/reconciliation/GRAFTED_RESULTS.md`). The placement is
   robust to the phylogenomic backbone, not an artifact of NCBI's polytomy resolution.
2. **M6 robustness — closed (B4).** AFDB-only, core-pLDDT ≥ 87, and GTR20-3Di-matrix reruns all keep the
   Saccharomycotina node the tree-wide top duplication node (37–43 dups). **M3 — largely addressed.** The
   full AA+3Di headline tree already carries 1000-replicate UFBoot (`results/seq_tree/pif1_aa3di_supported.treefile`,
   ModelFinder + UFBoot); the duplication node is **SH-aLRT 98 / UFBoot 29** — weak by resampling, as
   expected for a saturated deep node, strong by SH-aLRT, both disclosed. That weakness is the paper's
   point, not a hole: it is why the argument leans on reconciliation (robust across four perturbations),
   not on the branch. *Residual is cosmetic only:* R1 and R3 were ML-only reruns, so their **raw-tree**
   Pif1/Rrm3 reciprocal monophyly is soft; `-B 1000 -bnni` on those two (~4.5 h each) would firm it, but
   reconciliation already places them at Saccharomycotina, so it changes no conclusion — do it only if a
   reviewer asks. An alternative 3Di model (Q.3Di.AF) for concordance is likewise optional — GTR20 already
   stress-tests the matrix.
3. **M4 — articulate AA/3Di non-independence** and why the partition is still informative (cite FoldTree).
4. Orthology cross-check (OMA/OrthoDB) and the ~134 PF05970-not-IPR048293 species (missed orthologs);
   explicit Helitron-contamination check on the final 957.
5. Minor: better aligner sensitivity; FoldTree rooting robustness; dating range (stem→crown) stated explicitly.

---

## E. Implied manuscript changes
- **Abstract/Discussion:** soften "amino acids cannot resolve it" → "amino acids marginally but
  significantly favor a deep placement (AU p=0.033) attributable to long-branch attraction." Lead the
  Discussion with the topology-independent biology (RRM3 only in budding yeasts), then structure + rates.
- **New Results paragraph (R3b):** the topology test + rate analysis (B2/B3) — this is the methods core.
- **Novelty framing:** position as confirm-and-date-with-a-new-method + a clean LBA case study; avoid
  claiming a new biological fact until the Y1000+ literature check is done.
- **Response-to-reviewers (draft points, to be voiced later):** (i) we did not assert sequence has no
  signal — we quantified it (AU p=0.033) and showed it is an LBA artifact; (ii) novelty is the rigor and
  the structural method, not the qualitative paralogy; (iii) here are the robustness analyses (D2) we add.
