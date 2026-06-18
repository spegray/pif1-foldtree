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

This round addressed **M1** and **M2** (the two that decide whether this is a methods paper, a
confirmatory note, or needs reframing). M3–M6 and the minors remain open (Section D).

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

---

## D. Open items NOT addressed this round (prioritized)
1. **M5 — substitute the published phylogenomic species tree.** *Largely addressed* (`workflow/17`):
   mapping the duplication's species set onto the Shen et al. 2018 time-calibrated Y1000+ tree
   (`data/species_tree/shen2018_timetree.newick`) places it on a node subtending ~96% of budding
   yeasts (317/330 tips) dated **~330–383 Mya** — so the Saccharomycotina placement is *robust to the
   phylogenomic backbone*, not an artifact of NCBI's arbitrary polytomy resolution, and the date is now
   read directly off a calibrated tree. *Remaining:* a full GeneRax DL reconciliation against a grafted
   broad-fungi + Shen tree (for the per-branch loss pattern) — needs the Windows box (GeneRax/Rosetta).
2. **M3/M6 robustness:** bootstrap stability of the AA+3Di topology; an alternative 3Di model
   (Q.3Di.AF) for concordance; a pLDDT-filtered and predictor-balanced (AFDB vs ColabFold) rerun.
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
