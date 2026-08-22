# G4-engaging motif conservation — structural screen results (2026-07-08)

Scored 2,007 AF3 top models (956 PIF1 + 1,051 RecQ) with `workflow/22_score_all_g4.py`.
Table: `results/g4/interface_all.tsv`. Read against the Phase-0 caveat: AF3 reproduces G4 formation and the
motif→tetrad CONTACT, but the G4 pose is coarse (~8–10 Å off), so tight per-model distances are noisy — trust
the *aggregate* + the *identity of the contacting residue*, not any single model's geometry.

## PIF1 — the Arg-wedge (structure corroborates sequence) ✅
- **G4 folded in 94%** of models.
- Wedge→5′-tetrad distance: median 5.5 Å; **≤4 Å in 25%, ≤6 Å in 57%, ≤8 Å in 83%** — the wedge sits at/near
  the 5′ face in a majority; only tight stacking is undercounted (AF3 pose coarseness). Anchor ScPif1 R324 → 1.98 Å.
- **Among ≤4 Å contacts, 92% of wedges are basic (R/K)** (R 184, K 41). By clade: ascomycota 30% contact / 97%
  basic; basidiomycota 14% / 66%; mucoro 17% / 67%.
- **Verdict:** the structural screen corroborates the sequence signal (84% basic wedge). The PIF1 G4-engaging
  residue is a conserved BASIC wedge, tree-of-life-wide, strongest in Ascomycota — clean structure+sequence agreement.

## RecQ — the RQC β-wing (AF3 too coarse for a clean family-wide call) ⚠️
- **G4 folded in only 70%** (vs PIF1 94%) — AF3 folds the RecQ complex less reliably.
- The initial "closest RQC residue → 3′ tetrad" metric was **backbone-biased**: it returned mostly basic (K/R
  39%) and loop residues, NOT the stacking aromatic (anchor check: it picked RECQL1 **T562**, not the known
  **Y564**). So the raw "19% aromatic" was a metric artifact, not biology. RecQ's RQC grips the DNA backbone
  with basic residues; the β-wing aromatic STACK is a separate, subtler contact.
- Aromatic-specific metric (closest F/Y/W/H in the RQC → 3′ tetrad) correctly recovers **RECQL1 Y564 → 3.82 Å**,
  but family-wide aromatic engagement is weak/noisy: **≤5 Å in 33%** (euk 36%, bac 30%), median 7.7 Å.
  Anchors are mixed — RECQL1 Y564 3.82 Å (clean); BLM nearest-aromatic Y1160 7.3 Å (BLM's tip is the known
  non-aromatic Asn — expected); WRN F1037 11.5 Å and Sgs1 F1192 14.4 Å (G4 mis-docked / not folded in those models).
- **Verdict:** AF3's RecQ–G4 poses are too coarse to confirm the aromatic β-wing at family scale. The aromatic
  β-wing conservation rests on the ANCHOR structures (experimental 9I22 Y564; the Phase-0 validation) + the
  sequence signal (aromatic-enriched, alignment-fuzzy) + RQC fold conservation — the family-wide per-model
  structural statistic adds little beyond the best cases.

## RecQ β-wing by STRUCTURAL homology (the rigorous test) — reframes the hypothesis
`workflow/23_recq_bwing_struct.py`: Foldseek-superposed every RQC onto RECQL1 and read the residue at the
Y564-equivalent structural position (`results/g4/recq_bwing_struct.tsv`). Mapped 877/1051 (TM-to-RECQL1
median 0.74, all >=0.4). METHOD VALIDATED by the anchors — recovers RECQL1 **Y564**, BLM **N1164** (the known
non-aromatic Asn), WRN **F1037** (Phe), EcRecQ **H489** — exactly the literature's Tyr/Asn/Phe spread,
regardless of sequence divergence.
RESULT: the β-wing tip is **structurally-conserved in POSITION but chemically VARIABLE**. Identities:
Y 18%, then A/F/E/G/H/N/M... a long tail. **Aromatic (F/Y/W/H) = 34%** (F/Y/W = 26%) — Tyr-enriched (~3x over
baseline) but NOT a conserved chemistry; two-thirds of the family carry a non-aromatic residue there. By clade,
aromatic is MORE bacterial (bac 52%, euk 22%). The earlier sequence estimate (~74% FFT-NS-2 / 61% L-INS-i) was
inflated by alignment; the rigorous structural number is **34%**. So "conserved aromatic β-wing" was an
over-generalization from RECQL1/WRN — the family reality is a conserved structural SLOT with a plastic residue.

## Cross-family headline
Two unrelated helicase families engage a G-quadruplex with a single-residue protrusion of DIFFERENT chemistry —
PIF1 a **basic** wedge (5′ tetrad), RecQ an **aromatic** β-wing (3′ tetrad, opposite face). PIF1's rigid
single-residue wedge is captured cleanly by both sequence and AF3 structure; RecQ's flexible β-wing loop on a
harder-to-fold G4 is resolved cleanly only in the best models (RECQL1). The structure-over-sequence thesis lands
firmly for PIF1; for RecQ it stays anchor- and sequence-supported with limited AF3 corroboration.

## To strengthen RecQ (options, not yet done)
1. Per-protein STRUCTURAL homology: Foldseek-superpose each RQC onto RECQL1 and read the Y564-equivalent residue
   directly (rigorous; the pilot approach at scale — ~1,051 superpositions).
2. Restrict to the confident/well-folded subset (high ipTM + G4-folded) before aggregating.
3. Detect aromatic ring STACKING geometry (coplanarity), not just min distance.
