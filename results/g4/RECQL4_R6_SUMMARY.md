# R6 — G4 engagement in the RQC-less RecQ helicases (RECQL4/5-like) (2026-07-12)

Scored the full 82-job AlphaFold3 batch (41 RECQL4-like + 41 RECQL5-like, human RECQL4/RECQL5 anchors) on the
uniform 3′-tailed substrate (`GGTGGTGGTTGTTGTGGTGGTGGTGGT`+`T9`), top model each, with
`workflow/24_recql4_g4.py` → `results/g4/recql4_interface.tsv`. These helicases **lack the RQC/β-wing**, so
there is no known motif column to key off (unlike scripts 22/23); the scan is UNBIASED — for each model it
finds the protein residue whose **side chain** (not backbone) comes closest to the 3′ G-tetrad (the loading
face), records its identity/chemistry, the full ≤4 Å contact footprint, and where it sits in the chain.
Read against the Phase-0 caveat (AF3 reproduces G4 formation + the contact, but the pose is coarse ~8–10 Å).

## Engagement is real and about as good as canonical RecQ
- **G4 folded in 53/82 (65%)** — cf. canonical RecQ 70%, PIF1 94%. Folded models are the more confident ones
  (median ipTM 0.71 vs 0.67 unfolded; protein–DNA PAE 4.42 vs 5.85) — the fold calls track confidence, not noise.
- **89% (73/82) make a side-chain contact ≤4 Å at the 3′ face** (83% among folded). Median side-chain→3′-tetrad
  distance **1.90 Å** (≤2 Å in 23/44 folded) — genuine close engagement, not a marginal graze.
- Engagement is via the **helicase core, not disordered N-terminal arms**: only 6/73 engage through the
  N-terminal 15% of the chain.

## …but there is NO conserved engaging element (three independent axes)
1. **Identity / chemistry — heterogeneous.** The engaging residue is **42% basic / 32% aromatic / 12%
   hydrophobic / 8% polar / 5% acidic** (R 21, K 10, Y 9, F 7, H 4, S 3, W 3, E 3). No dominant type
   (contrast PIF1's 84% basic wedge). Same mix in **both** sub-clades — RECQL4-like 45% basic / 35% aromatic,
   RECQL5-like 39% / 27%. Robust to filtering: folded-only 39/30, confident (ipTM≥0.70, PAE≤6) 48/24 — and the
   aromatic fraction does **not** rise with confidence (falls to 24%), so no conserved aromatic is hiding in
   the noise.
2. **Structural position — scattered (confound-free).** Relative to the conserved Walker-A P-loop (`TGxGK[ST]`,
   located in 42/44 folded+contact models), the engaging residue sits anywhere from **613 aa N-terminal to
   907 aa C-terminal** of it (IQR 677 aa; 11/42 upstream, 19/42 within 0–500 aa downstream). No conserved
   slot — contrast canonical RecQ's single β-wing position. (Raw sequence position `pos_frac` spans 0.02–0.99,
   IQR 0.37–0.92; the Walker-A offset removes the N-terminal-length confound and shows the same scatter.)
3. **The human anchors themselves diverge.** RECQL4 engages via **R355 (basic)**, RECQL5 via **F373
   (aromatic)** — the two best-characterized RQC-less helicases do not even agree with each other.

## Verdict — completes the cross-family ladder
RQC-less RecQ helicases still dock the **correct (3′) face** of the G4 — the ancestral RecQ loading polarity
survives — but they engage it with **no conserved residue, chemistry, or structural position**: contact
reverts to generic, position-variable, basic-leaning core–DNA gripping. Losing the β-wing costs the family its
dedicated G4-engaging module, and **no conserved replacement evolves**. This closes a descending ladder of
G4-recognition specificity across the screen:

| Family | What is conserved | Element |
|---|---|---|
| **PIF1** | a **chemistry** (84% basic) | rigid 1A-domain Arg-wedge, fixed position |
| **RecQ (RQC⁺)** | a **structural position** (β-wing slot) | plastic residue; ancestral aromatic eroded paralog-specifically |
| **RecQ (RQC⁻ / RECQL4/5-like)** | only the **loading polarity** (3′ face) | no conserved residue/position — diffuse core contact |

## Honest caveats (R6 is exploratory)
- Coarse AF3 poses (Phase-0 gate); no experimental co-structure exists for **any** RQC-less RecQ–G4 complex to
  calibrate against — nothing anchors this arm the way 8XAK/9I22 anchor PIF1/RecQ.
- The min-distance metric favors long basic side chains reaching the DNA backbone, so the 42% "basic" most
  likely reflects **generic backbone-grip, not a PIF1-like tetrad-face wedge** — we do **not** claim these
  helicases use a basic wedge.
- 65% fold rate: a third of models never form a clean G4, and are excluded from the interface statistics.
- A formal Foldseek structural-superposition test (à la script 23) was not run: with coarse poses on
  full-length multidomain proteins it would be under-powered, and the chemistry + Walker-A-offset scatter
  already establish non-conservation.
