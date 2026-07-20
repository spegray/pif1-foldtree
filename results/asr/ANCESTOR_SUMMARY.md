# Pre-duplication Saccharomycotina PIF1 — ancestral reconstruction (2026-07-13)

**Question:** can we reconstruct the single ancestral PIF1 that existed on the Saccharomycotina stem just
before it duplicated into Pif1 and Rrm3, and is the Arg-wedge already present in it?

**Method.** Marginal (empirical-Bayes) ancestral sequence reconstruction with IQ-TREE 3.1.2 (`-asr`), on
the structure-informed AA+3Di topology held fixed (`-te pif1_aa3di.treefile`), under the same LG+I+G4 model
as the main sequence analysis (`workflow/19_asr_ancestor.py`; extraction `results/asr/extract_ancestor.R`).
The reconstructed node is the **duplication node = MRCA of ScPif1 (P07271) and ScRrm3 (P38766)** — IQ-TREE
`Node14`, subtending exactly the **197-tip, 197-species Saccharomycotina clade** (the same clade as R3/R4).
Scope: the trimmed **helicase core (209 columns)** only; the accessory N/C-terminal domains were trimmed
family-wide and are not reconstructed here.

## Result — the ancestor is well-resolved, and the wedge is ancestral

- **Confidence is high** (a shallow ancestor with 197 well-sampled descendants reconstructs well): mean
  per-site posterior **0.94**; **89%** of sites ≥ 0.8, **98%** ≥ 0.5. The earlier saturation caveat bites
  only at a minority of fast loop sites.
- **The Arg-wedge is retained.** ScPif1 R324 maps to alignment column 92 (position 89 of the fold-ready
  sequence); the reconstructed ancestor carries **R (Arg) there with posterior 0.961.** So the basic wedge
  was already present in the pre-duplication Saccharomycotina ancestor — it predates the Pif1/Rrm3 split
  rather than being a Pif1-specific innovation. This dovetails with the G4 side-project's finding of a
  conserved basic wedge across the family.

**Fold-ready ancestral core (206 aa; 3 high-gap columns dropped)** — `results/asr/ancestor_dupnode.fasta`:

```
TLSEEQQHVLDMVVQGKSIFFTGSAGTGKSVLLREIIKRLRKKYGPDSVAVTASTGLAACNIGGTTLHSFAGIGLGNESVEQLVKKIRRNKKSRQRWRNTKVLIIDEISMIDGELFDKLDQIARKIRKNDKPFGGIQLVITGDFFQLPPVSKDNNQPAKFCFESESWKECIKHTIVLTQVFRQKDNEFIDMLNEMRLGKLSPETEQ
```

The Walker-A P-loop (`GSAGTGKS`, ~res 23) and Walker-B (`IIDEISM`, ~res 100) reconstruct cleanly; the wedge
sits in a basic patch (`...VKKIRRNKKSR...`, the R at fold-position 89). Per-site posteriors and the kept/dropped
column flags are in `results/asr/ancestor_confidence.tsv`.

## Fold it on the AlphaFold3 server (± G4)

`results/asr/af3/ancestor_fold_jobs.json` is upload-ready for alphafoldserver.com — two jobs, using the
**same PIF1 substrate as the 956 modern models** for direct comparability:

1. **`ancestor_dupnode_G4`** — ancestor + parallel G4 (`TTTTTTTT`+`GGTGGTGGTTGTTGTGGTGGTGGTGGT`, 5′ loading
   tail) + 1 ATP (CCD_ATP) + 2 K⁺. *Does the ancestral wedge reach the 5′ G-tetrad?*
2. **`ancestor_dupnode_apo`** — ancestor alone. *Is the wedge structurally poised without substrate?*

What to look for: after superposing on a modern PIF1–G4 model, is the position-89 Arg side chain at/near the
5′ tetrad face (the modern moderns median ~5.5 Å, ≤ 4 Å in the tightest quartile)? Read the aggregate, not
one model's exact geometry (AF3's Phase-0 caveat: it reproduces G4 formation and the motif→tetrad *contact*
but places the quadruplex ~10 Å off its crystallographic pose).

**Caveats.** (i) The modern jobs folded full-length proteins; this ancestor is **core-only**, so the G4 may
fold less reliably without the rest of the protein — compare the core/wedge region, not global pose.
(ii) For a pose-robust read on wedge engagement, bump `modelSeeds` to `[1,2,3,4,5]` before submitting.

## Fold result — the wedge engages (AlphaFold3; `results/asr/score_ancestor_fold.py`)

Scored with the same metric as the 956 modern models (`workflow/22`):

- **+G4:** the quadruplex forms in **5/5 models**; the ancestral wedge (R89) reaches the 5′ G-tetrad at
  **median 5.57 Å (best 4.04 Å; per-model 5.57/6.82/4.52/4.04/10.02), ipTM 0.72–0.78** —
  indistinguishable from the modern PIF1 aggregate (median 5.5 Å, G4 folded 94%). Scores in
  `results/asr/ancestor_fold_scores.tsv`.
- **Apo (−G4):** the wedge folds as an ordered Arg in all five models (pLDDT 96, pTM 0.88).

So the reconstructed pre-duplication ancestor was already a G4-engaging PIF1 with the wedge in place; the
capability predates the Pif1/Rrm3 split. Written up as Results **R7** + **Table S3**.

## For the manuscript

- A clean, self-contained result: the ancestral core reconstructs at high confidence and **already carries
  the basic wedge**, which firms up the "conserved basic wedge" story with a reconstructed ancestor rather
  than an inference from extant tips.
- Partial bearing on the R5 polarity question (open item): the ancestor lets us ask, at functional sites,
  whether Pif1 or Rrm3 retains the ancestral state — a *character-level* test that is cleaner than the raw
  %-identity comparison (though the global identity-to-ancestor comparison still inherits Rrm3's faster
  rate, so it does not fully escape that confound). Per-site ancestral-state retention (Pif1 vs Rrm3) is the
  natural next analysis.
