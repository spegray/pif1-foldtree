# Two ways to crack a G-quadruplex: conservation of the G4-engaging motif in PIF1 and RecQ helicases across the tree of life

*Working manuscript draft. Methods and Results are written out in prose and figure-backed; Abstract,
Significance, Introduction, and Discussion are in Spencer's voice. Front/back matter (byline, affiliations,
funding, Zenodo DOI) is placeholder pending finalization. The RECQL4/5-like arm (Results R6) is complete but
exploratory — coarse predicted poses and no experimental co-structure in that clade — and flagged as such.*

**Spencer Gray**¹ *(co-authors TBD)*
¹ *Affiliation TBD*

---

## Significance
G-quadruplexes are knots of guanine-rich DNA that stall the replication fork, and cells keep specialized
helicases to unwind them. Two unrelated families do this job — PIF1 and RecQ — and each grips the
quadruplex with a single protruding residue: a basic "wedge" in PIF1, an aromatic "β-wing" in RecQ. We
asked how conserved those two motifs are across the whole tree of life, using ~2,000 AlphaFold3 models of
helicase–quadruplex complexes together with sequence and structure. The two families turn out to conserve
their solution in genuinely different ways: PIF1 holds a *chemistry* constant (a basic residue, everywhere),
while RecQ holds a *structural position* constant but lets the residue drift — and that drift is not random.
It tracks the RecQ paralogs, and it descends from an ancestral aromatic that specific eukaryotic lineages
independently abandoned. It is a concrete picture of two convergent DNA-repair machines conserving the same
function under two different rules.

## Working title
*One chemistry held constant, one position that drifted: structural conservation of the G4-engaging motif
in PIF1 and RecQ helicases across the tree of life.*

## Abstract
G-quadruplex (G4) DNA is a replication obstacle that dedicated helicases resolve, and two unrelated
families — PIF1 and RecQ — converge on a single-residue strategy for engaging it: PIF1 caps the 5′ G-tetrad
with a basic "wedge" (Arg324 in *Saccharomyces cerevisiae* Pif1), and RecQ caps the 3′ tetrad with an
aromatic "β-wing" residue projecting from the RQC winged-helix domain (Tyr564 in human RECQL1). How
conserved is each motif across the tree of life, in sequence and in three dimensions? We built a
standardized structural screen: 2,008 AlphaFold3 models of one helicase bound to a defined parallel G4 with
ATP and K⁺ (957 PIF1-family, 1,051 RecQ-family, sampled tree-of-life-wide), validated against the two
experimental complexes (PDB 8XAK and 9I22). AlphaFold3 reproduces G4 formation and the motif→tetrad contact
but places the quadruplex ~10 Å off its crystallographic pose, so we read a contact-level signal — whether
the motif engages the tetrad, and the identity of the residue that does — rather than atomic geometry. The
PIF1 wedge is a conserved basic residue family-wide (Arg or Lys in 84% of homologs), with a lineage-specific
Arg↔Lys swap; sequence and structure agree. The RecQ β-wing is not a conserved chemistry but a conserved
structural *position*: sequence alignment reports it as aromatic-enriched (74%), but a structure-based
homology mapping — Foldseek-superposing every RQC onto RECQL1 and reading the residue at the Y564 position —
returns 34%, the earlier figure inflated by averaging over paralogs. The β-wing residue is instead
paralog-diagnostic: RECQL1 and WRN keep an aromatic (Tyr, Phe), BLM converged on a conserved Asn, and the
fungal Sgs1 lineage uses small residues. Ancestral-state reconstruction on a RecQ tree recovers an aromatic
(Tyr) ancestor (maximum-likelihood P = 0.96), with the non-aromatic states arising as derived,
paralog-specific losses. The two families thus conserve the same function under two different rules — one
of chemistry, one of structural position — a distinction visible only once sequence, predicted structure,
and phylogeny are read together.

## Introduction
G-quadruplexes are four-stranded structures that guanine-rich DNA folds into wherever the sequence allows,
stacking planar G-tetrads around a central column of monovalent cations. They are thermodynamically stable
enough to block replication and transcription, and cells carry helicases whose job is to unfold them ahead
of the fork (Mendoza et al. 2016; Lerner and Sale 2019). Two structurally unrelated superfamily-2 helicase
families are central to that work in most organisms: the PIF1 family and the RecQ family. Both unwind G4s,
and, remarkably, both do it with a single amino-acid protrusion that pins the quadruplex.

In PIF1 the protrusion is an arginine in the "wedge" of the 1A subdomain. The *S. cerevisiae* Pif1–G4
co-crystal (Hu et al. 2024, PDB 8XAK) shows Arg324 wedged against the 5′-most G-tetrad, cation–π stacking on
the guanine face at the ssDNA–quadruplex junction where the helicase loads. In RecQ the protrusion is the
tip of the β-hairpin "wing" of the RQC (RecQ C-terminal) winged-helix domain. Structural and mutational work
on human RECQL1, BLM, and WRN identified this residue as the strand-separating element and, tellingly,
found it "structurally conserved but sequence-variable" — a tyrosine in RECQL1, an asparagine in BLM, a
phenylalanine in WRN (Pike et al. 2009; Kitano 2014; and the human RECQL1–G4 structures of [RECQL1-G4 2025],
PDB 9I22, where Tyr564 π-stacks the 3′-terminal tetrad). That single sentence in the literature — a
conserved position occupied by chemically different residues — is the question this paper takes to the whole
tree of life.

The question is worth asking from the substrate's side rather than the enzyme's. Both motifs are single
residues in short, divergent structural contexts (the PIF1 wedge in a fast-evolving helicase core; the RecQ
β-wing at the tip of a mobile loop), and asking "how is the residue that touches the G4 conserved?" is a
different question from "how is the protein conserved?". Sequence conservation alone cannot answer it
cleanly: the wedge sits in an alignable core, but the β-wing tip is exactly the kind of divergent loop that
multiple-sequence alignment cannot register to a single column (we show below that two standard aligners
disagree on it by ~11 columns). This is the regime where predicted structure earns its keep — a fold, and
the position of a residue within it, is conserved long after the sequence that specifies it has been
overwritten (Moi et al. 2025) — and where a co-modeled complex, even an imperfect one, can say whether a
homologous residue is positioned to touch the quadruplex at all.

Here we anchor on the substrate and screen both families structurally. We gather the PIF1- and
RecQ-family helicases tree-of-life-wide, model each on a defined parallel G4 with AlphaFold3, calibrate the
models against the two experimental complexes, and score whether — and with which residue — the G4-engaging
motif contacts its tetrad. Read alongside sequence conservation and, for RecQ, an ancestral-state
reconstruction, the two families give two different answers: PIF1 conserves a chemistry, RecQ conserves a
position.

---

## Methods

### Family sampling and homolog identification
- **PIF1 family (957 proteins).** Reused from the fungi-wide PIF1/RRM3 sample (Gray, companion manuscript):
  cellular PIF1-family helicases identified by **InterPro IPR048293** ("PIF1_RRM3_pfh1"), *not* Pfam
  PF05970, which also matches Helitron transposon helicases. Anchors: *S. cerevisiae* Pif1 (UniProt
  **P07271**, wedge **Arg324**), Rrm3 (**P38766**); *S. pombe* Pfh1 (**Q9UUA2**); human PIF1 (**Q9H611**).
- **RecQ family (37,658 candidates → 1,720 representatives).** RecQ-type helicases identified by **InterPro
  IPR004589** ("DNA helicase, ATP-dependent, RecQ type") across UniProt reference proteomes tree-of-life-wide
  (Eukaryota 14,026; Bacteria 23,500; Archaea 132 — archaea are RecQ-poor). Because RecQ is near-universal,
  raw counts track sequencing density rather than phylogenetic breadth; we subsampled to **one best
  representative per taxonomic family** (curated/annotated/AFDB-covered preferred), 1,720 proteins, with the
  experimental anchors forced in: human RECQL1 (**P46063**, β-wing **Tyr564**), BLM (**P54132**), WRN
  (**Q14191**), *S. cerevisiae* Sgs1 (**P35187**), *E. coli* RecQ (**P15043**).

### Core trimming and the RQC β-wing set
- PIF1 sequences trimmed to the PF05970 helicase core (HMMER `hmmsearch` envelope). The wedge (ScPif1 R324)
  falls within this core; its alignment column was located by mapping the anchor residue onto the untrimmed
  core MSA (MAFFT).
- RecQ sequences trimmed to the **RQC domain (Pfam PF09382**, the winged-helix module carrying the β-wing),
  at a relaxed threshold (E ≤ 1e-3) appropriate for the short domain. **1,051 of the 1,720** carry a
  detectable RQC and form the β-wing analysis set; **669 do not** — including human RECQL4 and RECQL5, whose
  RQCs are too divergent for the HMM — and are treated separately (R6). The RQC cores were aligned with MAFFT
  L-INS-i.

### Standardized AlphaFold3 complex
- Each job is one monomer helicase (full-length; the ~5,000-token server limit corresponds to ~4,900
  residues, above the largest helicase, so none is truncated) + a defined parallel-G4 DNA + one ATP + two
  K⁺, one model seed. **DNA:** the T7-AT11 construct of the ScPif1 co-crystal (Hu et al. 2024), a parallel
  G4 with an 8-nt ssDNA loading tail. Because PIF1 loads on a 5′ tail (5′→3′) and RecQ on a 3′ tail (3′→5′)
  and caps the opposite tetrad face, we keep the **same G4 core** but flip the tail: PIF1 jobs use `T7-AT11`
  (5′ tail), RecQ jobs use `AT11-T7` (3′ tail). A dual-tailed variant was tested and rejected — a dangling
  second tail degraded the ScPif1 model (wedge–tetrad distance 3.5 → 7.7 Å, ipTM 0.75 → 0.67). Job specs
  emitted programmatically (`workflow/19_make_af3_jobs.py`); non-standard residues (`X`, etc.) scrubbed to
  the nearest standard (X→G) before submission, flagged per protein.
- Two K⁺ ions match the two inter-tetrad cations of a three-tetrad parallel G4; ATP is the standard ligand.

### Validation gate
- AlphaFold3 models of the two experimental complexes were scored against them: ScPif1 + T7-AT11 vs **8XAK**
  (the co-crystal G4 core, differing by one tail thymine), and HsRECQL1 + c-Myc G4 vs **9I22**. Metrics:
  protein-core Cα RMSD after superposition, motif-residue displacement, motif→tetrad distance, G4 centroid
  offset, and nearest-neighbour guanine-O6 RMSD.

### Interface scoring
- For each top-ranked model (`model_0`), `workflow/20_g4_interface.py` / `22_score_all_g4.py`: G4 formation
  (a K⁺ coordinated by ≥ 4 guanine O6 within 3.5 Å), motif→engaged-tetrad distance (5′-most guanines for
  PIF1, 3′-most for RecQ), and interface confidence (ipTM, ranking score, protein–DNA chain-pair PAE) from
  the AlphaFold3 summary. The PIF1 wedge residue is the per-protein residue at the wedge alignment column;
  because the RecQ β-wing tip cannot be pinned by sequence, RecQ used the structure-based mapping below.

### RecQ β-wing by structural homology
- Every RecQ model's RQC domain was extracted and Foldseek-superposed onto the RECQL1 RQC; the residue
  aligning to RECQL1 **Tyr564** was read as that protein's β-wing residue (`workflow/23_recq_bwing_struct.py`;
  877/1,051 mapped, median TM to RECQL1 0.74). The method was validated by recovering the known anchor
  residues (RECQL1 Tyr564, BLM Asn1164, WRN Phe1037, *E. coli* RecQ His489). Paralog subfamily assignment
  used a Foldseek nearest-anchor classification against RECQL1, BLM, WRN, Sgs1, and *E. coli* RecQ.

### RecQ tree and ancestral-state reconstruction
- The RQC L-INS-i alignment (1,051 taxa) → **IQ-TREE** (LG+G4, fast search). The binary β-wing character
  (aromatic F/Y/W/H vs not; from the structural mapping) was reconstructed by **Fitch parsimony** under
  midpoint, bacterial-core, and single-outgroup rootings, and by a **maximum-likelihood symmetric two-state
  Mk model** (Felsenstein pruning with rate optimized by likelihood, rooted on the bacterial-core clade),
  reporting the root-state posterior probability. A distance-filtered variant (dropping β-wing calls > 15 Å
  from the tetrad, likely non-contacting) was run as a robustness check.

### Reproducibility
- All steps scripted (`workflow/18`–`23`) in a version-pinned conda environment; AlphaFold3 outputs, the
  interface table (`results/g4/interface_all.tsv`), the β-wing structural mapping
  (`results/g4/recq_bwing_struct.tsv`), and the validation records logged. Key tools: AlphaFold Server
  (AF3), Foldseek 10, MAFFT 7, HMMER 3.4, IQ-TREE 3.1.2, Biopython. Repo:
  github.com/spegray/pif1-foldtree.

---

## Results

### R1 — A tree-of-life AlphaFold3 screen, calibrated against the two experimental complexes
We modeled 2,008 helicase–G4 complexes — 957 PIF1-family and 1,051 RecQ-family proteins, each on a defined
parallel quadruplex with ATP and K⁺ — and first asked what such models can be trusted to report (Fig. 1).
Against the ScPif1–G4 co-crystal (8XAK), the top model reproduces the helicase fold (core Cα RMSD 1.4 Å over
501 residues) and places the wedge Arg324 within 2.2 Å of its crystallographic position; all five models
fold the quadruplex, with each K⁺ coordinated in the guanine-O6 channel, and the wedge contacts a 5′-tetrad
guanine at 3.5 Å, matching the 3.4 Å experimental contact. Against the human RECQL1–G4 structure (9I22) the
picture is the same in kind: the G4 folds, and the β-wing Tyr564 engages the 3′ tetrad (3.8 Å). But the
agreement is contact-level, not atomic: after superposing the proteins, the whole quadruplex sits ~10 Å off
its crystallographic pose (centroid offset 10.4 Å for ScPif1, 8.3 Å for RECQL1; nearest-neighbour O6 RMSD
~9.5 Å in both). AlphaFold3 gets the components and the fact of engagement right and the docking register
wrong. We therefore read a deliberately coarse signal throughout — *whether* the G4-engaging residue
contacts its tetrad, its identity, and the interface confidence — and make no claims about atomic binding
geometry. This is the honest ceiling of the method for these flexible G4 interfaces, and it is set by
experiment, not asserted.

### R2 — The PIF1 wedge is a conserved *basic* residue, tree-of-life-wide
Mapping the ScPif1 R324 wedge onto the 957-protein family alignment, the wedge column is occupied by a basic
residue — arginine or lysine — in **84%** of homologs (Arg 53%, Lys 31%), with the remainder scattered
across small and polar residues (Fig. 2). The conservation is of charge, not identity: the wedge is
Arg-dominant in Ascomycota and Lys-leaning in Basidiomycota, a lineage-specific swap between two residues
that present the same positive charge to the guanine face. The structural screen corroborates the sequence
signal on the confident subset: 94% of PIF1 models fold the G4, and where the wedge reaches the 5′ tetrad
its residue is basic in 92% of cases. (The absolute contact rate is modest — a quarter of wedges lie within
4 Å of the 5′ tetrad, more than half within 6 Å — as expected from the ~10-Å pose uncertainty of R1; the
robust readout is the *identity* of the contacting residue, not the tally.) The PIF1 solution to G4
engagement is a single conserved chemistry, held constant across the fungi and out to human.

### R3 — The RecQ β-wing is a conserved *position*, not a conserved residue — and sequence alignment overstates it
The RecQ β-wing behaves nothing like the PIF1 wedge. First, it cannot be pinned by sequence: two standard
alignments of the RQC domain (FFT-NS-2 and the more accurate L-INS-i) place the anchor residues in columns
~11 apart, and a column-based aromatic-conservation estimate swings from 74% to 61% between them — the
β-wing tip is exactly the divergent loop position that multiple-sequence alignment cannot register (Fig.
3A). We therefore let structure define homology: extracting every RQC domain and Foldseek-superposing it onto
RECQL1, we read the residue at the RECQL1 Tyr564 structural position directly (877/1,051 mapped; median TM to
RECQL1 0.74). The mapping is validated by recovering the known anchors without using their sequence — RECQL1
**Tyr564**, BLM **Asn1164**, WRN **Phe1037**, *E. coli* RecQ **His489** — precisely the "Tyr/Asn/Phe" spread
the literature describes.

Read this way, the β-wing tip is **aromatic in only 34%** of the family (F/Y/W 26%), Tyr-enriched but far
from conserved: two-thirds of RecQ homologs carry a non-aromatic residue there (Fig. 3B). The earlier
sequence figure (74%) was an alignment artifact of averaging a fuzzy column. So the RQC β-wing conserves a
*structural slot* — the winged-helix wing is present family-wide, and a residue occupies its tip — but the
chemistry of that residue is not conserved. This is a categorically different kind of conservation from the
PIF1 wedge, and it is why we anchored the analysis on structure rather than sequence.

### R4 — The β-wing residue is paralog-diagnostic
The 34% figure is a mixture, not a fact about any one lineage. Classifying each RQC-bearing RecQ to its
nearest subfamily anchor by structural similarity resolves the β-wing into five sharply different states
(Fig. 4; Table 1): the **RECQL1-like** proteins are 79% aromatic (Tyr-led), **WRN-like** 66% (Phe-dominant),
**bacterial RecQ** 51% (Tyr-led), the fungal **Sgs1-like** 16% (small residues — Ala, Gly, Ser), and the
metazoan **BLM-like** just 8% — but with a *conserved* asparagine at the β-wing position in 67% of them. The
β-wing is thus a paralog-diagnostic residue: RECQL1 and WRN hold an aromatic, BLM converged on a conserved
Asn as tightly as any aromatic signal, and the fungal lineage went to small residues. The apparent
"variability" of R3 was paralog pooling — the family-wide set is dominated by the 367 fungal Sgs1-like
proteins — and within a paralog the β-wing is far from random.

### R5 — The aromatic β-wing is ancestral, and the non-aromatic states are derived losses *(headline)*
The distribution across subfamilies has an evolutionary shape: the aromatic states sit in bacteria (the
ancestral, single-copy RecQ) and in the RECQL1/WRN paralogs, while the non-aromatic states (BLM's Asn,
Sgs1's small residues) sit in derived eukaryotic paralogs. Ancestral-state reconstruction on a RecQ RQC tree
makes this explicit. Under parsimony the β-wing reconstructs as **aromatic** at the root under both
robust rootings (midpoint and the bacterial-core clade), with the 414 non-aromatic eukaryotic tips falling
inside the BLM and Sgs1 clades as derived losses; only a single-tip outgroup rooting was ambiguous. A
maximum-likelihood two-state model puts a number on it: **root P(aromatic) = 0.96** across all 877 β-wing
calls, rising to **1.00** when non-contacting calls are removed (Fig. 5). The RecQ β-wing was ancestrally
aromatic — a tyrosine — and specific eukaryotic lineages independently abandoned it: BLM to asparagine, the
fungal Sgs1 line to small residues.

One caveat that arose here resolved in favour of this reading. Bacterial β-wings include a surprising
fraction of acidic residues (Glu/Asp, ~28% of bacterial calls), odd for a residue meant to contact a
negatively charged quadruplex. These turn out to be non-contacting: acidic-assigned bacterial β-wings sit
~28 Å from the tetrad (versus ~9 Å for aromatic ones; 2% versus 32% within 6 Å) despite comparable
domain-level alignment quality — the structural mapping landing on a surface residue in the most divergent
bacterial RQCs, not a real acidic engager. Excluding these non-contacting calls does not weaken the aromatic
ancestor; it sharpens it to P = 1.00.

### R6 — The RQC-less RecQ helicases (RECQL4/5-like) keep the loading face but lose the conserved element *(exploratory)*
A structural conservation story about the RecQ β-wing has an explicit boundary: 669 of the 1,720 sampled
RecQ-type helicases (39%) carry no detectable RQC/β-wing module at all. These are overwhelmingly eukaryotic
(540) and vertebrate-heavy, and they include human **RECQL4** and **RECQL5**, whose RQCs are too divergent to
detect — RECQL5-orthologs dominate the set (479) over RECQL4-orthologs (61). How, if at all, this large
minority engages a quadruplex without the canonical β-wing is a separate question, and one with no
experimental co-structure to calibrate against — nothing anchors it the way 8XAK and 9I22 anchor PIF1 and
canonical RecQ. We approached it as an open, unbiased screen. A balanced 82-protein batch (41 RECQL4-like +
41 RECQL5-like, human RECQL4/RECQL5 as anchors) was modelled on the same 3′-tailed substrate, and because
these proteins have no motif column to key off, we simply asked, for each top model, which residue's *side
chain* comes closest to the 3′ tetrad, how the protein sits, and whether anything about that engagement is
shared across the set (`workflow/24_recql4_g4.py`).

They do engage, and about as well as the canonical family. A quadruplex forms in 65% of models (versus 70%
for RQC-bearing RecQ), the folded models are the more confident ones (ipTM 0.71 vs 0.67; protein–DNA PAE 4.4
vs 5.9), and in 89% a protein side chain reaches the 3′ face — a median of 1.9 Å, 23 of 44 folded models
within 2 Å, genuinely close rather than a graze. That contact is made by the helicase core, not by the large
disordered N-terminal arms these proteins carry (only 6 of 73 engage through the N-terminal 15% of the
chain). So the ancestral RecQ instinct to load the 3′ end of the quadruplex survives the loss of the β-wing.

What does *not* survive is any conserved element making the contact. The engaging residue has no dominant
chemistry — 42% basic, 32% aromatic, the rest hydrophobic, polar or acidic — and the split is the same in
both sub-clades (RECQL4-like 45/35, RECQL5-like 39/27) and robust to filtering on fold and confidence; the
aromatic fraction, tellingly, does not rise in the confident subset but falls to 24%, so no conserved
aromatic is hiding under the noise. Nor is there a conserved *position*: measured against the one landmark
these divergent proteins share, the Walker-A P-loop, the engaging residue sits anywhere from 613 residues
upstream to 907 downstream (interquartile range 677 residues) — scattered across the whole architecture, not
parked in a slot the way the β-wing is. Even the two human anchors disagree with each other: RECQL4 engages
through R355, a basic residue; RECQL5 through F373, an aromatic. We read the 42% "basic" as generic
backbone-grip — long lysine and arginine side chains reaching the phosphate backbone at the junction, which
any minimum-distance metric will favour — not as a PIF1-style tetrad-face wedge, and we do not claim these
helicases have converged on one.

The reading that fits is a descending ladder of specificity. PIF1 conserves a
*chemistry* — a basic wedge at a fixed structural position. Canonical RecQ conserves a *position* — the
β-wing slot — while letting the residue drift. The RQC-less helicases conserve only a *polarity* — they still
find the 3′ face — while residue, chemistry and position all scatter. The β-wing, on this view, is the RecQ
family's dedicated quadruplex module: where it is present it is a conserved structural slot, and where it has
been lost, in the RECQL4/5 line, no conserved replacement has evolved to take its place (Fig. S6). This arm
is exploratory — coarse predicted poses, a third of models that never fold a clean G4, and no experimental
structure anywhere in the RQC-less clade to check it against — and we frame it as a hypothesis the ladder
predicts rather than a result the models prove. Full numbers in `results/g4/RECQL4_R6_SUMMARY.md`.

---

## Tables

**Table 1 — The RecQ β-wing residue is paralog-diagnostic** *(the R4 classification; Fig. 4).*

| Subfamily (nearest anchor) | n (RQC-bearing) | aromatic (F/Y/W/H) | dominant β-wing residue |
|---|---|---|---|
| RECQL1-like | 19 | 79% | Tyr |
| WRN-like | 62 | 66% | Phe |
| Bacterial RecQ | 350 | 51% | Tyr (+ non-contacting acidic) |
| Sgs1-like (fungal) | 367 | 16% | Ala / Gly / Ser |
| BLM-like (metazoan) | 79 | 8% | **Asn (67%)** |

---

## Figures

**Main:**
1. **Screen and validation** (`fig1_screen_validation`) — the pipeline (gather → standardized AF3 complex →
   validate → score) and the two calibration panels: ScPif1 vs 8XAK (core 1.4 Å, wedge 2.2 Å, contact 3.5 Å;
   G4 ~10 Å off) and HsRECQL1 vs 9I22 (β-wing 3.8 Å; G4 8.3 Å off). Establishes the contact-level ceiling.
2. **The PIF1 wedge** (`fig2_pif1_wedge`) — wedge-column residue frequencies across 957 homologs (Arg 53%,
   Lys 31%; 84% basic), with the Ascomycota-Arg / Basidiomycota-Lys clade split.
3. **Sequence vs structure for the RecQ β-wing** (`fig3_recq_seqvsstruct`) — (A) two aligners place the
   β-wing column ~11 positions apart (74% vs 61% aromatic); (B) the structure-based homology mapping onto
   RECQL1 Tyr564, anchor-validated (RECQL1 Y, BLM N, WRN F, EcRecQ H), returns 34% aromatic family-wide.
4. **Paralog-diagnostic β-wing** (`fig4_paralog`) — aromatic fraction by subfamily (RECQL1 79, WRN 66,
   bacterial 51, Sgs1 16, BLM 8%); BLM's conserved Asn. *(the two-kinds-of-conservation contrast with Fig. 2.)*
5. **An ancestral aromatic, eroded paralog by paralog** (`fig5_ancestral`) — RecQ RQC tree / schematic
   colored by β-wing chemistry; root reconstructs aromatic (ML P = 0.96); bacteria and RECQL1/WRN retain
   the ancestral Tyr, BLM (Asn) and fungal Sgs1 (small) are derived losses. *(headline)*

**Supplementary:**
- **S1** — the dual-tail substrate test (Pif1 wedge 3.5 → 7.7 Å, ipTM 0.75 → 0.67) justifying the tail-flip.
- **S2** — motif→tetrad distance distributions and G4-folding rates per family; the contact-rate/pose caveat.
- **S3** — the "closest-residue" metric artifact (returns RECQL1 Thr562, not the β-wing Tyr564) motivating
  the structural-homology mapping.
- **S4** — bacterial acidic β-wings are non-contacting (~28 Å vs ~9 Å); the acidic-caveat resolution.
- **S5** — ancestral-state robustness: parsimony rootings + ML posterior (0.96 all / 1.00 engaging-only).
- **S6** — the RQC-less clade (R6): (A) G4-fold and 3′-contact rates for RECQL4/5-like vs canonical RecQ;
  (B) engaging-residue chemistry (42% basic / 32% aromatic), unchanged across sub-clades and confidence
  filters; (C) engaging-residue position scattered across ±700 aa of the Walker-A P-loop (no conserved slot);
  (D) the descending-ladder schematic (chemistry → position → polarity-only) across PIF1 / RecQ / RECQL4-5.

---

## Discussion
Two unrelated helicase families crack the same knot of DNA with the same trick — a single residue jammed
against a G-tetrad — and the interesting result is that they conserve that trick under two different rules.
PIF1 conserves a *chemistry*: the wedge is basic (Arg or Lys) in 84% of homologs from fungi to human, and
the only thing that changes is which basic residue supplies the charge. RecQ conserves a *position*: the
β-wing is a structurally conserved slot at the tip of the RQC winged-helix, occupied by chemically different
residues that turn out to be diagnostic of the RecQ paralog — Tyr in RECQL1, Phe in WRN, Asn in BLM, small
residues in fungal Sgs1. Naming the difference this way required reading sequence, predicted structure, and
phylogeny together: sequence alone reports the β-wing as 74% aromatic, a number that dissolves once the
divergent loop is mapped by structure (34%) and then resolves into clean per-paralog states.

The evolutionary reading follows. Bacterial RecQ, the ancestral single-copy form, is aromatic (Tyr-led), and
the non-aromatic β-wings sit in derived eukaryotic paralogs; ancestral-state reconstruction recovers an
aromatic ancestor with maximum-likelihood posterior 0.96 (1.00 once non-contacting calls are removed). So
the "aromatic β-wing" of the textbooks is best read as the ancestral state, retained in RECQL1 and WRN and
independently abandoned in the BLM and Sgs1 lineages, rather than a universal feature of the family. That
BLM landed on a *conserved* asparagine, and fungal Sgs1 on small residues, suggests these are not neutral
drift but paralog-specific solutions worth their own functional study — the β-wing residue is one of the
few positions that cleanly separates the eukaryotic RecQ paralogs by structure.

The method's limits are set by experiment and stated plainly. AlphaFold3 reproduces G4 formation and the
fact of motif engagement but places the quadruplex ~10 Å off its crystallographic pose in both calibration
cases (8XAK, 9I22); the models are therefore a contact-level assay — is the motif positioned to touch the
tetrad, and which residue is it — not a source of atomic binding geometry. Every structural-conservation
statement here is an AlphaFold3-derived hypothesis in that sense, not an experimental contact, and the
family-wide *sequence*/identity signals (the PIF1 84%-basic, the RecQ paralog states, the ancestral
reconstruction) rest on the more robust structural-homology mapping rather than on per-model docking. Two
further caveats mark the edges. The RecQ ancestral reconstruction rests on a short-domain (RQC) tree that is
noisy at deep nodes, so the aromatic-ancestor result — though convergent across parsimony rootings, the ML
posterior, and the acidic-caveat resolution — would be firmed by a longer-domain (helicase + RQC) tree.
And the whole β-wing analysis excludes the ~40% of RecQ-type helicases that lack a detectable RQC, including
RECQL4 and RECQL5; our exploratory look at that clade (R6) finds they still dock the 3′ face but with no
conserved residue, chemistry, or position — consistent with the β-wing being the family's dedicated module
rather than one option among several — though coarse poses and the absence of any RQC-less co-structure keep
that a hypothesis rather than a settled result.

None of these caveats unsettle the central picture: two convergent DNA-repair machines conserve the same
function under two different rules, one holding a chemistry constant and one holding a structural position
constant while its residue drifts along paralog lines from an ancestral aromatic. It is a case where the
conservation of *function* and the conservation of *sequence* come apart in two distinct ways, visible only
because predicted structure lets us ask the question at the level of the single residue that does the work.

---

## References
*(Key references; to be completed and formatted on submission.)*

Hu X, et al. 2024. Structural basis for the recognition and unwinding of a G-quadruplex by the Pif1
helicase. *Nat Commun* 15:s41467-024-50575-8. https://doi.org/10.1038/s41467-024-50575-8 *(PDB 8XAK)*

Pike ACW, Shrestha B, Popuri V, Burgess-Brown N, Muzzolini L, Costantini S, et al. 2009. Structure of the
human RECQ1 helicase reveals a putative strand-separation pin. *Proc Natl Acad Sci USA* 106(4):1039–1044.
https://doi.org/10.1073/pnas.0806908106

Kitano K. 2014. Structural mechanisms of human RecQ helicases WRN and BLM. *Front Genet* 5:366.
https://doi.org/10.3389/fgene.2014.00366

[RECQL1-G4] et al. 2025. Structural mechanism of RECQ1 helicase in unfolding G-quadruplexes compared with
duplex DNA. *Nucleic Acids Res* 53(17):gkaf877. https://doi.org/10.1093/nar/gkaf877 *(PDB 9i1p/9i22/9i23)*

Bochman ML, Sabouri N, Zakian VA. 2010. Unwinding the functions of the Pif1 family helicases. *DNA Repair*
9(3):237–249. https://doi.org/10.1016/j.dnarep.2010.01.008

Mendoza O, Bourdoncle A, Boulé JB, Brosh RM Jr, Mergny JL. 2016. G-quadruplexes and helicases. *Nucleic
Acids Res* 44(5):1989–2006. https://doi.org/10.1093/nar/gkw079

Lerner LK, Sale JE. 2019. Replication of G quadruplex DNA. *Genes* 10(2):95.
https://doi.org/10.3390/genes10020095

Moi D, Bernard C, Steinegger M, Nevers Y, Langleib M, Dessimoz C. 2025. Structural phylogenetics unravels
the evolutionary diversification of communication systems in gram-positive bacteria and their viruses.
*Nat Struct Mol Biol* 32(12):2492–2502. https://doi.org/10.1038/s41594-025-01649-8

*(Software: AlphaFold3 / AlphaFold Server — Abramson et al. 2024, Nature; Foldseek — van Kempen et al. 2024;
MAFFT — Katoh & Standley 2013; IQ-TREE — Minh et al. 2020; HMMER — Eddy 2011; Biopython — Cock et al. 2009.)*

## Data and code availability
- Repo (scripts `18`–`23`, job specs, interface table, β-wing mapping, validation records):
  github.com/spegray/pif1-foldtree.
- To deposit on submission (Zenodo — TODO): the 2,008 AlphaFold3 job specifications and top models;
  `results/g4/interface_all.tsv`, `recq_bwing_struct.tsv`, `pif1_wedge_map.tsv`; the RecQ RQC alignment and
  tree; the validation superposition records; the two experimental references used (8XAK, 9I22, from the PDB).
- Sequences from UniProt; experimental structures from the PDB.

## Acknowledgments and funding
*(Placeholder — TBD.)*
