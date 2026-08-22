# Internal review — "Dating the PIF1/RRM3 duplication in budding yeasts"

*Adversarial, multi-axis review of `docs/MANUSCRIPT_SKELETON.md`, 2026-07-13. Nine independent
reviewers (two factual fact-checkers cross-checking every number against the `results/` files, two
citation auditors, claim-support, internal-consistency, a Spencer-voice editor, and two reviewer-2
lenses — methods and interpretation), each finding then re-checked by a skeptic against the same
ground-truth file before it survived. 40 raw findings; 32 verified. Line numbers refer to the
manuscript as reviewed.*

---

## Bottom line

The manuscript is fundamentally sound and close to submission. The headline result — 3Di structure
pulls the Pif1/Rrm3 duplication off a spurious base-of-Fungi placement onto the Saccharomycotina
stem, dated ~326–383 Mya — is backed by a genuinely strong robustness architecture, and every
load-bearing fact that could be checked against the repo verifies (AU *p* = 0.033, ΔlnL = 143; the
LG+C20/PMSF base-of-Fungi placement at ~916/917 tips; the copy-number tallies; the anchor
accessions). **There are no fabricated or wrong numbers of consequence.**

The dominant problem is not factual error but **overstatement of inferential strength**: two central
claims assert more than the analyses deliver. (1) GeneRax was run in fixed-topology EVAL mode, yet is
described as an "independent" method that "confirms" the placement and "integrates over gene-tree
error." (2) The long-branch-attraction diagnosis is *inferred* from rate asymmetry plus a
model-robustness result, not *demonstrated* by the outgroup-removal or site-stripping tests a
phylogenetics referee will demand. A second, cosmetic-but-pervasive class is bookkeeping drift
(three unreconciled Saccharomycotina counts, a 720-vs-719 backbone, a stale figure label, a seed
misattribution). Two voice items rise to major because a banned construction carries the thesis in
two of the three front-matter blocks.

After the interpretive claims are softened to match the methods, one or two LBA diagnostics are added
(or the thesis is reworded), the counts are reconciled, and the voice/citation fixes land, this is
submission-ready.

---

## Must-fix (a referee will raise these)

1. **Demonstrate LBA, don't just infer it.** *(major; Discussion 313–317, R2b 176–188)* The thesis
   needs the significant AA preference for the deep tree (AU *p* = 0.033) to be a *demonstrated*
   artifact, but the repo supports only rate asymmetry (Rrm3 1.22× faster; Saccharomycotina ~2×)
   plus PMSF — necessary, not sufficient. The PMSF argument is double-edged: a profile-mixture model
   is the standard LBA *remedy*, so its failure to move the node is equally consistent with real deep
   signal. **Decision:** either run one direct diagnostic — delete the human + early-diverging-fungal
   outgroups and re-infer the AA tree (if Saccharomycotina monophyly emerges from sequence alone,
   that is the clean proof), and/or slow-fast site-stripping — or reword the thesis to claim only
   that structure *overrides* sequence, without asserting the sequence signal is proven wrong.

2. **Stop calling the reconciliation independent / confirmatory / error-integrating.** *(major;
   Significance 18, Abstract 38, R4 200–205)* GeneRax ran `--strategy EVAL`, SPR radius 0 (Methods
   143–144) — on the *fixed* AA+3Di topology. It cannot independently confirm a topology it is
   handed, and it does not marginalize over the UFBoot uncertainty (node_448 / the Rrm3 daughter) it
   is invoked to overcome. Demonstrably input-driven: reconciling the AA-only gene tree puts node_448
   at only D = 8 (mushroom node_578 dominates), vs D = 43 for AA+3Di. **Fix:** reword R4 to
   "Because full gene-tree search was intractable, we fixed the ML topology (EVAL) and tested
   robustness instead by re-running reconciliation on three alternative gene trees and two
   species-tree backbones; all place the event on the Saccharomycotina ancestor." Drop "integrates
   over that gene-tree error and recovers the same event regardless" (204) and "independent and
   error-tolerant method" (201); soften "confirms" → "is consistent with / maps to" in the
   Significance and Abstract. To *genuinely* integrate over gene-tree error, reconcile the 1000
   UFBoot trees and report the fraction placing the duplication on the Saccharomycotina node.

3. **Remove the banned negative-parallelism from the front matter.** *(major; Significance 15,
   Abstract 31)* "the amino-acid signal at that depth is not merely weak but actively misleading"
   and "the amino-acid signal is not just weak but misleading" are the "not X but Y" correlative frame
   on the explicit banned list, and it carries the thesis in two of the three blocks a reader hits
   first. Recast to the asyndeton/colon form used in the Discussion itself (313: "the amino-acid core
   had not merely run quiet, it pointed the wrong way"). E.g. Significance: "…the amino-acid signal
   at that depth is worse than weak: it points the wrong way." Abstract: "…at the depth where the
   answer lies the amino-acid signal does not merely go quiet; it actively misleads."

4. **Reconcile the 103 / 104 / 105 Saccharomycotina counts.** *(minor but pervasive;
   R4 200, R3 192, Table 1 261)* Each is correct for a different object, but "all 104 sampled
   Saccharomycotina" (200) directly contradicts Table 1's 105 sampled. State the objects: 105 species
   sampled (manifest/Table 1); 104 in the species tree = node_448 clade (one dropped —
   *Wickerhamomyces ciferrii*, taxid 1206466, absent from both backbones); 103 species contribute a
   gene to the 197-tip AA+3Di clade. Add a one-clause reconciliation the first time they could confuse.

5. **Report the structure-only FoldTree placement** (currently the cleanest rebuttal you're not
   using). *(claim-support / adversarial-methods; R3, Methods 107–112, Fig 1)* The Methods promise
   "the gene tree three ways" but Results report only AA and AA+3Di. The fident FoldTree already
   resolves the Pif1/Rrm3 MRCA to 203 tips = 197 Saccharomycotina + only 6 others — i.e. **structure
   alone recovers the clade**, which directly answers the "the 3Di partition just pads and
   down-weights the significant AA signal" objection. Add one sentence (noting whether
   fident/alntmscore/lddt agree), or explicitly demote the structure-only tree to a supplement.

---

## By axis

### 1. Factual accuracy — clean; four minor slips (all verified against `results/`)
- **Seed misattribution** *(Methods 105)*: the main AA ML tree (LG+I+G4) and the LG+C20/PMSF check
  used `-seed 12345`, not 42 (`pif1.log`, `pif1_pmsf.log`); only `aa_bl`, both AU tests, and the
  AA+3Di tree used 42. Either correct line 105 to `-seed 12345` or soften the Reproducibility line
  (151) to "seeds fixed and recorded in each `.log`."
- **"median core pLDDT 88 in both predictors"** *(99–100; Fig S2 caption 304)*: AFDB core median is
  88.3 (n = 828) but ColabFold is 86.5 (n = 129, rounds to 87) — and ColabFold sits just *below* the
  pLDDT ≥ 87 robustness threshold used in R3. Write "88 (AFDB) and 87 (ColabFold)" or "~87–88 in
  both predictors."
- **Abstract "728 fungal species"** *(32)*: one of the 728 is human PIF1 (the outgroup named in the
  same parenthetical). Write "728 species (727 fungal plus human PIF1)." Methods (80) and Intro (66)
  are already correct.
- **node_578 = Agaricomyc*etes*, not Agaricomyc*otina*** *(Fig 4 caption 294; Fig 2 caption 286)*:
  the node subtends 61 leaves, all class Agaricomycetes; the 30 sampled Tremellomycetes of the
  subphylum are not under it. The Discussion (343) and all internal results docs already say
  "Agaricomycetes." D = 23 is correct.
- **Verified correct** (do not touch): AU *p* = 0.033 / ΔlnL = 143; PMSF base-of-Fungi at ~916/917
  tips; copy-number tallies (91/105, 430/447, 0/11, 58/92); anchor accessions.

### 2. Citations — all real and correctly attributed; coverage gaps to close
- **Attribution:** every one of the ~32 references is a real paper with correct authors/year/journal,
  DOIs resolve (spot-verified via PubMed). No fabrications, no mis-credits — 3DiPhy (Puente-Lelièvre
  2023), the FoldTree correction and structural-phylogenetics claim (Moi 2025), the AU test
  (Shimodaira 2002), and the IPR048293/PF05970 credits all check out.
- **Orphans (listed, never cited):** PhyML/Guindon 2010 (drop — IQ-TREE does the ML work); Varadi
  2024 (the in-text AFDB cite says only "2022" beside "model version v6," a post-2024 release — make
  it "Varadi et al. 2022, 2024"); six biology refs cited nowhere in the body — Crandall 2024, David
  2025, Harman & Manna 2016, Hérivaux 2018, Malone 2022, Opulente 2024 (cite in the Intro where they
  fit, or drop); thirdkind mentioned at Fig S3 without "(Penel et al. 2022)".
- **Missing primary citations:**
  - **MAJOR — Helitron / PF05970** *(82–84)*: the *sole* justification for the family filter (the
    "key decision") carries no cite. Add Thomas & Pritham 2015 (*Microbiol Spectr*) and/or
    Kapitonov & Jurka 2001 (*PNAS*) — the Helitron RepHel protein carries a PIF1-family helicase
    domain, which is exactly why PF05970 hits it.
  - **MAJOR — WGD "~100 Mya"** *(248, 298, 330)*: a dated claim central to the R6 comparison, uncited.
    Add Wolfe & Shields 1997 (*Nature*) at minimum; optionally Kellis et al. 2004 and/or
    Marcet-Houben & Gabaldón 2015.
  - **MINOR** — Tajima 1993 (relative-rate test); Le & Gascuel 2008 (LG matrix); Quang et al. 2008
    (C20) + Wang, Minh, Susko & Roger 2018 (PMSF — note the correct author list, not "von Haeseler");
    Huerta-Cepas et al. 2007 (the species-overlap *algorithm*, distinct from the ete3 software cite);
    Eddy 2011 (HMMER, which defines the core envelope).
  - **NIT / optional** — GTR20 has no canonical paper (just say "estimated from the data"); LBA
    (Felsenstein 1978, optional); recPhyloXML (Duchemin 2018, optional); UniProt / NCBI Taxonomy
    (very low).
  - **Adequately cited, do not flag:** the Intro Pif1/Rrm3 function statements are fully covered by
    Boulé & Zakian 2006 and Bochman et al. 2010.

### 3. Claim support
- The GeneRax overstatement (must-fix 2) is the main one.
- The structure-only FoldTree result is computed but unreported (must-fix 5).
- **Dating uncertainty caveat missing** *(R6 242–247; Discussion 326–327)*: 326/383/404 Mya are
  single point ages read straight off the Shen 2018 chronogram (no interval propagation in
  `workflow/17`), and the tree's per-node 95% HPD intervals are never mentioned. Add a clause that
  the ages are point estimates inheriting Shen's node-age CIs (not re-propagated here) and that
  326–383 is a crown-vs-stem *bracket*, not a statistical CI.

### 4. Internal consistency
- 103/104/105 species (must-fix 4).
- **720 vs 719 backbone** *(Methods 131)*: the NCBI tree has 720 leaves, but GeneRax pruned one
  further geneless leaf (taxid 1041607) → 719, matching the grafted backbone; both reconcile 941
  genes across 719 species. State the pruning so the two numbers don't read as an unexplained
  discrepancy.
- **Stale Fig 5 label** *(`fig5_timetree.R` line 54)*: reads "~330–383 Mya"; the crown age is 326.
  Change to "~326–383 Mya" and regenerate `fig5_timetree.svg` / `_outlined.pdf`.
- Abstract "728 fungal" (also a factual item).

### 5. Voice (against the spencer-voice guide)
- **Negative parallelism** "not X but Y" ×2 (must-fix 3).
- **Softer correlative** *(Discussion 322)*: "they are not an independent second dataset so much as
  the fold-level constraint…" — recast positive ("…carry fold-level constraint that outlasts the
  sequence once it has saturated; they are a second view of the same sites, not a second dataset").
  The Methods (120) already state the caveat as a plain negation, so two negations can collapse to one.
- **Em-dash parenthetical aside** *(99)*: "— median core pLDDT 88 in both predictors —" → parentheses
  (the only such construction in the manuscript).
- **"gold-standard" GeneRax** *(139)*: drop the epithet; if a ranking is wanted, say why plainly
  ("which models duplication and loss explicitly, unlike the species-overlap tally").

### 6. Adversarial — methods (reviewer-2, phylogenetics)
- **F1 — reconciliation non-independence / circularity** (feeds must-fix 2). Concrete path to keep
  the "integrates over error" claim honestly: reconcile the 1000 UFBoot trees and report the fraction
  on the Saccharomycotina node.
- **F2 — LBA demonstrated, not inferred** (must-fix 1).
- **F3 — root choice** *(111; already flagged as a to-do in `REVIEW_RESPONSE.md`)*: the FoldTree is
  rooted on the single, longest human PIF1 branch — itself the invoked LBA attractor — with no
  root-choice sensitivity test. Basidiomycota/Mucoromycota/Chytridiomycota outgroups are already in
  the sample; root on them or present an unrooted AU analysis, and note the reconciliation is robust
  to gene-tree rooting.
- **F5 — the structure-only FoldTree rebuts the "3Di dilutes AA" objection** (must-fix 5).

### 7. Adversarial — interpretation (reviewer-2, molecular evolution)
- **Rate/identity confound in the paralog-polarity reading** *(R5 233; Discussion 337)* — the sharpest
  interpretation critique. The identity gap (54.5% to Pif1 vs 52.3% to Rrm3; 442/536 outgroups closer
  to Pif1) and the longer Rrm3 branch are **not independent**: a faster-evolving Rrm3 (R2b, 1.22×)
  mechanically has lower % identity to *every* outgroup, so identity-to-outgroup is the *expected
  consequence of rate alone* — one rate signal double-counted, not separate evidence of
  ancestral-character retention. **Also:** `WINDOWS_RESULTS.md` line 79 ("core rates near-symmetric")
  contradicts R2b's 1.22×; reconcile that internal-doc statement before it reaches a reviewer. State
  the confound explicitly, drop the double-counting, and either defend polarity with a rate-corrected
  / ancestral-state argument or hedge "ancestral retention" further.
- **Gene conversion never mentioned** *(0 grep hits)* — a first-order alternative that would in
  principle corrupt both the duplication-node topology and the 442/536 asymmetry. But the protected
  quantities are largely *insulated* (the 442/536 outgroups are single-copy with no paralog to
  convert with; conversion would *shrink* the identity gap; the reciprocally monophyletic paralog
  clades are the opposite of a conversion signature), so a one-sentence name-and-dismiss likely
  suffices; a positive check (are within-species Pif1–Rrm3 pairs anomalously similar?) would be stronger.
- **HGT not testable under the DL model** *(140)*: UndatedDL forbids transfer by construction —
  every count shows T = 0 because the model can't infer it, not because it was tested. Add one
  sentence, and justify why HGT is implausible (nuclear, vertically inherited helicase).

---

## Strengths — preserve these
- **The robustness architecture is the real evidence.** Three AA+3Di reruns (AFDB-only D = 37,
  pLDDT ≥ 87 D = 39, GTR20 D = 43), two species-tree backbones (NCBI D = 9, grafted D = 43), and the
  predictor/quality/matrix controls. The fixes above ask you to *lean on these* rather than on the
  overstated "independent reconciliation" claim.
- **The family filter is a genuine methodological insight**, well documented (Pfam PF05970 →
  InterPro IPR048293, with the count deflation quantified in Table S1 / Fig S1). It just needs its
  transposon citation.
- **The hedging is honest and consistently placed** — UFBoot 29 disclosed up front, the asymmetry
  held "as an asymmetry worth naming rather than a settled fact," the inflated grafted-tree totals
  given "one honest qualification." Keep this candor.
- **Predictor-consistent structure handling** (ColabFold to match the predictor; identical residue
  coordinates for AA and 3Di; core-pLDDT gating) pre-empts obvious referee objections.
- **The Discussion already contains the in-voice version of the thesis** the front matter botches
  (313) — so the voice fix is copy-your-own-solution, not invention.

---

## Open questions — author decisions only Spencer can make
1. **LBA:** run outgroup-removal / slow-fast site-stripping to *demonstrate* it, or reword the thesis
   to claim only that structure overrides sequence? (New analysis vs. paragraph edit.)
2. **GeneRax:** reconcile the 1000 UFBoot trees so you can legitimately keep "integrates over
   gene-tree error," or drop that claim and rest on the reruns + backbones?
3. **Structure-only FoldTree:** promote its 197-Saccharomycotina placement to a standalone R3 result,
   or demote it to a supplementary input/robustness check? (The text currently promises three trees
   and delivers two.)
4. **Paralog polarity:** defend "Pif1 retains ancestral character" with a rate-corrected / ASR
   argument, or hedge/drop it given the rate confound?
5. **Gene conversion:** one-sentence name-and-dismiss, or run an actual homogenization check?
6. **Seed provenance:** correct line 105 to `12345`, or generalize the Reproducibility statement (and
   note that seed 42 covers only the downstream runs)?
7. **Internal doc:** reconcile `WINDOWS_RESULTS.md` "near-symmetric rates" vs the manuscript's 1.22×.

---

*Process note: the citation-coverage reviewer crashed mid-run on a transient server error and was
re-run separately to complete the missing-primary-citation sweep; all other axes completed in the
main pass. Findings marked "verified" were reproduced by a second agent against the same ground-truth
file. Full per-agent transcripts under
`.claude/.../subagents/workflows/wf_45e92671-cc9/`.*
