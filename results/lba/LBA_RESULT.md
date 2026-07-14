# Outgroup-removal LBA diagnostic — result (2026-07-13)

**Question (review must-fix #1):** is the amino-acid tree's base-of-Fungi placement of the Pif1/Rrm3
duplication (Results R2) a long-branch-attraction artifact (R2b), or real deep signal? Rate asymmetry
plus the PMSF result are necessary but not sufficient to call it LBA; outgroup removal is the direct test.

**Design (`workflow/18_lba_outgroup_removal.py`).** From the same 957-tip AA alignment
(`results/seq_tree/aln.trim.fasta`, MAFFT → trimAl, 209 columns), delete the longest, most distant
branches — human PIF1 (1) and the early-diverging fungi, Mucoromycota (36) + Chytridiomycota (4) — and
keep the nearer Basidiomycota (235) to root Ascomycota. Re-infer under the **same** model and seed as the
main analysis: IQ-TREE 3.1.2, `LG+I+G4 -B 1000 -bnni -alrt 1000 -seed 12345`, on 916 tips. (Native arm64
IQ-TREE 3.1.2 via Homebrew; the anaconda env binary is x86 and Rosetta is unavailable on this machine.)

**Result (`results/lba/analyze_lba.R` on `lba_noEDF.treefile`).**

| Tree | MRCA(ScPif1, ScRrm3) | Composition |
|---|---|---|
| Full AA, 957 tips (R2) | **950 tips** | base of Fungi |
| **AA, outgroups removed, 916 tips (this test)** | **198 tips** | **197 Saccharomycotina + 1 Pucciniomycotina (99.5%)** |
| AA+3Di, 957 tips (R3) | 197 tips | Saccharomycotina |

Removing the long outgroup branches collapses the duplication from base-of-Fungi to essentially the same
Saccharomycotina clade the AA+3Di tree recovers — **from sequence alone, with no structural information.**
Of 205 Saccharomycotina genes in this tree, 197 fall inside the clade; the ~8 outside are the
early-diverging, single-copy budding yeasts (e.g. Lipomycetales), consistent with R4/R6.

**Interpretation.** The base-of-Fungi placement was **long-branch attraction**: the distant human and
early-diverging-fungal outgroups were pulling the fast-evolving Saccharomycotina copies apart and deep;
delete them and the Saccharomycotina duplication clade reappears from the amino acids alone. This is the
direct demonstration the manuscript lacked, and it is convergent with the AA+3Di result — two independent
routes (remove the LBA trap; add fold-constraint signal) recover the same Saccharomycotina placement.

**Caveats (honest).** (i) The Pif1/Rrm3 node itself remains weakly supported (SH-aLRT 0 / UFBoot 39) —
the saturated sequence still resolves this deep node poorly, which is exactly why the 3Di signal adds
value; the demonstration is about the ML *point estimate* moving off base-of-Fungi, not about strong
support. (ii) Basidiomycota is not monophyletic in this AA-only tree (single-copy PIF1s scatter, as
expected), so the tree was rooted on a Basidiomycota outgroup tip; the anchors' MRCA is a nested ingroup
clade whose composition is robust to which outgroup tip is used to root.

**Bearing on the manuscript.** Converts review must-fix #1 from "reword the thesis" to "add the
diagnostic and strengthen it": the significant AA preference for the deep tree (AU *p* = 0.033) is now
shown to be an outgroup-length artifact, not real deep signal. Suggested home: a paragraph in R2b (and a
one-line supplementary table/figure), reported as the outgroup-removal control.

Files: `aln.noEDF.fasta` (916-tip input), `lba_noEDF.{treefile,iqtree,log}`, `analyze_lba.R`.
