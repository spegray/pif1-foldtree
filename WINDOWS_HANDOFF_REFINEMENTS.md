# Windows / WSL2 handoff #2 — robustness refinements

**Fresh Claude Code session on the Windows desktop (WSL2, native x86-64).** You have no memory of the
prior work; everything you need is in this repo. A Mac session did the analysis and pre-built all
inputs; your job is the *compute* (IQ-TREE + GeneRax run natively here, fast and crash-free) and to
push results back. Ask before anything destructive; otherwise proceed.

## What this is
We dated a gene duplication: the budding-yeast PIF1/RRM3 helicase pair arose in the **Saccharomycotina
ancestor (~330–383 Mya)**. The key result is that the protein's **3Di structure** resolves a deep node
that amino-acid sequence gets wrong (long-branch attraction). An adversarial review flagged robustness
items; this run closes them. **For each rerun the test is the same:** rebuild the AA+3Di tree on a
perturbed input, reconcile it, and confirm the **Pif1/Rrm3 duplication still maps to Saccharomycotina**
(robust) rather than "Fungi" (broke). All inputs are pre-built and committed.

## Setup
```bash
cd ~/pif1-foldtree && git pull          # get the pre-built inputs + scripts
conda activate pif1recon                 # reuse the env from handoff #1 (IQ-TREE 3.1.2, ete3, dendropy, biopython)
#   if missing: conda create -y -n pif1recon -c conda-forge -c bioconda iqtree generax dendropy ete3 biopython pandas numpy
```
Use `iqtree2` or `iqtree3` — whichever your `pif1recon` env provides (`which iqtree2 iqtree3`). Below I write `iqtree2`.
Run everything **from the repo root** (partition files use repo-relative paths). ete3 downloads the NCBI
taxdump on first `12_reconcile.py` run (a few minutes, once).

---

## R1 — Predictor-batch test (review M6)  ·  *primary*
Does the 828-AFDB vs 129-ColabFold split drive the result? Re-run on **AFDB structures only** (828 tips):
```bash
iqtree2 -s results/seq_tree/aln_AFDBonly.fasta -p results/seq_tree/aa3di_partition.nex \
        -t results/seq_tree/start_AFDBonly.treefile -n 50 -T AUTO -seed 42 \
        -pre results/seq_tree/rerun_AFDBonly -redo
python workflow/12_reconcile.py --tree results/seq_tree/rerun_AFDBonly.treefile
```
**Pass = the anchor MRCA(Pif1,Rrm3) line still reads "→ Saccharomycotina".**

## R2 — pLDDT-stress test (review M6)
The helicase cores are already uniformly high-confidence (median core pLDDT **88**; 5th–95th pct 84–91),
so this worry is largely moot — but to be thorough, re-run on the **top ~¾ by core pLDDT** (≥87, 750 tips):
```bash
iqtree2 -s results/seq_tree/aln_pLDDT.fasta -p results/seq_tree/aa3di_partition.nex \
        -t results/seq_tree/start_pLDDT.treefile -n 50 -T AUTO -seed 42 \
        -pre results/seq_tree/rerun_pLDDT -redo
python workflow/12_reconcile.py --tree results/seq_tree/rerun_pLDDT.treefile
```

## R3 — Matrix-robustness test (3Di model)
Is the result specific to the 3DiPhy matrix? Re-run the **full** alignment with the 3Di partition under
**GTR20+FO+G** (IQ-TREE estimates the 3Di exchange rates from our own data instead of using the fixed matrix):
```bash
iqtree2 -s results/seq_tree/alnAA_3di.fasta -p results/seq_tree/aa3di_partition_gtr20.nex \
        -t results/seq_tree/pif1_aa3di.treefile -n 50 -T AUTO -seed 42 \
        -pre results/seq_tree/rerun_gtr20 -redo
python workflow/12_reconcile.py --tree results/seq_tree/rerun_gtr20.treefile
```
(Orthogonal support already exists: the foldseek FoldTree — a totally different structural method —
also recovered the Saccharomycotina clade. R3 just confirms it isn't a 3DiPhy-matrix artifact.)

## R4 — GeneRax loss pattern on a phylogenomic tree (review M5)  ·  *optional / advanced*
The placement is already shown robust to the phylogenomic backbone (Mac `workflow/17` on Shen 2018) and
the NCBI-tree GeneRax DL run (handoff #1, `results/reconciliation/aa3di/`) already gives a loss pattern.
The only thing left is the loss pattern on a *phylogenomic* species tree, which needs a grafted tree:
1. Graft the Shen 2018 budding-yeast topology (`data/species_tree/shen2018_timetree.newick`) into the
   non-Saccharomycotina backbone of `data/species_tree/ncbi_species.nwk`, matching tips by genus+species
   (dendropy; ~54–205 of our Saccharomycotina species match Shen — keep matched, leave the rest at their
   NCBI position). Write `data/species_tree/grafted_species.nwk`.
2. `python workflow/14_prep_generax.py --genetree results/seq_tree/pif1_aa3di.treefile --label aa3di_phylo`
   then point its families.txt species tree at the grafted tree and run GeneRax 2.0.4
   (`--strategy SPR -r UndatedDL`, env `gx204` — NOT 2.1.3, which segfaults). Compare the loss pattern.

Skip R4 if short on time — it's a refinement of an already-robust result.

## Report back
Write `results/refinements/WINDOWS_REFINEMENTS_RESULTS.md` with, for R1–R3 (and R4 if run):
- the anchor MRCA(Pif1,Rrm3) mapping (Saccharomycotina? or did it move?),
- duplication/speciation event counts from `12_reconcile.py`,
- one line: does the Saccharomycotina placement hold? (expected: yes for all three).
Then push:
```bash
git add results/seq_tree/rerun_*.treefile results/seq_tree/rerun_*.iqtree results/refinements WINDOWS_HANDOFF_REFINEMENTS.md
git commit -m "Robustness refinements (Windows): predictor, pLDDT, 3Di-matrix reruns"
git push
```
Tell the user: *"Refinements done — push complete; go back to the Mac and pull."*

---
**One-line mission:** rebuild the AA+3Di tree three ways (AFDB-only, high-pLDDT, GTR20 3Di model),
confirm each still places the PIF1/RRM3 duplication in **Saccharomycotina**, and push the results.
