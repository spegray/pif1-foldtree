# PIF1 / RRM3 — dating the duplication in Ascomycota

**Question.** When did the two budding-yeast PIF1-family helicases — *PIF1* and *RRM3* —
split from a single ancestral gene? Concretely: on which branch of the fungal species tree
did the duplication producing the PIF1 and RRM3 lineages occur, and ~how long ago?

This repo is a documented, re-runnable pipeline answering that question.

---

## The single most important methods decision: IPR048293, not PF05970

The PIF1 helicase Pfam domain **PF05970** is carried by *two unrelated things*: (1) the cellular
PIF1/RRM3/Pfh1 helicases we care about, and (2) **Helitron** rolling-circle transposons, which
encode a PIF1-family helicase and are rampant in TE-rich fungal genomes. Searching by PF05970
therefore drowns the cellular signal in transposons (3–10× inflation; see table). We instead use
the InterPro entry **IPR048293 "PIF1_RRM3_pfh1"**, which is *specifically the cellular group* and
keeps every anchor (Pif1, Rrm3, Pfh1, human PIF1) while dropping the Helitrons.

| Taxon (reference proteomes) | PF05970 | **IPR048293 (used)** |
|---|---|---|
| Ascomycota (ingroup) | 2,255 | **681** (563 sp) |
| Basidiomycota | 2,393 | **235** (133 sp) |
| Mucoromycota | 742 | **36** (27 sp) |
| Chytridiomycota | 33 | **4** (4 sp) |
| *S. cerevisiae* | 2 | **2** (Pif1+Rrm3 ✓) |
| *S. pombe* | — | **1** (Pfh1 ✓) |
| Human | 2 | **1** (PIF1, Q9H611 ✓) |

---

## Verified facts (checked 2026-06-15)

| Item | Value |
|---|---|
| Seed: *S. cerevisiae* **Pif1** | UniProt **P07271** (PIF1_YEAST, 859 aa) |
| Seed: *S. cerevisiae* **Rrm3** | UniProt **P38766** (RRM3_YEAST, 723 aa) |
| Outgroup anchors | *S. pombe* **Q9UUA2** (Pfh1) · human **Q9H611** (PIF1) |
| **Family filter (primary)** | InterPro **IPR048293 "PIF1_RRM3_pfh1"** (cellular PIF1/RRM3/Pfh1) |
| Looser models (not primary) | Pfam PF05970 (Helitron-contaminated), PF21530; InterPro IPR010285, IPR051055 |
| Cellular PIF1 set | **957 proteins / ~728 species** (Asco 681 ingroup + Basidio 235 + Mucoro 36 + Chytrid 4 + human 1) |
| AFDB structures | prediction API `https://alphafold.ebi.ac.uk/api/prediction/<acc>`; current model version **v6** |

### Early copy-number signal (≥2 vs 1 cellular paralog per species, by subphylum)

| Subphylum / class | 1-copy sp | ≥2-copy sp |
|---|---:|---:|
| **Saccharomycotina** (budding yeasts) | 14 | **91** |
| Pezizomycotina (filamentous asco) | **430** | 17 |
| Taphrinomycotina (incl. fission yeast) | **11** | 0 |
| Agaricomycotina (mushrooms) | 34 | 58 |
| Ustilaginomycotina | 19 | 3 |
| Pucciniomycotina | 14 | 3 |
| Mucoromycotina / Glomeromycotina | 21 | 5 |

**Read:** budding yeasts are overwhelmingly two-copy while their Ascomycota sisters
(filamentous + fission yeast) are overwhelmingly single-copy → the PIF1/RRM3 duplication
looks **restricted to Saccharomycotina** (stem or early crown). The two-copy mushrooms are
almost certainly an **independent** duplication — *copy-counting cannot tell shared ancestry
from convergence; the gene tree + reconciliation must.* (Hence comprehensive Basidiomycota
sampling.) This is a hypothesis the pipeline tests, not a conclusion.

---

## Directory layout

```
pif1-foldtree/
  README.md            # this file
  environment.yml      # conda env (Stages 3-6); data scripts need no env
  manifest.csv         # ONE row per protein: id, taxonomy, structure provenance, paralog label
  workflow/
    01_gather_homologs.py        # UniProt -> cellular PIF1 (--interpro IPR048293) per taxon  
    02_combine_cellular.py       # combine per-taxon cellular tables -> selected.tsv          
    02_select_representatives.py # (older) stratified subsample; unused now the set is small
    03_fetch_structures.py       # AFDB prediction API -> structures + manifest              
    04_fetch_fasta.py            # UniProt -> selected.faa (sequences for alignment)        
  data/
    seqs/
      cellular/    per-taxon IPR048293 tables (ascomycota.tsv, basidiomycota.tsv, ...)
      selected.tsv all 957 cellular proteins, tagged group/role  (THE master list)
      selected.faa all 957 sequences (from 04)
      _superseded_PF05970/  earlier Pfam-based candidate files (kept for the audit trail)
    structures/  afdb/ (.cif), af3/ (gap predictions), cores/ (corecut)
    species_tree/  published fungal tree, pruned to our taxa
  results/
    seq_tree/ struct_tree/ reconciliation/ dating/ figures/
```

---

## Runbook

### Environment  (Apple Silicon → Rosetta osx-64; see `environment.yml`)
```bash
CONDA_SUBDIR=osx-64 conda env create -f environment.yml   # invoke conda's bin directly in scripts
conda activate pif1
conda config --env --set subdir osx-64
conda list > env.lock.txt      # record exact versions (audit trail)
```

### Stage 1 — gather cellular PIF1 per taxon
```bash
for spec in 4890:ascomycota 5204:basidiomycota 1913637:mucoromycota 451435:chytridiomycota 9606:human; do
  python3 workflow/01_gather_homologs.py --taxon ${spec%%:*} --interpro IPR048293 \
      --out data/seqs/cellular/${spec##*:}.tsv
done
```

### Stage 0/1 — combine into the master list
```bash
python3 workflow/02_combine_cellular.py --indir data/seqs/cellular --out data/seqs/selected.tsv
# selected.tsv is THE auditable taxon/protein list; group/role columns flag ingroup vs outgroup.
```

### Stage 1b — fetch sequences 
```bash
python3 workflow/04_fetch_fasta.py --in data/seqs/selected.tsv --out data/seqs/selected.faa
```

### Stage 2 — fetch structures + manifest
```bash
python3 workflow/03_fetch_structures.py --in data/seqs/selected.tsv \
    --manifest manifest.csv --outdir data/structures/afdb --fmt cif
# proteins with no AFDB model are tagged MISSING_predict_with_AF3 in manifest.csv -> AF3/ColabFold.
```

### Stage 3 — corecut to the helicase core
PIF1 proteins have disordered N/C termini (full-length AFDB mean pLDDT ~60–68; the folded
helicase core scores higher). Trim BOTH sequences and structures to the shared core so
architecture differences don't confound either tree (the paper's "corecut", Fig. 4).
- Sequences: `hmmsearch` the PIF1 HMM (PF05970 model is fine for *locating* the domain envelope) → take the aligned core per protein.
- Structures: slice each `.cif` to the matching residues (Biopython) → `data/structures/cores/`.
- Or run the **fold_tree** corecut step (it does this automatically).

### Stage 4a — sequence ML gene tree
```bash
mafft --maxiterate 1000 --localpair data/seqs/cores.faa > results/seq_tree/aln.fasta
trimal -in results/seq_tree/aln.fasta -out results/seq_tree/aln.trim.fasta -automated1
iqtree  -s results/seq_tree/aln.trim.fasta -m MFP -B 1000 -alrt 1000 -pre results/seq_tree/pif1
# env ships IQ-TREE v3.1.2; the binary is `iqtree` (or `iqtree3`), NOT `iqtree2`
```

### Stage 4b — FoldTree structural tree
FoldTree Colab (zero install) — upload `data/structures/cores/`:
- Colab: https://colab.research.google.com/github/DessimozLab/fold_tree/blob/main/notebooks/FoldTree.ipynb
- Repo:  https://github.com/DessimozLab/fold_tree
Foldseek all-vs-all → Fident (3Di+AA, statistically corrected) → QuickTree → MAD root.
Best-of-both: concatenate AA+3Di, IQ-TREE partitioned (helps the deep nodes).

### Stage 5 — species tree, pruned to our taxa
Prune a published genome-scale fungal tree (Li et al. 2021 *Curr Biol*; Y1000+/Shen et al.
2018 *Cell* for the budding-yeast portion) to our taxids (dendropy / R `ape`).
Fallback: BUSCO `fungi_odb10` single-copy orthologs + IQ-TREE.

### Stage 6 — reconciliation
```bash
generax --families families.txt --species-tree data/species_tree/species.nwk \
        --rec-model UndatedDL --prefix results/reconciliation/generax
```
Also run **NOTUNG** (Java GUI). Cross-check GeneRax vs NOTUNG, and sequence-tree vs structure-tree.

### Stage 7 — absolute dating 
1. **Bracket.** Map the duplication branch onto a published *time-calibrated* fungal tree
   and read its age window (stem ↔ crown age of the clade it maps to). TimeTree / Shen 2018 / Li 2021.
2. **Our own relaxed clock.** Date the Stage-4a gene tree (fixed topology) with
   **treePL** (penalized likelihood, fast) and/or **MCMCtree** (PAML, Bayesian, approx. likelihood),
   placing calibrations on speciation nodes (e.g. Asco/Basidio split; Saccharomycotina crown; the
   *Saccharomyces* whole-genome-duplication landmark ~100 Mya); the duplication-node age falls out
   with a credible interval.
Report the node, support values, sequence-vs-structure agreement, the paralog retention/loss pattern,
and both date estimates.

---

## QC checklist
1. Cross-check IPR048293 membership against an independent orthology source (OMA/OrthoDB) — and
   check the ~134 Ascomycota species that have a PF05970 protein but *no* IPR048293 call, in case
   a real cellular ortholog is merely under-annotated.
2. UFBoot/SH-aLRT support OK on the duplication node and its two bracketing nodes.
3. Sequence- and structure-derived gene trees place the duplication on the **same** branch.
4. GeneRax and NOTUNG agree on that branch.
5. Paralog loss/retention pattern is biologically coherent (fission-yeast single Pfh1, and the
   single-copy filamentous fungi, resolved as pre-duplication vs loss *by the tree*, not assumed).
6. Answer is stable to outgroup choice and to MAD vs outgroup rooting.
7. The two-copy mushrooms (Agaricomycotina) resolve as an *independent* duplication, not as
   orthologs of yeast Pif1/Rrm3.
8. Both dating approaches (bracket + relaxed clock) give overlapping ranges.

## Platform notes
- Apple Silicon (arm64): use the Rosetta (osx-64) conda env (full bioconda coverage; Rosetta 2 present).
  In scripts call conda's binary directly (`$(conda info --base)/bin/conda`) — the shell *function*
  isn't available in non-interactive shells (caused an early `__conda_exe: permission denied`).
- NOTUNG is a Java JAR (download from the Notung site); needs a JRE (`conda install -c conda-forge openjdk`).
- Run FoldTree via Colab to avoid installing the structural toolchain locally.

