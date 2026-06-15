# PIF1 / RRM3 — dating the duplication in Ascomycota

**Question.** When did the two budding-yeast PIF1-family helicases — *PIF1* and *RRM3* —
split from a single ancestral gene? Concretely: on which branch of the fungal species tree
did the duplication producing the PIF1 and RRM3 lineages occur, and ~how long ago?

This repo is a documented, re-runnable pipeline answering that question. It is the running
lab notebook: every command, version, and decision is recorded here and in `manifest.csv`.

---

## The approach in one paragraph

A gene-duplication date is a **reconciliation** result, not a "look at the tree" result.
We build one **gene tree** of all *cellular* PIF1-family proteins across an even, fungi-wide
taxon sample (Ascomycota ingroup + Basidiomycota / early-diverging-fungi / human outgroups),
then **reconcile** it against a trusted fungal **species tree** (GeneRax / NOTUNG):
reconciliation maps the duplication node onto a specific species-tree branch — *that branch is
the answer*. A **sequence-based maximum-likelihood tree is the workhorse** (the PIF1/RRM3 split
is recent and high-identity, the regime where sequence ML matches or beats structure — see the
FoldTree paper, Moi et al. *Nat Struct Mol Biol* 2025, doi:10.1038/s41594-025-01649-8, Fig. 2d
& Discussion). A **structural (FoldTree) tree is the confirmatory layer**, valuable for the deep
PIF1-family backbone and rooting. Structures come from the **AlphaFold Database** (no mass AF3).

See the full plan: `~/.claude/plans/users-spencergray-desktop-s41594-025-01-lexical-snowglobe.md`.

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
  README.md            # this file — the lab notebook
  environment.yml      # conda env (Stages 3-6); data scripts need no env
  manifest.csv         # ONE row per protein: id, taxonomy, structure provenance, paralog label
  workflow/
    01_gather_homologs.py        # UniProt -> cellular PIF1 (--interpro IPR048293) per taxon  [DONE]
    02_combine_cellular.py       # combine per-taxon cellular tables -> selected.tsv          [DONE]
    02_select_representatives.py # (older) stratified subsample; unused now the set is small
    03_fetch_structures.py       # AFDB prediction API -> structures + manifest               [running]
    04_fetch_fasta.py            # UniProt -> selected.faa (sequences for alignment)          [DONE]
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

## Status

- [x] **Filter decided** — InterPro **IPR048293** (cellular PIF1), replacing PF05970 (Helitron-contaminated).
- [x] **Stage 1 — gather** (`01 --interpro IPR048293`) → `data/seqs/cellular/*.tsv`.
- [x] **Outgroups added** — comprehensive: all Basidiomycota + Mucoromycota + Chytridiomycota cellular PIF1 + human PIF1 (Q9H611).
- [x] **Stage 0/1 — combine** (`02_combine_cellular`): **957 proteins** → `data/seqs/selected.tsv` (tagged ingroup/outgroup; keep-all).
- [x] **Stage 1b — sequences** (`04`): 957 → `data/seqs/selected.faa`.
- [x] **Stage 2 — structures** (`03`): **828/957 from AFDB** (median pLDDT 64.6); **129 lack an AFDB model** (extracted to `data/structures/af3/to_predict_129.faa`).
- [ ] **Stage 2b — fold the 129 gaps**: **localcolabfold on the home RTX 4090 via WSL2** (AF2, to match AFDB; lean `--num-models 1 --num-recycle 3`, no Amber; ~3–6 h). 115 ≤1500 aa fold fine; 14 >1500 aa (two ~2.7 kb *Candidozyma*) may OOM at 24 GB → `--max-msa 512:1024` or core-only. Steps in `data/structures/af3/HOWTO_predict_129.md`; integrate with `workflow/06_integrate_predictions.py` → manifest.
- [x] **Conda env** (`pif1`, osx-64/Rosetta): built; `env.lock.txt` written. IQ-TREE binary is `iqtree` (v3.1.2), not `iqtree2`.
- [x] **Stage 3 — corecut** (`05`): all 957 trimmed to the PF05970 helicase core → `data/seqs/cores.faa` + `data/seqs/tip_map.tsv` (tree-safe labels + taxids).
- [~] **Stage 4a — sequence ML tree**: MAFFT FFT-NS-i + trimAl → 209-col alignment (`results/seq_tree/aln.trim.fasta`); **IQ-TREE running** (`-m MFP -B 1000 -alrt 1000`, 957 taxa).
- [ ] **Stage 4b** — FoldTree structural tree (Colab).
- [ ] **Stage 5** — prune a published fungal species tree to our taxa.
- [ ] **Stage 6** — reconcile → the duplication branch. **(the answer)**
- [ ] **Stage 7** — absolute dating: bracket from a published clock **+** our own relaxed clock.

---

## Runbook

### Environment  (Apple Silicon → Rosetta osx-64; see `environment.yml`)
```bash
CONDA_SUBDIR=osx-64 conda env create -f environment.yml   # invoke conda's bin directly in scripts
conda activate pif1
conda config --env --set subdir osx-64
conda list > env.lock.txt      # record exact versions (audit trail)
```

### Stage 1 — gather cellular PIF1 per taxon  *(done)*
```bash
for spec in 4890:ascomycota 5204:basidiomycota 1913637:mucoromycota 451435:chytridiomycota 9606:human; do
  python3 workflow/01_gather_homologs.py --taxon ${spec%%:*} --interpro IPR048293 \
      --out data/seqs/cellular/${spec##*:}.tsv
done
```

### Stage 0/1 — combine into the master list  *(done)*
```bash
python3 workflow/02_combine_cellular.py --indir data/seqs/cellular --out data/seqs/selected.tsv
# selected.tsv is THE auditable taxon/protein list; group/role columns flag ingroup vs outgroup.
```

### Stage 1b — fetch sequences  *(done)*
```bash
python3 workflow/04_fetch_fasta.py --in data/seqs/selected.tsv --out data/seqs/selected.faa
```

### Stage 2 — fetch structures + manifest  *(running)*
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

### Stage 4a — sequence ML gene tree (workhorse)
```bash
mafft --maxiterate 1000 --localpair data/seqs/cores.faa > results/seq_tree/aln.fasta
trimal -in results/seq_tree/aln.fasta -out results/seq_tree/aln.trim.fasta -automated1
iqtree  -s results/seq_tree/aln.trim.fasta -m MFP -B 1000 -alrt 1000 -pre results/seq_tree/pif1
# env ships IQ-TREE v3.1.2; the binary is `iqtree` (or `iqtree3`), NOT `iqtree2`
```

### Stage 4b — FoldTree structural tree (confirmatory)
FoldTree Colab (zero install) — upload `data/structures/cores/`:
- Colab: https://colab.research.google.com/github/DessimozLab/fold_tree/blob/main/notebooks/FoldTree.ipynb
- Repo:  https://github.com/DessimozLab/fold_tree
Foldseek all-vs-all → Fident (3Di+AA, statistically corrected) → QuickTree → MAD root.
Best-of-both: concatenate AA+3Di, IQ-TREE partitioned (helps the deep nodes).

### Stage 5 — species tree, pruned to our taxa
Prune a published genome-scale fungal tree (Li et al. 2021 *Curr Biol*; Y1000+/Shen et al.
2018 *Cell* for the budding-yeast portion) to our taxids (dendropy / R `ape`).
Fallback: BUSCO `fungi_odb10` single-copy orthologs + IQ-TREE.

### Stage 6 — reconciliation = the answer
```bash
generax --families families.txt --species-tree data/species_tree/species.nwk \
        --rec-model UndatedDL --prefix results/reconciliation/generax
```
Also run **NOTUNG** (Java GUI). Cross-check GeneRax vs NOTUNG, and sequence-tree vs structure-tree.

### Stage 7 — absolute dating (both approaches) + synthesis
1. **Bracket (primary).** Map the duplication branch onto a published *time-calibrated* fungal tree
   and read its age window (stem ↔ crown age of the clade it maps to). TimeTree / Shen 2018 / Li 2021.
2. **Our own relaxed clock (confirmation).** Date the Stage-4a gene tree (fixed topology) with
   **treePL** (penalized likelihood, fast) and/or **MCMCtree** (PAML, Bayesian, approx. likelihood),
   placing calibrations on speciation nodes (e.g. Asco/Basidio split; Saccharomycotina crown; the
   *Saccharomyces* whole-genome-duplication landmark ~100 Mya); the duplication-node age falls out
   with a credible interval. *Finalize specific calibrations at implementation.*
Report the node, support values, sequence-vs-structure agreement, the paralog retention/loss pattern,
and both date estimates.

---

## QC checklist (how we'll trust the answer)
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

## Decision log
- 2026-06-15: Project scaffolded. Strategy: sequence-ML-primary + reconciliation as the dating
  step; structure (FoldTree) confirmatory; AFDB-first structures; fungi-wide sampling on laptop+Colab.
- 2026-06-15: **Switched family filter PF05970 → IPR048293** after discovering PF05970 conflates
  cellular PIF1 with Helitron transposon helicases (3–10× inflation). `01` gained `--interpro`.
- 2026-06-15: **Comprehensive outgroups** chosen (all cellular Basidiomycota + Mucoromycota +
  Chytridiomycota + human PIF1). New `02_combine_cellular.py` keeps the whole 957-protein set
  (subsampling unnecessary at this size); `02_select_representatives.py` retired but kept.
- 2026-06-15: **Absolute dating = both** (bracket from a published clock + our own treePL/MCMCtree
  relaxed clock). Added Stage 7 plan.
- 2026-06-15: Sequences fetched (`04`, 957). Structures (`03`) and conda env building.
- 2026-06-15: **Gap-folding route = localcolabfold on the home RTX 4090 via WSL2** (not Colab/AF3).
  Folds only the 129 AFDB-absent proteins, AF2 to match AFDB, lean settings; ~one evening. ColabFold
  has no native-Windows GPU, hence WSL2. Steps written to `HOWTO_predict_129.md`; plan Stage 2b.
- 2026-06-15: Data acquisition complete — `pif1` env built (IQ-TREE v3.1.2 = binary `iqtree`);
  828/957 AFDB structures, 129 → ColabFold (user chose full AF2 coverage; `05` corecut → 957
  helicase cores; MAFFT FFT-NS-i + trimAl → 209-col alignment; IQ-TREE gene-tree search launched.
  New scripts `05_corecut.py`, plus `data/structures/af3/HOWTO_predict_129.md` for the predictions.
