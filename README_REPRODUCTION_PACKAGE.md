# FoldTree reproduction package — structure phylogeny of a protein family

This package reproduces the **structural** half of the PIF1/RRM3 dating project:
build a phylogeny of a protein family from predicted 3D structures rather than
from sequence alone, using **Foldseek 3Di** + **FoldTree**. It is meant to be
read and re-run end to end by someone applying the same recipe to *their own*
protein family. Every command below is the real command that produced the files
in this folder, with its actual inputs and outputs.

## 1. Overview — why a structure phylogeny

Protein 3D structure is conserved long after sequence similarity has decayed, so
for an ancient or fast-evolving family a tree built from structure can resolve
relationships that a sequence tree cannot. The approach here has four moves:

1. **Select the family** from UniProt by an InterPro signature (not just any
   Pfam domain) so you get the *cellular/functional* group and not transposon
   look-alikes, restricted to reference proteomes so AlphaFold models exist.
2. **Trawl the AlphaFold Database (AFDB)** for a predicted structure per protein
   (here: 828 of 957 had one), and **fold the gaps locally with ColabFold**
   (the remaining 129) — deliberately with AlphaFold2 to match AFDB's predictor.
3. **Core-cut** every sequence *and* every structure to the same conserved
   helicase-core residues, so the sequence tree and the structure tree describe
   the identical residues (a fair concordance check).
4. **Foldseek** does an all-vs-all 3Di+AA structural alignment; **FoldTree**
   converts those similarities to distances and infers a distance tree per
   metric, and a combined **AA+3Di partitioned IQ-TREE** is the tree that
   actually resolves the recent duplication.

For PIF1 the dataset is 957 cellular PIF1-family helicases across an even
fungi-wide taxon sample (681 Ascomycota ingroup + 235 Basidiomycota + 36
Mucoromycota + 4 Chytridiomycota + 1 human outgroup).

> **Honest caveat up front:** the three pure-distance FoldTrees
> (`fident`/`alntmscore`/`lddt`) leave the anchor PIF1/RRM3 node unresolved at
> "Fungi". Only the combined **AA+3Di ML** tree resolves the recent split. For a
> recent, high-identity split, treat the distance FoldTrees as a backbone/QC
> layer and lead with the AA+3Di tree. See §6.

## 2. What's in this folder

```
pif1-foldtree/
├── README_REPRODUCTION_PACKAGE.md   # this file
├── README.md                        # original project README (note: its Stage-4b
│                                     #   status line is STALE — see §6)
├── environment.yml                  # conda spec for the 'pif1' env
├── env.lock.txt                     # pinned versions (foldseek/fastme/iqtree/trimal/...)
│
├── workflow/                        # the numbered pipeline scripts
│   ├── 01_gather_homologs.py        # UniProt -> per-taxon cellular tables
│   ├── 02_combine_cellular.py       # concatenate -> selected.tsv (+ group/role)
│   ├── 03_fetch_structures.py       # AFDB trawl -> manifest.csv + .cif
│   ├── 04_fetch_fasta.py            # UniProt -> selected.faa
│   ├── 05_corecut.py                # hmmsearch PF05970 -> cores.faa + tip_map.tsv
│   ├── 06_integrate_predictions.py  # promote ColabFold gap rows in manifest.csv
│   ├── 07_structcorecut.py          # slice each structure to the core span
│   ├── 08_foldtree.py               # Foldseek all-vs-all + FastME -> distance trees
│   └── 13_build_aa3di.py            # map 3Di onto AA cols -> AA+3Di partition
│
├── data/
│   ├── seqs/
│   │   ├── cellular/                # 5 per-taxon tables from script 01 (asco/basidio/...)
│   │   ├── selected.tsv             # MASTER 957-protein list (18 cols incl group/role)
│   │   ├── selected.faa             # 957 full-length sequences
│   │   ├── cores.faa                # 957 core sequences, tree-safe tip labels
│   │   └── tip_map.tsv              # KEYSTONE: tip_label/accession/core_from/core_to/...
│   ├── hmm/
│   │   └── PF05970.hmm              # HMM used by 05 to find the core envelope
│   ├── 3di/                         # Foldseek 3Di DB of the cores + substitution matrix
│   │   ├── coresdb*  coresdb_ss  coresdb.lookup
│   │   └── 3di_substmat.txt         # 3Di rate matrix for the IQ-TREE 3Di partition
│   └── structures/
│       ├── afdb/                    # 828 AFDB .cif (all model_v6)
│       ├── af3/                     # 129 ColabFold .pdb  (DIR IS MISNAMED — these are AF2,
│       │                            #   not AF3) + to_predict_129.faa, HOWTO_predict_129.md,
│       │                            #   _cores_note.txt
│       └── cores/                   # 957 core-only .pdb  (THE Foldseek input; regenerable)
│
├── manifest.csv                     # 957 rows: per-protein structure provenance + path
│
├── results/
│   ├── seq_tree/
│   │   ├── aln.fasta  aln.trim.fasta        # MAFFT + trimAl AA alignment (957 x 209 trimmed)
│   │   ├── 3di.trim.fasta                   # 3Di mapped onto the trimmed AA columns
│   │   ├── alnAA_3di.fasta                  # concatenated AA+3Di (957 x 418)
│   │   ├── aa3di_partition.nex             # IQ-TREE partition (LG+G | 3Di matrix +G)
│   │   ├── pif1.treefile                    # AA-only ML tree (start tree for AA+3Di)
│   │   ├── pif1_aa3di_supported.treefile    # FINAL combined tree (1000 UFBoot + SH-aLRT)
│   │   ├── pif1_aa3di_supported.contree
│   │   ├── pif1_aa3di_supported.iqtree      # IQ-TREE report (models/supports/params)
│   │   └── corecut.domtbl                   # raw hmmsearch audit table from 05
│   └── struct_tree/
│       ├── fident_foldtree.rooted.nwk       # PRIMARY structural distance tree
│       ├── alntmscore_foldtree.rooted.nwk   # robustness check
│       └── lddt_foldtree.rooted.nwk         # robustness check
│
└── docs/
    ├── PLAN.md                      # detailed plan + the Stage-2b ColabFold recipe
    └── REVIEW_RESPONSE.md           # WHY the distance trees don't resolve the split
└── WINDOWS_HANDOFF.md, WINDOWS_HANDOFF_REFINEMENTS.md  # external-box context
```

**Not shipped (regenerable):** `results/struct_tree/res.m8` (~1.15 GB Foldseek
output), the `*_distmat.phy` / unrooted intermediate trees, and the heavy
`results/seq_tree/` bootstrap/checkpoint files. Re-running scripts 08 and 13
regenerates them. The `data/structures/cores/` folder (123 MB) is also
deterministically regenerable from `afdb/`+`af3/`+`tip_map.tsv`+`manifest.csv`
via script 07 — drop it first if you need to shrink the package.

## 3. Software prerequisites

Two execution contexts:

| Scripts | Needs | How to run |
| --- | --- | --- |
| `01`, `02`, `03`, `04`, `06` | **python3 stdlib only** + network to UniProt/EBI. No conda env, no pip, no API key. | base `python3` |
| `05`, `07`, `08`, `13` + IQ-TREE | The **`pif1` conda env**: HMMER, Foldseek, FastME, trimAl, IQ-TREE 3, biopython, dendropy, numpy, pandas. | `conda run -n pif1 ...` |
| `03` ColabFold gap-fold | **localColabFold** (`colabfold_batch`, AF2 weights, CUDA jaxlib) + an NVIDIA GPU. **Run off-repo** — see §4 step 3. | external GPU box |

Create the analysis env from the shipped spec (Apple Silicon: this env runs
osx-64 under Rosetta):

```bash
conda env create -n pif1 -f environment.yml
# or pin exactly:
conda create -n pif1 --file env.lock.txt
conda activate pif1
foldseek version        # expect 10.x
fastme --version        # expect 2.1.6.3
iqtree3 --version       # expect 3.1.2
hmmsearch -h | head -1   # HMMER on PATH
```

The repo also bundles an IQ-TREE binary at
`tools/iqtree-3.1.2-macOS-arm/bin/iqtree3` if your env's `iqtree3` is missing.
localColabFold is **not** captured in `env.lock.txt` — the gap-folding was run on
an external Windows/WSL2 RTX 4090 box and its colabfold/JAX/CUDA versions were
never recorded (a known provenance gap; see §6).

Run all commands from the repo root: `/Users/spencergray/Desktop/pif1-foldtree`.

## 4. Step-by-step walkthrough

### Step 1 — Gather the family from UniProt → per-taxon tables

```bash
for spec in 4890:ascomycota 5204:basidiomycota 1913637:mucoromycota 451435:chytridiomycota 9606:human; do
  python3 workflow/01_gather_homologs.py --taxon ${spec%%:*} --interpro IPR048293 \
    --out data/seqs/cellular/${spec##*:}.tsv
done
```

- **In:** live UniProt REST API. Filter is
  `xref:interpro-IPR048293 AND taxonomy_id:<taxon> AND keyword:KW-1185`.
- **Out:** `data/seqs/cellular/{ascomycota,basidiomycota,mucoromycota,chytridiomycota,human}.tsv`
  (681 / 235 / 36 / 4 / 1 rows).
- **Critical:** pass `--interpro IPR048293`, **not** the script default
  `--pfam PF05970`. IPR048293 (`PIF1_RRM3_pfh1`) is the *cellular* group; bare
  PF05970 also matches Helitron transposon helicases and inflates hits 3–10×.
  The `keyword:KW-1185` reference-proteome restriction (on by default) is what
  guarantees AFDB coverage downstream; `--no-refproteome` drops it but then
  structure coverage is no longer guaranteed.

### Step 2 — Combine → the master list

```bash
python3 workflow/02_combine_cellular.py --indir data/seqs/cellular --out data/seqs/selected.tsv
```

- **In:** the 5 per-taxon tables.
- **Out:** `data/seqs/selected.tsv` — **the auditable 957-protein list**, 18
  columns, with added `group` (asco/basidio/…) and `role` (ingroup/outgroup)
  columns that 03, 04, 05 all read. Don't skip this step even though it isn't a
  "focus" script.

### Step 3 — Fetch sequences

```bash
python3 workflow/04_fetch_fasta.py --in data/seqs/selected.tsv --out data/seqs/selected.faa
```

- **In:** the `Entry` column of `selected.tsv`; UniProt `stream` endpoint
  (batches of 100, 0.3 s sleep).
- **Out:** `data/seqs/selected.faa` — 957 full-length sequences with verbatim
  UniProt headers. The script verifies every requested accession came back.

### Step 4 — Sequence corecut (HMM envelope) → cores + tip_map

```bash
conda run -n pif1 python workflow/05_corecut.py
```

- **In (defaults):** `selected.faa`, `selected.tsv`, `data/hmm/PF05970.hmm`.
- **Out:** `data/seqs/cores.faa` (957 core sequences, tree-safe tip labels
  `<accession>_<Genus>`), `data/seqs/tip_map.tsv` (the **keystone** file:
  `tip_label, accession, taxid, …, core_from, core_to, full_len, core_len`),
  and `results/seq_tree/corecut.domtbl` (raw hmmsearch audit).
- **Note:** selection used the *narrow* InterPro signature; core *location* uses
  the *broader* PF05970 HMM. That's fine — finding the envelope on an already-
  filtered set can't re-introduce Helitrons. **Watch the stderr:** 05 silently
  drops any protein with no PF05970 hit at E<1e-5 (0 dropped on this run, 957 in
  → 957 out). The `core_from/core_to` columns are reused verbatim to slice the
  *structures* in step 7, which is what keeps the two trees on identical residues.

### Step 5 — Trawl AFDB for structures → manifest

```bash
python3 workflow/03_fetch_structures.py --in data/seqs/selected.tsv \
  --manifest manifest.csv --outdir data/structures/afdb --fmt cif
```

- **In:** `selected.tsv`; the AFDB prediction API
  `https://alphafold.ebi.ac.uk/api/prediction/<accession>` (one GET each).
- **Out:** `data/structures/afdb/AF-<acc>-F1-model_v6.cif` (828 hits on this run)
  + `manifest.csv` (957 rows: accession, gene, taxonomy, `structure_source`,
  `afdb_version`, `mean_plddt`, `structure_path`, …). The 129 misses are written
  with `structure_source=MISSING_predict_with_AF3`.
- **Notes:** the AFDB version is read live from the API (`latestVersion`/`cifUrl`),
  **not hardcoded** — don't assume v6 for your own run. Only HTTP 404 means
  "no model"; any other error aborts (re-run is safe — existing files are
  skipped). **Update the hardcoded `USER_AGENT` email (~line 39) to your own
  before hammering EBI.**

### Step 6 — Fold the AFDB-missing proteins with ColabFold (off-repo)

The 129 gaps were folded **outside this repo** on a Windows/WSL2 RTX 4090, with
**AF2 (ColabFold)** specifically to match AFDB's predictor — mixing AF3 would
risk FoldTree clustering by *predictor* instead of biology. The exact recipe
lives in `data/structures/af3/HOWTO_predict_129.md` and `docs/PLAN.md`:

```bash
# 1) bare-accession headers so each output is named by accession
awk '/^>/{split($0,a,"|"); print ">"a[2]; next} {print}' \
  data/structures/af3/to_predict_129.faa > data/structures/af3/to_predict_129.clean.faa

# 2) lean fold: 1 model, 3 recycles, no Amber relax, remote-MMseqs2 MSA
colabfold_batch --num-models 1 --num-recycle 3 to_predict_129.clean.faa af3_out/
#   (add --max-msa 512:1024 for >1500 aa; 3 ~2.7 kb giants were folded core-only)

# 3) rename rank_001 per protein -> data/structures/af3/<ACC>.pdb
```

The 129 resulting `.pdb` (B-factor column = per-residue pLDDT) ship in
`data/structures/af3/` — they are **irreplaceable** because the run was never
logged. Then promote them in the manifest (the **only in-repo** command for this
stage):

```bash
python3 workflow/06_integrate_predictions.py
```

- **In:** `manifest.csv` (gap rows) + `data/structures/af3/<ACC>.pdb` +
  `_cores_note.txt`.
- **Out:** `manifest.csv` rewritten in place — 129 rows become
  `structure_source=AF2_ColabFold` with `structure_path` + computed
  `mean_plddt`; a new `structure_note` column flags the 3 core-only giants
  (`core_only_AF2`). Final distribution: **828 AFDB + 129 AF2_ColabFold = 957,
  0 MISSING**. Idempotent and git-tracked (revert via git).

> The folder is named `af3/` for historical reasons only — its contents are
> **AF2/localColabFold**. Trust `structure_source` in the manifest, not the
> folder name.

### Step 7 — Structure corecut → uniform per-tip core .pdb

```bash
conda run -n pif1 python workflow/07_structcorecut.py
```

- **In (defaults):** `data/seqs/tip_map.tsv`, `manifest.csv`,
  `data/structures/afdb/` (.cif) + `data/structures/af3/` (.pdb).
- **Out:** `data/structures/cores/<tip_label>.pdb` — **957 core-only .pdb**, one
  per tip, 1:1 with `tip_map`. Each keeps only standard polymer residues whose
  author number falls in `[core_from, core_to]`, preserving original full-length
  numbering. This normalizes the mixed .cif/.pdb input into one uniform `.pdb`
  folder for Foldseek.
- **Notes:** the join key is **accession**, and the slice bounds are inherited
  from the sequence corecut (`tip_map.core_from/core_to`) — *not* recomputed —
  which is what guarantees the sequence and structure trees use identical
  residues. The 3 core-only giants (matched by `"core_only"` in
  `structure_note`) are kept **whole** because their .pdb is already renumbered
  from 1. Failures are logged to stderr only (0 on this run) — watch the
  `had no structure on disk` / `failed to parse/slice` lines on a new family.
  This script is the only authoritative documentation of the structure-slicing
  step (the original README omits it by name).

### Step 8a — Foldseek all-vs-all + FastME → structural distance FoldTrees

```bash
conda run -n pif1 python workflow/08_foldtree.py \
  --foldseek $CONDA_PREFIX/bin/foldseek --fastme $CONDA_PREFIX/bin/fastme
```

Internally this runs (you can also run by hand):

```bash
foldseek easy-search data/structures/cores data/structures/cores \
  results/struct_tree/res.m8 results/struct_tree/tmp \
  --format-output query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,lddt,lddtfull,alntmscore \
  --exhaustive-search -e 100.0 --max-seqs 2000
# then per metric:
fastme -i results/struct_tree/fident_distmat.phy -o results/struct_tree/fident_fastme.nwk -n
```

- **In:** `data/structures/cores/` (query dir == target dir → all-vs-all);
  `tip_map.tsv` to resolve `--outgroup-acc Q9H611` (human PIF1) to tip
  `Q9H611_Homo`. `easy-search` defaults to alignment-type 2 = **3Di+AA**.
- **Out:** `results/struct_tree/{fident,alntmscore,lddt}_foldtree.rooted.nwk` —
  three structural distance trees (957 tips). For each metric, similarity is
  converted to distance by the FoldTree/Tajima correction
  `d = -b·ln(1 − (1−sim)/b)` (b = 0.93 for fident, 0.95 for the others), FastME
  builds a BIONJ+NNI tree, negative branch lengths are clamped to 1e-4, and the
  tree is outgroup-rooted on `Q9H611_Homo`. The huge `res.m8` (~1.15 GB) is a
  regenerable intermediate; add `--redo` to force Foldseek to re-run.

### Step 8b — Combined AA+3Di tree (the one that resolves the split)

First build the concatenated alignment + partition:

```bash
conda run -n pif1 python workflow/13_build_aa3di.py
```

- **In:** `results/seq_tree/aln.fasta` + `aln.trim.fasta` (your MAFFT+trimAl AA
  alignment — see note below), `data/3di/coresdb` (Foldseek 3Di states of the
  cores), `data/3di/3di_substmat.txt`.
- **Out:** `results/seq_tree/3di.trim.fasta` (3Di mapped onto the trimmed AA
  columns, identical gap pattern), `alnAA_3di.fasta` (957 × 418 = 209 AA + 209
  3Di), `aa3di_partition.nex` (charset `AA=1-209` LG+G, `TDi=210-418` with the
  3Di matrix +G).

Then infer the tree (these IQ-TREE commands are **not scripted** — they were run
by hand; recorded in `pif1_aa3di_supported.iqtree`/`.log`):

```bash
# final supported tree: 1000 UFBoot + 1000 SH-aLRT
iqtree3 -s results/seq_tree/alnAA_3di.fasta -p results/seq_tree/aa3di_partition.nex \
  -B 1000 -alrt 1000 -bnni -T AUTO -seed 42 \
  -pre results/seq_tree/pif1_aa3di_supported -redo
```

- **Out:** `results/seq_tree/pif1_aa3di_supported.treefile` (+ `.contree`,
  `.iqtree`) — **the final combined tree** carried into reconciliation, and the
  only tree that resolves the PIF1/RRM3 duplication.

> **Provenance gaps in step 8b you must fill for your own family:** no numbered
> script builds `data/3di/coresdb` or stages `3di_substmat.txt` — they were made
> off-repo. Build the DB with
> `foldseek createdb data/structures/cores data/3di/coresdb`, and get the 3Di
> substitution matrix from the DessimozLab `fold_tree` repo / FoldTree paper.
> The AA alignment (`aln.fasta`/`aln.trim.fasta`) comes from running MAFFT on
> `cores.faa` then trimAl (`-automated1`, `-colnumbering`) — that AA-tree
> sub-pipeline lives upstream of this structural handoff.

## 5. Adapting this to your own protein family

The recipe is family-agnostic; you change a small number of identifiers:

1. **Family selection signature (step 1)** — replace `--interpro IPR048293`
   with the InterPro entry specific to your cellular/functional group. This is
   the #1 decision: pick an InterPro signature that excludes paralog/transposon
   look-alikes, not just any Pfam domain. Keep `keyword:KW-1185` on if you want
   guaranteed AFDB coverage.
2. **Core-cut HMM (step 4)** — replace `data/hmm/PF05970.hmm` with the Pfam HMM
   for your family's conserved domain, and **record where/when you fetched it**
   (the PF05970 provenance here is under-documented — don't repeat that mistake).
   The HMM only locates the envelope, so the broader Pfam model is fine.
3. **Rooting anchor (steps 7–8)** — replace the human-PIF1 outgroup accession
   `Q9H611` (`--outgroup-acc`) with your family's outgroup, and make sure that
   accession is in your `selected.tsv` so it has a tip.
4. **Input column / manifest schema** — keep the UniProt-accession column named
   `Entry` (hardcoded in 03), and `manifest.csv` keyed by accession with
   `structure_path` + `structure_note`. The taxonomy columns
   (`subphylum/class/order/family/genus`) are read with `.get()`, so they're
   optional, but `Entry` is required.
5. **Tip labels** — labels are `<accession>_<Genus>`, uniquified on genus
   collision. Join on the **accession** column, never the label.
6. **Predictor consistency** — fold your gaps with the same predictor AFDB used
   (AF2 via ColabFold), for the reason in §6.

Everything else (scripts 02, 04, 05, 06, 07, 08, 13) runs unchanged once those
identifiers are set.

## 6. Caveats

- **AFDB-vs-ColabFold predictor confound.** 828 structures are AFDB (AF2), 129
  are local ColabFold. They were *both* folded with AlphaFold2 on purpose, so
  the tree can't cluster proteins by predictor. If you fold gaps with a
  *different* model (e.g. AF3), you reintroduce that confound — don't.
- **The ColabFold run is not reproducible from this repo.** It ran on an
  external RTX 4090/WSL2 box; the repo has the recipe, the input FASTA, the 129
  output PDBs, and the integration step, but **no run log and no
  colabfold/JAX/CUDA versions**.
- **Distance FoldTrees don't resolve the recent split.** All three
  `fident`/`alntmscore`/`lddt` → FastME trees leave the PIF1/RRM3 anchor at
  "Fungi". Only the **AA+3Di ML partition** tree resolves the duplication
  (`docs/REVIEW_RESPONSE.md`). For a recent, high-identity split, lead with the
  combined tree and use the distance trees as a QC/backbone layer.
- **No explicit pLDDT filtering.** Every structure is used regardless of model
  confidence (AFDB pLDDT here: min 52 / median 65 / max 88). Low-confidence
  predictions still feed Foldseek. If your family has many poorly-modeled
  members, consider filtering on `mean_plddt` in `manifest.csv` before step 7.
- **Single-outgroup rooting.** The distance trees are rooted on the *one* human
  PIF1 tip (`Q9H611_Homo`); the AA+3Di tree is inferred unrooted and rooted
  downstream. A single outgroup is fragile — a misplaced or long-branch outgroup
  can mislead the root. Add more outgroups if you can.
- **`data/3di/` and the AA alignment have no generating script in this package.**
  `coresdb` and `3di_substmat.txt` were staged off-repo (build/obtain them as in
  §4 step 8b), and the two AA+3Di IQ-TREE commands were run by hand (recorded in
  `pif1_aa3di_supported.iqtree`).
- **The original `README.md` is stale on Stage 4b** — it says "FoldTree (Colab),
  not done". It *was* done locally via scripts 08/13 + manual IQ-TREE. Trust the
  scripts and `results/`, not that status line.
- **Live-API drift.** Steps 1, 3, 5 hit live UniProt/AFDB; results drift as those
  databases update. Repo snapshot date: **2026-06-15**. Record your own run date.

---

**Expected counts (verified on disk):** `selected.tsv` 957 rows · `selected.faa`
957 · `cores.faa` 957 · `tip_map.tsv` 957 · `manifest.csv` 957 (828 AFDB + 129
AF2_ColabFold) · `data/structures/afdb/` 828 `.cif` · `data/structures/af3/` 129
`.pdb` · `data/structures/cores/` 957 `.pdb`.
