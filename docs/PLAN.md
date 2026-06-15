# Plan: Dating the PIF1 / RRM3 duplication in Ascomycota

## Context

**Question:** In budding yeast, when did the two PIF1-family helicases — *PIF1* and *RRM3* —
split from a single ancestral gene? Equivalently: on which branch of the fungal species tree
did the gene duplication that produced the PIF1 and RRM3 lineages occur (and, optionally, how
long ago)?

**Why this plan reshapes the original approach.** The original method was: gather Ascomycota
genomes → pull every PIF1-family helicase → predict all structures with AlphaFold 3 → build
traditional + structural phylogenetic trees and "hone in." After reading the attached FoldTree
paper (Moi et al., *Nat Struct Mol Biol* 2025, doi:10.1038/s41594-025-01649-8) and confirming
the `pif1-foldtree` working directory is empty (clean start), three corrections make the project
far more tractable and the answer far more defensible:

1. **This is a gene-duplication *dating* problem, which is a reconciliation task — not a
   "build a tree and look" task.** The answer comes from reconciling a PIF1-family **gene tree**
   against a trusted fungal **species tree**: reconciliation maps the duplication node onto a
   specific species-tree branch. That branch *is* the answer, and reconciliation guards against
   the copy-counting trap (a 1-copy species may have *lost* a paralog, not predated the split).

2. **Sequence ML is the workhorse; structure is confirmatory — not the reverse.** The paper's own
   result (Fig. 2d, Discussion): structural phylogenetics wins for *deep, twilight-zone* (<~20–25%
   identity) relationships, but for *closely related, high-identity* homologs a sequence ML tree is
   as good or better and FoldTree "may suffer in quality closer to the leaves." The PIF1/RRM3 split
   inside budding yeast is recent and high-identity → lead with sequence ML, use the structural tree
   to confirm and to stabilize the *deep* PIF1-family backbone and rooting.

3. **We barely need AlphaFold 3.** FoldTree consumes the **AlphaFold Database** directly. Sampling
   from UniProt **reference proteomes** guarantees precomputed AFDB structures we can bulk-download
   by accession. AF3 (via ColabFold / AF3 server) is reserved only for the handful of genome-mined
   proteins absent from UniProt.

**Decisions locked with the user:** Ascomycota-wide taxon sampling · AFDB-first structures · the
**129 AFDB gaps folded locally with AF2/ColabFold (localcolabfold on the home RTX 4090 via WSL2)** to
match AFDB's predictor · sequence ML + reconciliation primary, structural FoldTree confirmatory.

**Primary deliverable:** the species-tree branch (ancestral node) on which the PIF1→{PIF1, RRM3}
duplication is placed, with support values and a sequence-vs-structure agreement check.
**Optional secondary deliverable:** an absolute-time estimate for that branch via a calibrated
fungal time tree.

---

## Project layout (to be created in `~/Desktop/pif1-foldtree`)

A single auditable tree with a manifest spreadsheet as the source of truth:

```
pif1-foldtree/
  README.md                # running lab-notebook: every command, version, date, decision
  environment.yml          # conda/mamba env, version-pinned (reproducibility)
  manifest.csv             # ONE row per protein: species, clade, genome source, UniProt acc,
                           #   AFDB path/AF3?, gene name, paralog label, domain coords
  data/
    species_tree/          # published fungal species tree, pruned to our taxa
    proteomes/             # reference proteomes (or links) per species
    seqs/                  # gathered PIF1-family protein FASTAs
    structures/            # AFDB .pdb/.cif (+ any AF3 predictions), and corecut cores
  results/
    seq_tree/              # MSA, trimmed alignment, IQ-TREE outputs
    struct_tree/           # fold_tree outputs (Foldseek matrix, NJ tree)
    reconciliation/        # GeneRax / NOTUNG outputs → the duplication node
    figures/
  workflow/                # Snakemake / scripts wrapping each stage
```

---

## Stage 0 — Define the taxon set and the species-tree backbone

Goal: an even, Ascomycota-wide sample with guaranteed AFDB coverage and a trusted species tree.

- **Backbone tree:** prune a published genome-scale fungal phylogeny to our taxa rather than
  inferring one from scratch (defensible + laptop-friendly). Candidates: Li et al. 2021
  *Curr Biol* "A genome-scale phylogeny of the kingdom Fungi"; for the budding-yeast portion,
  the Y1000+ / Shen et al. 2018 *Cell* budding-yeast tree. *Verify current best reference at
  implementation time.*
- **Even Ascomycota sampling** across the three subphyla, prioritizing species that are UniProt
  **reference proteomes** (→ AFDB structures exist):
  - *Saccharomycotina* (budding yeasts) — dense, since the split is hypothesized here: e.g.
    *S. cerevisiae*, *S. uvarum/paradoxus*, *Naumovozyma*, *Kluyveromyces*, *Lachancea*,
    *Zygosaccharomyces*, *Candida/Debaryomyces* (CTG clade), *Yarrowia* (early-diverging).
  - *Pezizomycotina* (filamentous) — *Aspergillus*, *Neurospora*, *Fusarium*, *Botrytis*,
    *Zymoseptoria*, etc., spanning the major classes.
  - *Taphrinomycotina* — *Schizosaccharomyces pombe/japonicus* (single PIF1 = Pfh1),
    *Pneumocystis*, *Taphrina*.
  - **Outgroups for rooting:** ≥2 Basidiomycota + ≥1 early-diverging fungus (Mucoromycota /
    chytrid). Optionally one animal/plant PIF1 to anchor the deep root.
- Record every species + genome source + UniProt proteome ID in `manifest.csv`.
- **Scale control:** target ~80–150 taxa. Density matters most *around the candidate node*
  (Saccharomycotina); the rest provides even bracketing + rooting.

## Stage 1 — Gather PIF1-family proteins (automated, with an orthology cross-check)

- **Seeds:** *S. cerevisiae* Pif1 (UniProt **P07271**) and Rrm3 (UniProt **P38766**).
- **Primary search:** `hmmsearch` with the PIF1 helicase HMM against each proteome.
  *Verify the exact Pfam/InterPro model at implementation time — likely Pfam PF05970 "PIF1"
  / InterPro IPR010285 "DNA helicase Pif1-like"; do NOT assume.* Alternative if HMM is awkward:
  `jackhmmer`/`phmmer` or BLASTp from the two seeds, then confirm each hit carries the PIF1
  domain via `hmmscan` against Pfam.
- **Independent cross-check:** pull precomputed orthogroups from **OMA** (same lab as FoldTree),
  **OrthoDB**, or **eggNOG** for Pif1/Rrm3. Agreement between de-novo search and precomputed
  orthology is a built-in QC. Disagreements get flagged in `manifest.csv`.
- **De-duplicate:** one representative protein per gene (collapse isoforms; note that yeast *PIF1*
  itself produces nuclear + mito forms from one gene via alternative start — that is a single gene).
- **Label paralogs provisionally** (PIF1-like vs RRM3-like) but treat labels as hypotheses to be
  confirmed by the tree, not inputs to it.

## Stage 2 — Map to UniProt and fetch structures (AFDB first)

- Map every gathered protein to a UniProt accession via the **UniProt ID-mapping** service.
- Bulk-download AFDB structures by accession via the **prediction API**
  (`https://alphafold.ebi.ac.uk/api/prediction/<acc>` → current model version **v6**, not the v4 in
  older docs). `workflow/03_fetch_structures.py` already did this: **828/957 from AFDB**.
- **Gaps:** the **129** proteins with no AFDB model (recent accessions AFDB hasn't computed) are
  tagged `MISSING_predict_with_AF3` in `manifest.csv` and their sequences are already extracted to
  `data/structures/af3/to_predict_129.faa`. We fold *only these 129* locally (Stage 2b) — the 828
  AFDB structures are reused untouched.
- Record structure source (AFDB version vs AF2_ColabFold) per protein in `manifest.csv` — provenance for audit.

## Stage 2b — Fold ONLY the 129 AFDB gaps locally (WSL2 + localcolabfold, RTX 4090)

Goal: complete structural coverage for the confirmatory FoldTree (Stage 4b) **without** re-folding
anything AFDB already has, and with the **same predictor family (AlphaFold2)** as AFDB so the
structural tree never clusters proteins by predictor. ColabFold = AF2 with a fast remote-MMseqs2 MSA.

**Why local on the 4090, not Colab/AF3:** 129 monomers exceed a single free-Colab session (~12 h cap,
disconnects); the AF3 server is AF2-inconsistent and capped at ~30 jobs/day (~5 days). The 4090 runs
the whole batch unattended in **one evening (~3–6 h)**. (See the AF2-vs-AF3 reasoning recorded in
`data/structures/af3/HOWTO_predict_129.md`.)

**One-time WSL2 + localcolabfold setup** (ColabFold has *no native-Windows GPU* — JAX needs Linux;
WSL2 gives the 4090 a Linux environment with full CUDA passthrough):

1. **WSL2 + Ubuntu** — in an *admin* PowerShell: `wsl --install` (installs Ubuntu), reboot. Needs
   Windows 11 or Windows 10 21H2+.
2. **GPU driver** — install the latest NVIDIA Windows driver (Studio/Game-Ready ≥ R525). Do **not**
   install a CUDA driver *inside* WSL; the Windows driver exposes the GPU to WSL automatically.
3. **Verify GPU in WSL** — open Ubuntu, run `nvidia-smi` → it must list the RTX 4090. (If not: update
   the Windows driver and `wsl --update` in PowerShell.)
4. **Install localcolabfold** (handles the CUDA-enabled JAX + downloads AF2 weights):
   ```bash
   sudo apt update && sudo apt install -y wget git
   wget https://raw.githubusercontent.com/YoshitakaMo/localcolabfold/main/install_colabbatch_linux.sh
   bash install_colabbatch_linux.sh
   export PATH="$HOME/localcolabfold/colabfold-conda/bin:$PATH"   # use the exact path the installer prints
   colabfold_batch --help      # smoke test
   ```
   The 4090 is Ada (sm_89) → needs a recent CUDA-12 jaxlib, which the current installer provides. A
   `no kernel image is available` error = stale install → re-run the installer / update the driver.

**Fold the gaps (only the 129):**

5. **Bare-accession headers** so every output is named by accession (UniProt headers are
   `>tr|ACC|...`). One-liner:
   ```bash
   awk '/^>/{split($0,a,"|"); print ">"a[2]; next} {print}' \
       data/structures/af3/to_predict_129.faa > data/structures/af3/to_predict_129.clean.faa
   ```
6. **Run the batch** (lean settings — FoldTree's Foldseek 3Di only needs one good backbone; 5 models
   and Amber relax are wasted cost). Working on the WSL filesystem (`~/…`) is faster than `/mnt/c/…`:
   ```bash
   colabfold_batch --num-models 1 --num-recycle 3 \
       to_predict_129.clean.faa  af3_out/
   # default msa_mode = remote MMseqs2 server (no local 2 TB DB); no --amber (relax off)
   ```
   - **115 gaps are ≤1500 aa → fold fine** on 24 GB (~1–3 min each).
   - **14 gaps are >1500 aa**; the **two ~2.7 kb *Candidozyma* outliers (A0ABX8I2X5, A0A2V1AIT0)**
     will likely OOM at 24 GB. For the long ones, add `--max-msa 512:1024` (cuts memory); if any
     still OOM, fold **only their HMM helicase-core region** (slice with the Stage-3 envelope) and set
     `structure_note=core_only` in `manifest.csv`. 14/957 with a documented caveat is negligible.
7. **Keep the top model per protein** — colabfold writes `<ACC>_unrelaxed_rank_001_*.pdb`; copy each to
   `data/structures/af3/<ACC>.pdb`.
8. **Integrate** — add `workflow/06_integrate_predictions.py`: scan `data/structures/af3/<ACC>.pdb`,
   normalize names, and update each gap's `manifest.csv` row (`structure_source=AF2_ColabFold`,
   `structure_path`, `mean_plddt` = mean Cα B-factor). These `.pdb` then join the 828 AFDB `.cif` for
   the **structural** corecut feeding Stage 4b (Foldseek reads both `.cif` and `.pdb`, so the mix is fine).

**Audit:** record `colabfold_batch --version` + JAX/CUDA versions alongside `env.lock.txt`, and the run
date, in the README decision log.

## Stage 3 — Trim to the conserved helicase core ("corecut")

Goal: stop variable N/C-terminal extensions from confounding either tree (the paper built the
"corecut" pipeline, Fig. 4, for exactly this).

- **Sequences:** trim each protein to the PIF1/SF1 helicase-core boundaries from the HMM alignment.
- **Structures:** trim each structure to the matching residues (fold_tree's corecut step, or a
  short Biopython script using the HMM-derived coordinates). Cluster N/C termini separately if you
  want to annotate architecture changes (optional, paper-style).

## Stage 4a — Sequence ML tree (the workhorse)

- Align cores: **MAFFT** `--maxiterate 1000 --localpair` (L-INS-i) or **MUSCLE5**.
- Trim: **trimAl** `--automated1`.
- Infer: **IQ-TREE 2** with ModelFinder (`-m MFP`; expect LG+G/LG+I+G, consider a profile-mixture
  model for the deeper nodes), `-B 1000` ultrafast bootstrap + `-alrt 1000` SH-aLRT.
- Runs comfortably on a laptop for ~100–300 single-domain sequences.

## Stage 4b — Structural tree (confirmatory) via FoldTree

- Use the **fold_tree** pipeline (Snakemake or the Colab notebook):
  - GitHub: `https://github.com/DessimozLab/fold_tree`
  - Colab: `https://colab.research.google.com/github/DessimozLab/fold_tree/blob/main/notebooks/FoldTree.ipynb`
- Pipeline: Foldseek all-vs-all → **Fident** (3Di+AA, alignment mode 1) with the paper's statistical
  correction → distance tree (QuickTree/FastME) → **MAD** rooting.
- **Best-of-both option** (recommended for the deep backbone): extract 3Di sequences, concatenate
  AA+3Di, run **IQ-TREE partitioned** (LG matrix for AA partition + 3Di substitution matrix). Per
  the paper this combines sequence + structure signal and helps the deepest nodes.

## Stage 5 — Species tree, pruned to our taxa

- Prune the published backbone (Stage 0) to exactly our species (Python `ete3` or R `ape`).
- *Fallback if a published tree lacks some taxa:* build from concatenated single-copy orthologs
  (**BUSCO** `fungi_odb10`) with IQ-TREE — heavier, only if needed.

## Stage 6 — Reconciliation = the actual answer

Reconcile the PIF1 gene tree (Stage 4a, confirmed by 4b) against the species tree (Stage 5):

- **NOTUNG** (primary; GUI, approachable for a non-bioinformatician): duplication/loss parsimony,
  rearranges weakly-supported gene-tree branches, visually maps the duplication node onto the
  species tree → read off the branch.
- **GeneRax** (robustness check; CLI): joint ML gene-tree + reconciliation under a
  duplication-loss(-transfer) model; accounts for gene-tree error.
- **ALE** (optional): integrates over gene-tree uncertainty from a bootstrap/posterior sample.
- **Output:** the species-tree branch carrying the PIF1→{PIF1, RRM3} duplication. Cross-check that
  NOTUNG and GeneRax agree, and that the sequence- and structure-derived gene trees place it on the
  same branch. Examine the loss pattern (which lineages kept both vs one paralog).

## Stage 7 — Optional absolute dating + synthesis

- **Relative answer (primary):** "the duplication maps to the ancestor of clade X."
- **Absolute time (optional):** map that branch onto a calibrated fungal time tree (**TimeTree**,
  or published fungal molecular-clock estimates) to attach an approximate MYA range. Full
  node-dating (e.g. BEAST/MCMCtree) is out of scope unless you want a defended absolute age.
- **Synthesis:** report the node, support values on the duplication and bracketing nodes,
  sequence-vs-structure agreement, and the paralog-retention/loss pattern.

---

## Scalability / automation levers (addresses "surprisingly large")

- **Reuse precomputed orthology** (OMA/OrthoDB/eggNOG) instead of ad-hoc all-vs-all BLAST.
- **Sample UniProt reference proteomes** → AFDB coverage is automatic → AF3 nearly eliminated.
- **Wrap each stage in Snakemake** (fold_tree already is) + a `manifest.csv` source of truth →
  one-command, re-runnable, auditable; adding taxa later = add rows + re-run.
- **Pin versions** in `environment.yml`; log every command/date in `README.md`.

## Verification / QC (how we'll know the answer holds)

1. **Orthology agreement:** de-novo HMM search vs OMA/OrthoDB membership match (flag mismatches).
2. **Tree support:** UFBoot/SH-aLRT ≥ conventional thresholds on the duplication node and the two
   bracketing nodes.
3. **Sequence vs structure concordance:** both gene trees place the duplication on the same
   species-tree branch (this is the key robustness signal given the paper's leaf-level caveat).
4. **Reconciliation concordance:** NOTUNG and GeneRax agree on the duplication branch.
5. **Loss-pattern sanity:** retention/loss across lineages is biologically coherent (e.g. fission
   yeast single Pfh1 = pre-duplication or independent loss, resolved by the tree, not assumed).
6. **Rooting robustness:** answer is stable to outgroup choice and to MAD vs outgroup rooting.

## Things to verify at implementation time (do not assume)

- Exact PIF1 Pfam/InterPro model ID (likely PF05970 / IPR010285 — confirm).
- Current best published fungal species tree to prune (Li et al. 2021 vs newer; Y1000+ for yeasts).
- fold_tree repo path / Colab notebook URL still current; Foldseek + IQ-TREE current versions.
- Seed accessions P07271 (Pif1) / P38766 (Rrm3) still current in UniProt.
- AFDB file version (`_v4` vs newer) at download time.

## Open follow-ups for the user (not blocking the plan)

- Do you want the optional **absolute-time** estimate, or is the phylogenetic placement enough?
- Any specific Ascomycota lineages you want guaranteed in the sample (beyond an even spread)?
