# Predicting the 129 missing structures

129 of the 957 cellular PIF1 proteins have **no AlphaFold DB model** (mostly very recent
UniProt accessions AFDB hasn't computed). Input sequences are in **`to_predict_129.faa`**
(this folder). The sequence ML tree (Stage 4a) already uses all 957; these predictions are
only needed for the confirmatory FoldTree (Stage 4b) at full taxon coverage.

The 828 AFDB structures are **AlphaFold2**. Predict the gaps with the *same* method so the
FoldTree doesn't cluster proteins by predictor. ColabFold = AF2 with a fast MMseqs2 MSA.

## Recommended: localcolabfold on the home RTX 4090 (via WSL2)

129 monomers exceed one free-Colab session (~12 h cap, disconnects), and the AF3 server is
AF2-inconsistent and capped at ~30 jobs/day (~5 days). The 4090 runs the whole batch unattended
in **one evening (~3–6 h)**. ColabFold has **no native-Windows GPU support** (JAX needs Linux), so
run it under **WSL2**, which gives the 4090 a Linux env with full CUDA passthrough.

### One-time setup
1. **WSL2 + Ubuntu** — admin PowerShell: `wsl --install`, reboot (Windows 11, or Win10 21H2+).
2. **GPU driver** — latest NVIDIA *Windows* driver (≥ R525). Do **not** install a CUDA driver inside
   WSL; the Windows driver exposes the GPU automatically.
3. **Verify** — in Ubuntu: `nvidia-smi` must list the RTX 4090. (If not: update driver; `wsl --update`.)
4. **Install localcolabfold** (sets up CUDA-enabled JAX + downloads AF2 weights):
   ```bash
   sudo apt update && sudo apt install -y wget git
   wget https://raw.githubusercontent.com/YoshitakaMo/localcolabfold/main/install_colabbatch_linux.sh
   bash install_colabbatch_linux.sh
   export PATH="$HOME/localcolabfold/colabfold-conda/bin:$PATH"   # exact path is printed by the installer
   colabfold_batch --help          # smoke test
   ```
   The 4090 is Ada (sm_89) → needs a recent CUDA-12 jaxlib (the current installer provides it). A
   `no kernel image is available` error means a stale install → re-run the installer / update the driver.

### Fold the 129 gaps
5. **Bare-accession headers** (so every output is named by its accession):
   ```bash
   awk '/^>/{split($0,a,"|"); print ">"a[2]; next} {print}' \
       to_predict_129.faa > to_predict_129.clean.faa
   ```
6. **Run the batch** (lean: one good backbone is all Foldseek's 3Di needs). Copy the FASTA into the
   WSL filesystem (`~/…`) — much faster I/O than `/mnt/c/…`:
   ```bash
   colabfold_batch --num-models 1 --num-recycle 3  to_predict_129.clean.faa  out/
   # default msa_mode = remote MMseqs2 server (no local DB); no --amber (relaxation off)
   ```
   - **115 gaps ≤1500 aa fold fine** on 24 GB (~1–3 min each).
   - **14 gaps >1500 aa**; the **two ~2.7 kb *Candidozyma* outliers (`A0ABX8I2X5`, `A0A2V1AIT0`)**
     will likely OOM at 24 GB. For the long ones add `--max-msa 512:1024`; if any still OOM, fold
     **only the helicase-core region** (the Stage-3 HMM envelope) and set `structure_note=core_only`
     in `manifest.csv`.
7. Keep the **top model** per protein (`<ACC>_unrelaxed_rank_001_*.pdb`) — see naming below.

## Fallback A: ColabFold batch notebook (Colab, no install)
Use only if WSL2 is unavailable. https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/batch/AlphaFold2_batch.ipynb
→ Runtime → GPU (T4). Upload `to_predict_129.faa`; settings `msa_mode=mmseqs2_uniref_env`,
`num_models=1`, `num_recycles=3`, `use_amber=False`. **Split into ~3 chunks** (the notebook resumes,
skipping entries already in the output dir).

## Naming so the pipeline can ingest the results

Drop the top model for each protein in **this folder** named **`<ACCESSION>.pdb`**
(e.g. `A0ABN8WKL1.pdb`). A small `workflow/06_integrate_predictions.py` (added once files
exist) will scan this folder, copy/standardize the files, and update `manifest.csv`
(`structure_source = AF2_ColabFold`, `structure_path`, pLDDT from the b-factor column).

## Alternative: AlphaFold3 web server (alphafoldserver.com)

Works, but: requires login, ~**30 jobs/day** (≈5 days for 129), and **AF3 ≠ AF2** — mixing
AF3 predictions with the AF2 AFDB set is a mild methodological inconsistency for FoldTree.
If you go this route, save each result as `<ACCESSION>.cif` here. (Optional sanity check:
re-predict a handful of AFDB proteins with AF3 too, to confirm predictor choice doesn't move
their placement in the structural tree.)

## Note
Full structural coverage is a *confirmatory* layer. If predicting all 129 stalls, the FoldTree
on the 828 we already have is scientifically sufficient; the primary dating answer comes from
the all-957 sequence tree + reconciliation.
