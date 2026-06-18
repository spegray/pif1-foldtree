# Windows / WSL2 handoff — GeneRax reconciliation + supported AA+3Di tree

**You are a fresh Claude Code session running on a Windows desktop (via WSL2, native x86-64 Linux).**
You have no memory of the prior work; everything you need is in this repo and this file. A parallel
session on a Mac did all the upstream work and is waiting for your results. Read this whole file,
then execute the tasks. Ask the user before anything destructive; otherwise proceed.

---

## The project in one paragraph
We are dating the duplication that produced the two budding-yeast PIF1-family helicases, **PIF1**
and **RRM3** — i.e. on which branch of the fungal species tree the ancestral gene duplicated. Method:
gather all cellular PIF1-family proteins across a fungi-wide taxon sample (957 proteins / 728 species,
Ascomycota ingroup + Basidiomycota / early-diverging-fungi / human outgroups), build a **gene tree**,
and **reconcile** it against a **species tree** — reconciliation maps the duplication node onto a
species-tree branch, and *that branch is the answer*. (Background: the FoldTree paper, Moi et al.
*Nat Struct Mol Biol* 2025.)

## What's already done (on the Mac) and the answer so far
- 957 proteins, structures (828 AlphaFold DB + 129 ColabFold), all in `manifest.csv`.
- **Sequence gene trees built:** plain `LG+I+G4` (`results/seq_tree/pif1.treefile`) and an **AA+3Di
  partitioned tree** (`results/seq_tree/pif1_aa3di.treefile`) that adds the 3Di structural alphabet.
- **NCBI species tree** (`data/species_tree/ncbi_species.nwk`, 720 spp, binary) + gene→species map
  (`data/species_tree/gene_species.map`, 941 genes).
- **Key finding:** a fast species-overlap reconciliation places the PIF1/RRM3 duplication at
  **Saccharomycotina**, *but only on the AA+3Di tree* — the AA-only trees couldn't resolve the deep
  node (the ScPif1/ScRrm3 ancestor mis-mapped to "Fungi" because the 209-aa core has saturated). The
  3Di structural signal fixed that. The two-copy **mushrooms (Agaricomycetes) are an independent
  duplication**, not shared ancestry.

## Why you (Windows) are running this
**GeneRax** — the gold-standard maximum-likelihood duplication–loss reconciliation tool — **segfaults
under Rosetta on the Mac** (it's an x86 binary emulated on Apple Silicon). On your native x86-64 Linux
it runs fine. Your job is to run GeneRax to **rigorously confirm** the Saccharomycotina placement,
give the **gene-loss pattern** (which lineages kept both paralogs vs one), and pin **stem vs crown**.
Optionally, rebuild the AA+3Di tree **with bootstrap support** (fast on this hardware).

---

## Step 1 — environment (WSL2, native linux-64 conda)
In a WSL2 Ubuntu shell. If you don't have conda/mamba, install Miniforge first:
```bash
# (skip if `conda` already works in WSL2)
wget -q https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
bash Miniforge3-Linux-x86_64.sh -b -p "$HOME/miniforge3" && eval "$($HOME/miniforge3/bin/conda shell.bash hook)"

# native-linux env with the tools that were broken/slow under Rosetta:
conda create -y -n pif1recon -c conda-forge -c bioconda \
    generax iqtree dendropy ete3 biopython pandas numpy
conda activate pif1recon
generax --help | head -3        # smoke test — should print, not segfault
```

## Step 2 — get the data
The repo is **github.com/spegray/pif1-foldtree** (private). It was cloned here earlier but is stale:
```bash
cd ~/pif1-foldtree 2>/dev/null && git pull        # if it exists
# else: gh auth login   (or use your token), then:
#   git clone https://github.com/spegray/pif1-foldtree.git ~/pif1-foldtree && cd ~/pif1-foldtree
```
Confirm the reconciliation inputs are present (pre-built on the Mac, paths are repo-relative):
```bash
ls results/reconciliation/aa3di/   # families.txt, genetree.pruned.nwk, aln.pruned.fasta
ls data/species_tree/ncbi_species.nwk data/3di/3di_substmat.txt
```

## Step 3 — TASK A (primary): GeneRax reconciliation
**Run from the repo root** (the families.txt paths are relative). Do the AA+3Di tree first (the one
that resolved the deep node), then the LG+I+G4 tree for comparison. `EVAL` reconciles the tree as-is;
`SPR` additionally lets GeneRax correct weakly-supported branches.
```bash
# AA+3Di tree, reconcile as-is:
generax -f results/reconciliation/aa3di/families.txt -s data/species_tree/ncbi_species.nwk \
        -r UndatedDL --unrooted-gene-tree --geneSearchStrategy EVAL -p results/reconciliation/aa3di/eval
# AA+3Di tree, with gene-tree correction:
generax -f results/reconciliation/aa3di/families.txt -s data/species_tree/ncbi_species.nwk \
        -r UndatedDL --unrooted-gene-tree --geneSearchStrategy SPR  -p results/reconciliation/aa3di/spr
# LG+I+G4 tree (comparison):
generax -f results/reconciliation/lgig/families.txt  -s data/species_tree/ncbi_species.nwk \
        -r UndatedDL --unrooted-gene-tree --geneSearchStrategy SPR  -p results/reconciliation/lgig/spr
```
(If a bare `generax` misbehaves, try `mpiexec -np 1 generax ...`.)

**Interpret the result** — two complementary reads:
1. `*/reconciliations/pif1_speciesEventCounts.txt` and `pif1_eventCounts.txt` give per-species-tree-node
   D/S/L counts. Find the node(s) with a Duplication (D) — that's where PIF1/RRM3 split — and note the
   loss (L) pattern.
2. Run our species-overlap reporter on the GeneRax-corrected gene tree to name the clade in plain terms:
   ```bash
   python workflow/12_reconcile.py --tree results/reconciliation/aa3di/spr/results/pif1/geneTree.newick
   ```
   (ete3 will download the NCBI taxdump on first run — a few minutes, once.) **Expected:** the
   anchor MRCA(Pif1,Rrm3) maps to **Saccharomycotina**, confirming the Mac result with the gold-standard tool.

## Step 4 — TASK B (optional, fast here): AA+3Di tree WITH support values
The Mac built the AA+3Di tree topology-only (UFBoot is incompatible with its `-n` search cap). On this
hardware you can afford a full search + 1000 ultrafast bootstraps. **Run from repo root:**
```bash
iqtree2 -s results/seq_tree/alnAA_3di.fasta -p results/seq_tree/aa3di_partition.nex \
        -B 1000 -alrt 1000 -bnni -T AUTO -seed 42 -pre results/seq_tree/pif1_aa3di_supported -redo
# then reconcile the supported tree:
python workflow/14_prep_generax.py --genetree results/seq_tree/pif1_aa3di_supported.treefile --label aa3di_supp
generax -f results/reconciliation/aa3di_supp/families.txt -s data/species_tree/ncbi_species.nwk \
        -r UndatedDL --unrooted-gene-tree --geneSearchStrategy SPR -p results/reconciliation/aa3di_supp/spr
```

## Step 5 — report back to the Mac session
Write a short `results/reconciliation/WINDOWS_RESULTS.md` capturing, for each run:
- the species-tree node / clade the duplication maps to (Saccharomycotina? stem or crown?),
- the D and L event totals (`pif1_eventCounts.txt`),
- whether AA+3Di and LG+I+G4 agree, and whether GeneRax confirms the species-overlap result,
- the UFBoot support on the duplication node (Task B), if run.

Then push everything so the Mac can pull it:
```bash
git add results/reconciliation WINDOWS_HANDOFF.md results/seq_tree/pif1_aa3di_supported.* 2>/dev/null
git commit -m "GeneRax reconciliation + supported AA+3Di (Windows/WSL2 native run)"
git push
```
Tell the user: *"Done — push complete; go back to the Mac session and pull."*

---
**One-line summary of your mission:** run GeneRax (which the Mac can't) to confirm with a
gold-standard tool that the PIF1/RRM3 duplication maps to **Saccharomycotina**, get the gene-loss
pattern and stem-vs-crown precision, optionally add bootstrap support, and push the results back.
