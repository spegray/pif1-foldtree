#!/usr/bin/env python3
"""
18_lba_outgroup_removal.py  --  LBA diagnostic: does the amino-acid tree recover the Saccharomycotina
Pif1/Rrm3 duplication once the long-branch outgroups are removed?

Reviewer-2 (and the internal review, docs/MANUSCRIPT_REVIEW_2026-07-13.md) asked whether the
base-of-Fungi amino-acid placement (Results R2) is a long-branch-attraction artifact (R2b) or real
deep signal. Rate asymmetry + the PMSF result are necessary but not sufficient to call it LBA; the
clean, direct test is outgroup removal.

Design: delete the longest, most distant branches -- the human PIF1 outgroup (1) and the
early-diverging fungi, Mucoromycota (36) + Chytridiomycota (4) -- and KEEP the nearer Basidiomycota
(235) as the outgroup that roots Ascomycota. Re-infer the AA-only tree under the SAME model as the
main analysis (LG+I+G4, 1000 UFBoot + SH-aLRT, seed 12345). Interpretation:
  * if the ScPif1 (P07271) / ScRrm3 (P38766) orthologs now collapse to a Saccharomycotina-only clade
    (as they do in the AA+3Di tree), the deep AA placement was long-branch attraction -- the clean
    demonstration the manuscript currently lacks;
  * if they stay deep (base of the remaining tree), the deep signal is NOT an outgroup-length artifact,
    and the thesis should be reworded to "structure overrides sequence" without asserting the AA
    signal is proven wrong.

Input : data/seqs/tip_map.tsv (tip_label -> group), results/seq_tree/aln.trim.fasta
        (957 x 209, the main AA alignment: MAFFT -> trimAl -automated1).
Output: results/lba/aln.noEDF.fasta (Ascomycota + Basidiomycota = 916 tips), then IQ-TREE ->
        results/lba/lba_noEDF.{treefile,iqtree,log,...}.

Usage (from repo root):
    python3 workflow/18_lba_outgroup_removal.py build   # build + verify the reduced alignment only
    python3 workflow/18_lba_outgroup_removal.py run     # build, then run IQ-TREE (long; background it)

Date: 2026-07-13.
"""
import csv
import os
import subprocess
import sys

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TIPMAP = os.path.join(REPO, "data/seqs/tip_map.tsv")
ALN = os.path.join(REPO, "results/seq_tree/aln.trim.fasta")
OUTDIR = os.path.join(REPO, "results/lba")
OUTALN = os.path.join(OUTDIR, "aln.noEDF.fasta")
PRE = os.path.join(OUTDIR, "lba_noEDF")
# arm64-native IQ-TREE 3.1.2 (Homebrew); the anaconda env's iqtree3 is x86 and needs Rosetta, which is
# unavailable on this machine. Same 3.1.2 release as the main analysis, so results are comparable.
IQTREE = "/opt/homebrew/bin/iqtree3"

KEEP_GROUPS = {"ascomycota", "basidiomycota"}  # drop: human, mucoromycota, chytridiomycota
EXPECT_KEEP = 916

mode = sys.argv[1] if len(sys.argv) > 1 else "run"
os.makedirs(OUTDIR, exist_ok=True)

# tip_label -> group; build the keep set
keep = set()
with open(TIPMAP) as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        if row["group"] in KEEP_GROUPS:
            keep.add(row["tip_label"])


def read_fasta(path):
    name, seq = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq)
                name, seq = line[1:].split()[0], []
            else:
                seq.append(line)
    if name is not None:
        yield name, "".join(seq)


n_in = n_out = 0
with open(OUTALN, "w") as out:
    for name, seq in read_fasta(ALN):
        n_in += 1
        if name in keep:
            out.write(">%s\n%s\n" % (name, seq))
            n_out += 1

sys.stderr.write("[18] alignment %d -> %d tips (dropped %d = human + Mucoromycota + Chytridiomycota)\n"
                 % (n_in, n_out, n_in - n_out))
if n_out != EXPECT_KEEP:
    sys.exit("[18] ERROR: expected %d kept tips, got %d -- check group labels" % (EXPECT_KEEP, n_out))

if mode == "build":
    sys.stderr.write("[18] build-only; reduced alignment ready at %s\n" % OUTALN)
    sys.exit(0)

# Homebrew's iqtree3 is the single-core (sequential) build, which requires -T 1.
cmd = [IQTREE, "-s", OUTALN, "-m", "LG+I+G4", "-B", "1000", "-bnni", "-alrt", "1000",
       "-T", "1", "-seed", "12345", "-pre", PRE, "-redo"]
sys.stderr.write("[18] running: %s\n" % " ".join(cmd))
subprocess.run(cmd, check=True)
sys.stderr.write("[18] done -> %s.treefile\n" % PRE)
