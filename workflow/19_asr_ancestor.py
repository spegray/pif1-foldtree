#!/usr/bin/env python3
"""
19_asr_ancestor.py -- reconstruct the pre-duplication Saccharomycotina PIF1 helicase core.

Marginal (empirical-Bayes) ancestral sequence reconstruction at the DUPLICATION NODE -- the MRCA of the
Pif1 clade and the Rrm3 clade (i.e. the single ancestral gene immediately before it split into Pif1 and
Rrm3). We fix the structure-informed AA+3Di topology (-te) and reconstruct ancestral AMINO ACIDS under the
same LG+I+G4 model as the main sequence analysis; IQ-TREE writes per-site posterior probabilities for
every internal node to results/asr/asr.state.

Scope + caveat: this reconstructs the trimmed HELICASE CORE (209 columns) only -- the accessory N/C-terminal
domains were trimmed away (they cannot be aligned family-wide) and are not reconstructed here. Confidence is
site-dependent: high at conserved motifs (incl. the Arg-wedge position), low at fast/saturated loops. See
results/asr/extract_ancestor.R for the dup-node extraction, the per-site confidence profile, and the wedge check.

Input : results/seq_tree/aln.trim.fasta (957 x 209 AA core), results/seq_tree/pif1_aa3di.treefile (topology).
Output: results/asr/asr.{state,iqtree,log,treefile}.
Run from repo root:  python3 workflow/19_asr_ancestor.py    Date: 2026-07-13.
"""
import os
import subprocess
import sys

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ALN = os.path.join(REPO, "results/seq_tree/aln.trim.fasta")
TREE = os.path.join(REPO, "results/seq_tree/pif1_aa3di.treefile")
OUTDIR = os.path.join(REPO, "results/asr")
PRE = os.path.join(OUTDIR, "asr")
IQTREE = "/opt/homebrew/bin/iqtree3"  # arm64-native 3.1.2 (single-core build -> -T 1)

os.makedirs(OUTDIR, exist_ok=True)
cmd = [IQTREE, "-s", ALN, "-te", TREE, "-m", "LG+I+G4", "-asr",
       "-T", "1", "-seed", "12345", "-pre", PRE, "-redo"]
sys.stderr.write("[19] running: %s\n" % " ".join(cmd))
subprocess.run(cmd, check=True)
sys.stderr.write("[19] done -> %s.state\n" % PRE)
