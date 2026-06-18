#!/usr/bin/env python3
"""
18_filter_alignments.py  --  build robustness-rerun inputs for the structural signal (review M6).

Two worries about the AA+3Di result: (1) 3Di from low-confidence AlphaFold regions could be noise;
(2) the 828 AFDB vs 129 ColabFold predictor split could induce a batch artifact. We test both by
re-running the AA+3Di tree on filtered subsets of the SAME alignment and checking the duplication
still maps to Saccharomycotina. These subsets are built here (on the Mac, where the structures live)
so the Windows box only has to run IQ-TREE.

Outputs (results/seq_tree/):
  aln_pLDDT.fasta     alnAA_3di.fasta keeping only tips whose core mean Cα pLDDT >= threshold
  aln_AFDBonly.fasta  alnAA_3di.fasta dropping the 129 ColabFold structures (AFDB predictor only)
The partition file (aa3di_partition.nex) is reused unchanged for both (columns are identical).

Run:  conda run -n pif1 python workflow/18_filter_alignments.py --plddt 70
"""
import argparse
import csv
import os
import sys

import numpy as np
from Bio import SeqIO


def core_mean_plddt(pdb):
    vals = []
    with open(pdb) as fh:
        for line in fh:
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                try:
                    vals.append(float(line[60:66]))
                except ValueError:
                    pass
    return float(np.mean(vals)) if vals else None


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--aln", default="results/seq_tree/alnAA_3di.fasta")
    ap.add_argument("--cores", default="data/structures/cores")
    ap.add_argument("--tip-map", default="data/seqs/tip_map.tsv")
    ap.add_argument("--manifest", default="manifest.csv")
    ap.add_argument("--plddt", type=float, default=70.0)
    args = ap.parse_args()

    tip2acc = {}
    for r in csv.DictReader(open(args.tip_map), delimiter="\t"):
        tip2acc[r["tip_label"]] = r["accession"]
    src = {}
    for r in csv.DictReader(open(args.manifest)):
        src[r["accession"]] = r.get("structure_source", "")

    recs = list(SeqIO.parse(args.aln, "fasta"))
    plddt = {}
    for rec in recs:
        pdb = os.path.join(args.cores, f"{rec.id}.pdb")
        plddt[rec.id] = core_mean_plddt(pdb) if os.path.exists(pdb) else None

    vals = np.array([v for v in plddt.values() if v is not None])
    pct = np.percentile(vals, [5, 25, 50, 75, 95])
    sys.stderr.write(f"[18] core mean pLDDT over {len(vals)} tips: "
                     f"5/25/50/75/95% = {pct[0]:.0f}/{pct[1]:.0f}/{pct[2]:.0f}/{pct[3]:.0f}/{pct[4]:.0f}\n")

    keep_p = [r for r in recs if plddt.get(r.id) is not None and plddt[r.id] >= args.plddt]
    keep_a = [r for r in recs if src.get(tip2acc.get(r.id, ""), "") == "AFDB"]

    with open("results/seq_tree/aln_pLDDT.fasta", "w") as fh:
        for r in keep_p:
            fh.write(f">{r.id}\n{r.seq}\n")
    with open("results/seq_tree/aln_AFDBonly.fasta", "w") as fh:
        for r in keep_a:
            fh.write(f">{r.id}\n{r.seq}\n")

    # sanity: the two anchors must survive both filters (else the test can't place the node)
    anchors = {"P07271": None, "P38766": None}
    acc2tip = {v: k for k, v in tip2acc.items()}
    for acc in anchors:
        t = acc2tip.get(acc)
        anchors[acc] = (t, plddt.get(t), src.get(acc, ""))
    sys.stderr.write(f"[18] anchors: " + "; ".join(
        f"{acc}={t} (pLDDT {p:.0f}, {s})" for acc, (t, p, s) in anchors.items() if t) + "\n")
    sys.stderr.write(f"[18] wrote aln_pLDDT.fasta ({len(keep_p)}/{len(recs)} tips, >= {args.plddt:.0f}) "
                     f"and aln_AFDBonly.fasta ({len(keep_a)} tips, AFDB only)\n")
    sys.stderr.write("[18] reuse results/seq_tree/aa3di_partition.nex for both (columns unchanged).\n")


if __name__ == "__main__":
    main()
