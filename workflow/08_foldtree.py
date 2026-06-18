#!/usr/bin/env python3
"""
08_foldtree.py  --  Stage 4b: structural (FoldTree) phylogeny from the helicase-core structures.

Reproduces the DessimozLab FoldTree method (Moi et al., Nat Struct Mol Biol 2025) faithfully,
using only tools already in the pif1 env (foldseek, fastme, dendropy, numpy, pandas) -- no
Snakemake, no Colab. Pipeline:

  1. foldseek easy-search, all-vs-all, exhaustive, 3Di+AA (alignment-type 2 = default).
  2. For each similarity metric (fident primary; alntmscore, lddt as robustness checks):
       - symmetric average of the two directional scores,
       - similarity -> distance via the FoldTree correction  d = -b * ln(1 - (1-sim)/b)
         (this is exactly their Tajima_dist series; b=0.93 for fident, 0.95 otherwise),
       - FastME (NNI) on the distance matrix -> tree,
       - negative branch lengths clamped to a small delta,
       - outgroup-rooted on the human PIF1 tip (Q9H611) for orientation.

Tips are named by the SAME tip_label as the IQ-TREE sequence tree (07 named the core .pdb that
way), so the two trees are directly comparable -> the sequence-vs-structure concordance QC.

Run in the pif1 env:
    conda run -n pif1 python workflow/08_foldtree.py \
        --foldseek $CONDA_PREFIX/bin/foldseek --fastme $CONDA_PREFIX/bin/fastme
"""
import argparse
import os
import subprocess
import sys

import numpy as np
import pandas as pd
import dendropy

COLS = ("query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,"
        "evalue,bits,lddt,lddtfull,alntmscore").split(",")
METRICS = {"fident": 0.93, "alntmscore": 0.95, "lddt": 0.95}   # metric -> b (FoldTree bfactor)


def run(cmd):
    sys.stderr.write("[08] $ " + " ".join(cmd) + "\n")
    subprocess.run(cmd, check=True)


def run_foldseek(foldseek, coresdir, out_m8, tmp, max_seqs, evalue):
    run([foldseek, "easy-search", coresdir, coresdir, out_m8, tmp,
         "--format-output", ",".join(COLS),
         "--exhaustive-search", "-e", str(evalue), "--max-seqs", str(max_seqs)])
    return out_m8


def foldtree_distance(sim, b):
    """FoldTree correction: d = -b * ln(1 - (1-sim)/b), clamped to stay finite."""
    raw = 1.0 - sim                                   # raw distance, 0 (identical) .. 1
    x = np.clip(raw / b, 0.0, 1.0 - 1e-9)             # keep the log argument in (0,1]
    d = -b * np.log(1.0 - x)
    d[d <= 0] = 0.0                                    # kill IEEE -0.0 / tiny negs (FastME rejects '-0.0000')
    np.fill_diagonal(d, 0.0)
    return d


def write_phylip(ids, dm, path):
    with open(path, "w") as fh:
        fh.write(f"{len(ids)}\n")
        for i, name in enumerate(ids):
            fh.write(name + " " + " ".join(f"{v:.4f}" for v in dm[i]) + "\n")
    return path


def clamp_negatives(in_nwk, out_nwk, delta=1e-4):
    t = dendropy.Tree.get(path=in_nwk, schema="newick", preserve_underscores=True)
    for e in t.preorder_edge_iter():
        if e.length is not None and e.length < 0:
            e.length = delta
    t.write(path=out_nwk, schema="newick", unquoted_underscores=True)
    return out_nwk


def outgroup_root(in_nwk, out_nwk, og_label):
    t = dendropy.Tree.get(path=in_nwk, schema="newick", preserve_underscores=True)
    node = t.find_node_with_taxon_label(og_label)
    if node is None:
        sys.stderr.write(f"[08] WARN: outgroup {og_label} not found; leaving unrooted\n")
        return None
    t.reroot_at_edge(node.edge, update_bipartitions=False)
    t.write(path=out_nwk, schema="newick", unquoted_underscores=True)
    return out_nwk


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cores", default="data/structures/cores")
    ap.add_argument("--outdir", default="results/struct_tree")
    ap.add_argument("--tip-map", default="data/seqs/tip_map.tsv")
    ap.add_argument("--foldseek", default="foldseek")
    ap.add_argument("--fastme", default="fastme")
    ap.add_argument("--outgroup-acc", default="Q9H611", help="root on this accession's tip (human PIF1)")
    ap.add_argument("--max-seqs", type=int, default=2000)
    ap.add_argument("--evalue", type=float, default=100.0)
    ap.add_argument("--redo", action="store_true", help="re-run foldseek even if res.m8 exists")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    out_m8 = os.path.join(args.outdir, "res.m8")
    tmp = os.path.join(args.outdir, "tmp")

    if args.redo or not os.path.exists(out_m8):
        run_foldseek(args.foldseek, args.cores, out_m8, tmp, args.max_seqs, args.evalue)
    else:
        sys.stderr.write(f"[08] reusing existing {out_m8}\n")

    res = pd.read_table(out_m8, header=None, names=COLS,
                        usecols=["query", "target", *METRICS.keys()])  # skip the huge lddtfull col
    res["query"] = res["query"].str.replace(".pdb", "", regex=False)
    res["target"] = res["target"].str.replace(".pdb", "", regex=False)
    ids = sorted(set(res["query"]) | set(res["target"]))
    pos = {p: i for i, p in enumerate(ids)}
    sys.stderr.write(f"[08] {len(res)} foldseek rows over {len(ids)} structures\n")

    qi = res["query"].map(pos).to_numpy()
    ti = res["target"].map(pos).to_numpy()

    # outgroup tip_label from accession
    og_label = None
    for line in open(args.tip_map):
        f = line.rstrip("\n").split("\t")
        if len(f) > 1 and f[1] == args.outgroup_acc:
            og_label = f[0]
            break

    n = len(ids)
    for metric, b in METRICS.items():
        S = np.zeros((n, n))
        vals = res[metric].to_numpy(dtype=float)
        np.add.at(S, (qi, ti), vals)          # both directions reported by foldseek
        np.add.at(S, (ti, qi), vals)          # symmetrize
        S /= 2.0
        dm = foldtree_distance(S, b)
        phy = write_phylip(ids, dm, os.path.join(args.outdir, f"{metric}_distmat.phy"))
        raw_nwk = os.path.join(args.outdir, f"{metric}_fastme.nwk")
        run([args.fastme, "-i", phy, "-o", raw_nwk, "-n"])
        clean = clamp_negatives(raw_nwk, os.path.join(args.outdir, f"{metric}_foldtree.nwk"))
        if og_label:
            outgroup_root(clean, os.path.join(args.outdir, f"{metric}_foldtree.rooted.nwk"), og_label)
        sys.stderr.write(f"[08] {metric}: tree -> {clean}\n")

    sys.stderr.write(f"\n[08] done. Primary structural tree: "
                     f"{args.outdir}/fident_foldtree.rooted.nwk\n")
    sys.stderr.write("[08] (alntmscore / lddt trees written too -> metric-robustness check)\n")


if __name__ == "__main__":
    main()
