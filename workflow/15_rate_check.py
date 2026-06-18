#!/usr/bin/env python3
"""
15_rate_check.py  --  is RRM3 long-branched / faster-evolving than PIF1? (the LBA defense)

The AA-only tree confidently places the Pif1/Rrm3 duplication deep (base of Fungi); structure + copy
number + reconciliation say it is in the Saccharomycotina ancestor. If the AA tree is a long-branch-
attraction (LBA) / saturation artifact, we expect the Saccharomycotina copies -- especially the
"more-derived" Rrm3 -- to sit on longer branches (faster amino-acid evolution), making them prone to
being pulled toward the distant outgroups. This script tests that two ways:

  1. RATE asymmetry from branch lengths (needs a tree where Saccharomycotina is monophyletic, with
     AA branch lengths): from the duplication node (MRCA of the two anchors) -- equal divergence time
     for both daughter clades -- compare the path length to Pif1-subclade tips vs Rrm3-subclade tips.
     Also report terminal-branch lengths for Sacch-Pif1, Sacch-Rrm3, and non-Saccharomycotina tips.
     A permutation test gives a p-value for the Pif1-vs-Rrm3 difference. (numpy only, no scipy.)

  2. Tajima's (1993) relative-rate test -- model-free -- on ScPif1 vs ScRrm3 against an outgroup
     (default S. pombe Pfh1): counts sites unique to each lineage; chi-square (1 df, crit 3.84).

Run:  conda run -n pif1 python workflow/15_rate_check.py --tree results/seq_tree/aa_bl.treefile
"""
import argparse
import csv
import sys

import numpy as np
import dendropy
from Bio import SeqIO


def load_maps(tip_map, manifest):
    tip2acc, acc2tip = {}, {}
    for r in csv.DictReader(open(tip_map), delimiter="\t"):
        tip2acc[r["tip_label"]] = r["accession"]
        acc2tip[r["accession"]] = r["tip_label"]
    sub = {}
    for r in csv.DictReader(open(manifest)):
        sub[r["accession"]] = r.get("subphylum", "")
    return tip2acc, acc2tip, sub


def tajima_rrt(aln, a, b, out):
    seqs = {r.id: str(r.seq) for r in SeqIO.parse(aln, "fasta")}
    if not all(x in seqs for x in (a, b, out)):
        return None
    sa, sb, so = seqs[a], seqs[b], seqs[out]
    m1 = m2 = 0  # m1: a differs from out, b matches out; m2: b differs, a matches
    for ca, cb, co in zip(sa, sb, so):
        if "-" in (ca, cb, co) or "X" in (ca, cb, co):
            continue
        if ca == co and cb != co:
            m2 += 1
        elif cb == co and ca != co:
            m1 += 1
    chi2 = (m1 - m2) ** 2 / (m1 + m2) if (m1 + m2) else 0.0
    return m1, m2, chi2


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tree", help="tree with Saccharomycotina monophyletic + AA branch lengths")
    ap.add_argument("--aln", default="results/seq_tree/aln.trim.fasta")
    ap.add_argument("--tip-map", default="data/seqs/tip_map.tsv")
    ap.add_argument("--manifest", default="manifest.csv")
    ap.add_argument("--pif1-acc", default="P07271")
    ap.add_argument("--rrm3-acc", default="P38766")
    ap.add_argument("--outgroup-acc", default="Q9UUA2", help="Tajima outgroup (default S. pombe Pfh1)")
    args = ap.parse_args()

    tip2acc, acc2tip, sub = load_maps(args.tip_map, args.manifest)
    pif1_tip, rrm3_tip = acc2tip[args.pif1_acc], acc2tip[args.rrm3_acc]

    print("\n=== Tajima's relative-rate test (model-free) ===")
    og = acc2tip.get(args.outgroup_acc)
    res = tajima_rrt(args.aln, pif1_tip, rrm3_tip, og) if og else None
    if res:
        m1, m2, chi2 = res
        print(f"  outgroup: {og}")
        print(f"  sites unique to ScPif1: {m1}   unique to ScRrm3: {m2}")
        print(f"  chi2 = {chi2:.2f} (1 df; >3.84 => rates differ at p<0.05; >6.63 at p<0.01)")
        print(f"  => {'Rrm3 faster' if m2 > m1 else 'Pif1 faster'}; "
              f"{'SIGNIFICANT' if chi2 > 3.84 else 'not significant'}")
    else:
        print("  (outgroup or anchors missing from alignment)")

    if not args.tree:
        print("\n[15] no --tree given; skipping branch-length rate comparison "
              "(re-run with the AA-branch-length Saccharomycotina-monophyletic tree).")
        return

    print(f"\n=== branch-length rate comparison ({args.tree}) ===")
    t = dendropy.Tree.get(path=args.tree, schema="newick", preserve_underscores=True)
    hum = acc2tip.get("Q9H611")
    node = t.find_node_with_taxon_label(hum) if hum else None
    if node is not None and len(t.seed_node.child_nodes()) > 2:
        t.reroot_at_edge(node.edge, update_bipartitions=False)
    for e in t.preorder_edge_iter():
        if e.length is None:
            e.length = 0.0
    t.calc_node_root_distances()

    dup = t.mrca(taxon_labels=[pif1_tip, rrm3_tip])
    sacch_tips = {l.taxon.label for l in dup.leaf_iter()}
    kids = dup.child_nodes()
    pif_clade = next(c for c in kids if pif1_tip in {l.taxon.label for l in c.leaf_iter()})
    rrm_clade = next(c for c in kids if rrm3_tip in {l.taxon.label for l in c.leaf_iter()})

    def node_to_tip(clade):
        return np.array([l.root_distance - dup.root_distance for l in clade.leaf_iter()])

    dP, dR = node_to_tip(pif_clade), node_to_tip(rrm_clade)
    print(f"  Saccharomycotina clade under duplication node: {len(sacch_tips)} tips "
          f"(Pif1 subclade {len(dP)}, Rrm3 subclade {len(dR)})")
    print(f"  mean root(dup)->tip   Pif1 {dP.mean():.3f}   Rrm3 {dR.mean():.3f}   "
          f"ratio Rrm3/Pif1 = {dR.mean()/dP.mean():.2f}")
    print(f"  median                Pif1 {np.median(dP):.3f}   Rrm3 {np.median(dR):.3f}")

    # permutation test on the difference of means
    obs = dR.mean() - dP.mean()
    pool = np.concatenate([dP, dR]); n = len(dP)
    rng = np.random.default_rng(42)
    perm = np.array([(lambda s: s[:n].mean() - s[n:].mean())(rng.permutation(pool)) for _ in range(10000)])
    p = (np.abs(perm) >= abs(obs)).mean()
    print(f"  permutation p (Pif1 vs Rrm3 rate difference): {p:.4f} "
          f"({'significant' if p < 0.05 else 'n.s.'})")

    # terminal branch lengths by group
    def tip_edges(labels):
        return np.array([l.edge.length for l in t.leaf_node_iter()
                         if l.taxon.label in labels and l.edge.length is not None])
    sp = {x for x in sacch_tips if x == pif1_tip or x in {l.taxon.label for l in pif_clade.leaf_iter()}}
    sr = {l.taxon.label for l in rrm_clade.leaf_iter()}
    nonsacch = {l.taxon.label for l in t.leaf_node_iter()} - sacch_tips
    print(f"  mean terminal branch  Pif1 {tip_edges(sp).mean():.3f}   "
          f"Rrm3 {tip_edges(sr).mean():.3f}   non-Sacch {tip_edges(nonsacch).mean():.3f}")
    print("  (Saccharomycotina copies longer than non-Sacch => post-duplication acceleration / "
          "LBA-prone, supporting that the AA deep placement is artifactual)")


if __name__ == "__main__":
    main()
