#!/usr/bin/env python3
"""
11_reconcile_lca.py  --  Stage 6 (fast pass): LCA reconciliation -> which clade the duplication maps to.

Transparent gene-tree/species-tree reconciliation by the LCA (last-common-ancestor) rule -- the core
of what NOTUNG/GeneRax do, but dependency-light and fully auditable. For the gene-tree node that is
the MRCA of the two S. cerevisiae anchors (Pif1 P07271, Rrm3 P38766) -- i.e. the PIF1/RRM3
duplication node -- we collect the species of every gene descended from it and ask: what is the
NCBI clade that is their last common ancestor? THAT clade's stem is the branch the duplication maps
onto. We also report the species/subphylum composition so "Saccharomycotina-restricted" can be read
off directly, and the sister context (stem vs crown).

This reproduces -- and quantifies -- the sequence-vs-structure tension: run it on both gene trees.
GeneRax (ML reconciliation, accounts for gene-tree error) is the rigorous confirmation (next step).

CAVEAT: LCA reconciliation takes the gene tree at face value, so it is only as good as the tree
(one misplaced deep tip drags the LCA deeper). That is exactly why we cross-check sequence vs
structure here, and follow with GeneRax-SPR (which can correct weakly-supported branches).

Usage:
    conda run -n pif1 python workflow/11_reconcile_lca.py --tree results/seq_tree/pif1.treefile
    conda run -n pif1 python workflow/11_reconcile_lca.py --tree results/struct_tree/fident_foldtree.rooted.nwk
"""
import argparse
import csv
import sys
from collections import Counter

import dendropy
from ete3 import NCBITaxa


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tree", required=True)
    ap.add_argument("--tip-map", default="data/seqs/tip_map.tsv")
    ap.add_argument("--manifest", default="manifest.csv")
    ap.add_argument("--pif1-acc", default="P07271")
    ap.add_argument("--rrm3-acc", default="P38766")
    ap.add_argument("--root-on", default="Q9H611")
    args = ap.parse_args()

    ncbi = NCBITaxa()
    tip2tax, acc2tip, sub = {}, {}, {}
    for r in csv.DictReader(open(args.tip_map), delimiter="\t"):
        tip2tax[r["tip_label"]] = int(r["taxid"]) if r.get("taxid") else None
        acc2tip[r["accession"]] = r["tip_label"]
    for r in csv.DictReader(open(args.manifest)):
        sub[r["accession"]] = r.get("subphylum", "") or "(none)"
    tip2acc = {v: k for k, v in acc2tip.items()}

    pif1_tip, rrm3_tip = acc2tip[args.pif1_acc], acc2tip[args.rrm3_acc]
    t = dendropy.Tree.get(path=args.tree, schema="newick", preserve_underscores=True)

    # root on human outgroup if the tree is unrooted
    if args.root_on in acc2tip:
        og = t.find_node_with_taxon_label(acc2tip[args.root_on])
        if og is not None and len(t.seed_node.child_nodes()) > 2:
            t.reroot_at_edge(og.edge, update_bipartitions=False)

    dup = t.mrca(taxon_labels=[pif1_tip, rrm3_tip])
    leaves = [l.taxon.label for l in dup.leaf_iter()]
    taxids = [tip2tax[x] for x in leaves if tip2tax.get(x)]

    # NCBI last-common-ancestor clade of all species under the duplication node
    lineages = {tx: ncbi.get_lineage(tx) for tx in set(taxids) if ncbi.get_lineage(tx)}
    common = set.intersection(*[set(v) for v in lineages.values()]) if lineages else set()
    ref = next(iter(lineages.values()))
    mrca_tx = max(common, key=lambda x: ref.index(x)) if common else None
    name = ncbi.get_taxid_translator([mrca_tx])[mrca_tx] if mrca_tx else "?"
    rank = ncbi.get_rank([mrca_tx]).get(mrca_tx, "?") if mrca_tx else "?"

    subc = Counter(sub.get(tip2acc.get(x, ""), "(none)") for x in leaves)
    nsp = len(set(taxids))

    print(f"\n=== {args.tree} ===")
    print(f"duplication node = MRCA(Pif1,Rrm3): {len(leaves)} genes / {nsp} species")
    print(f">> the PIF1/RRM3 duplication maps to the stem of:  **{name}**  (NCBI rank: {rank})")
    print(f"   subphylum composition under the node: "
          + ", ".join(f"{k}={v}" for k, v in subc.most_common(8)))

    # sister context: what is immediately OUTSIDE the duplication node (-> stem vs crown read)
    if dup.parent_node is not None:
        sibs = [c for c in dup.parent_node.child_nodes() if c is not dup]
        sib_tax = [tip2tax[l.taxon.label] for s in sibs for l in s.leaf_iter()
                   if tip2tax.get(l.taxon.label)]
        sib_lin = {tx: ncbi.get_lineage(tx) for tx in set(sib_tax) if ncbi.get_lineage(tx)}
        if sib_lin:
            sc = set.intersection(*[set(v) for v in sib_lin.values()])
            sref = next(iter(sib_lin.values()))
            sx = max(sc, key=lambda x: sref.index(x)) if sc else None
            sname = ncbi.get_taxid_translator([sx])[sx] if sx else "?"
            print(f"   immediate sister group is within: {sname} "
                  f"({len(set(sib_tax))} species)  -> read stem-vs-crown from this")


if __name__ == "__main__":
    main()
