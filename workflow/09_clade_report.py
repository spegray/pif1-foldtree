#!/usr/bin/env python3
"""
09_clade_report.py  --  preliminary read on WHERE the PIF1/RRM3 duplication falls.

NOT the final answer (that comes from Stage 6 reconciliation against the species tree). This is
a fast sanity look at a single gene tree: find the MRCA of the two S. cerevisiae anchors
Pif1 (P07271) and Rrm3 (P38766) -- in a gene tree, the MRCA of two same-species paralogs IS
their duplication node -- split it into its two child subtrees (the PIF1 lineage and the RRM3
lineage), and tabulate the taxonomic composition of each. If the RRM3 lineage is Saccharomycotina-
only, the duplication is Saccharomycotina-restricted (the README hypothesis).

Usage:
    conda run -n pif1 python workflow/09_clade_report.py --tree results/seq_tree/pif1.treefile
    conda run -n pif1 python workflow/09_clade_report.py --tree results/struct_tree/fident_foldtree.rooted.nwk
"""
import argparse
import csv
import sys
from collections import Counter

import dendropy


def load_maps(tip_map, manifest):
    acc2tip, tip2acc = {}, {}
    for r in csv.DictReader(open(tip_map), delimiter="\t"):
        acc2tip[r["accession"]] = r["tip_label"]
        tip2acc[r["tip_label"]] = r["accession"]
    sub, grp, org = {}, {}, {}
    for r in csv.DictReader(open(manifest)):
        sub[r["accession"]] = r.get("subphylum", "") or "(none)"
    # coarse group from tip_map (ascomycota/basidiomycota/mucoromycota/chytridiomycota/human)
    for r in csv.DictReader(open(tip_map), delimiter="\t"):
        grp[r["tip_label"]] = r.get("group", "") or "(none)"
        org[r["tip_label"]] = r.get("organism", "")
    return acc2tip, tip2acc, sub, grp, org


def describe(node, tip2acc, sub, grp):
    tips = [l.taxon.label for l in node.leaf_iter()]
    subc = Counter(sub.get(tip2acc.get(t, ""), "(none)") for t in tips)
    grpc = Counter(grp.get(t, "(none)") for t in tips)
    return tips, subc, grpc


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tree", required=True)
    ap.add_argument("--tip-map", default="data/seqs/tip_map.tsv")
    ap.add_argument("--manifest", default="manifest.csv")
    ap.add_argument("--pif1-acc", default="P07271")
    ap.add_argument("--rrm3-acc", default="P38766")
    ap.add_argument("--root-on", default="Q9H611", help="accession to outgroup-root on if tree is unrooted")
    args = ap.parse_args()

    acc2tip, tip2acc, sub, grp, org = load_maps(args.tip_map, args.manifest)
    pif1_tip, rrm3_tip = acc2tip[args.pif1_acc], acc2tip[args.rrm3_acc]

    t = dendropy.Tree.get(path=args.tree, schema="newick", preserve_underscores=True)

    # root on the human outgroup if present and the tree looks unrooted (>2 children at seed)
    if args.root_on in acc2tip:
        og = t.find_node_with_taxon_label(acc2tip[args.root_on])
        if og is not None and len(t.seed_node.child_nodes()) > 2:
            t.reroot_at_edge(og.edge, update_bipartitions=False)
            sys.stderr.write(f"[09] rooted on {acc2tip[args.root_on]} ({args.root_on})\n")

    mrca = t.mrca(taxon_labels=[pif1_tip, rrm3_tip])
    print(f"\n=== {args.tree} ===")
    print(f"anchors: PIF1={pif1_tip}  RRM3={rrm3_tip}")
    print(f"duplication node = MRCA(Pif1,Rrm3): subtends {len(list(mrca.leaf_iter()))} leaves "
          f"of {len(t.leaf_nodes())} total")

    kids = mrca.child_nodes()
    if len(kids) < 2:
        print("  (MRCA has <2 children -- unresolved; inspect manually)")
        return

    for side, anchor in (("PIF1", pif1_tip), ("RRM3", rrm3_tip)):
        # the child subtree containing this anchor
        child = next((c for c in kids if anchor in {l.taxon.label for l in c.leaf_iter()}), None)
        if child is None:
            print(f"  [{side}] anchor not under either child?! (gene-tree messiness)")
            continue
        tips, subc, grpc = describe(child, tip2acc, sub, grp)
        nsp = len({org[x] for x in tips if x in org})
        print(f"\n  --- {side} lineage (the child holding {anchor}) ---")
        print(f"      {len(tips)} genes / ~{nsp} species")
        print(f"      by broad group : " + ", ".join(f"{k}={v}" for k, v in grpc.most_common()))
        print(f"      by subphylum   : " + ", ".join(f"{k}={v}" for k, v in subc.most_common(8)))
        if side == "RRM3":
            non_sacch = {k: v for k, v in subc.items()
                         if "saccharomycotina" not in k.lower() and k not in ("(none)",)}
            verdict = ("Saccharomycotina-RESTRICTED (hypothesis supported on this tree)"
                       if not non_sacch else
                       f"NOT Saccharomycotina-only -- also: {non_sacch}")
            print(f"      >> RRM3-lineage verdict: {verdict}")


if __name__ == "__main__":
    main()
