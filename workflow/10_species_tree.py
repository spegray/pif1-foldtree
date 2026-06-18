#!/usr/bin/env python3
"""
10_species_tree.py  --  Stage 5 (fast pass): NCBI-taxonomy species tree for our 728 species.

Reconciliation (Stage 6) maps the PIF1/RRM3 duplication onto a branch of a trusted SPECIES tree.
This builds the quick first-pass backbone straight from NCBI taxonomy (via ete3 NCBITaxa) using
the taxids already in tip_map.tsv -- no manual tree wrangling. A published phylogenomic tree
(Y1000+/Li et al.) will replace it later as the rigorous confirmation; the rest of Stage 6 is
identical, so swapping backbones is a one-file change.

Outputs
-------
  data/species_tree/ncbi_species.nwk   rooted, binary species tree; leaves = taxid strings
  data/species_tree/gene_species.map   GeneRax mapping, one line per gene:  <tip_label> <taxid>

Binary: GeneRax's DL models want a binary species tree, so NCBI polytomies (e.g. many genera in a
family) are arbitrarily resolved -- harmless for placing a DEEP duplication node, and the published
backbone fixes the fine structure later.

First run downloads the NCBI taxdump (~once, a few minutes) into ~/.etetoolkit/taxa.sqlite.

Run in the pif1 env (needs ete3):
    conda run -n pif1 python workflow/10_species_tree.py
"""
import argparse
import csv
import os
import sys


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tip-map", default="data/seqs/tip_map.tsv")
    ap.add_argument("--out-tree", default="data/species_tree/ncbi_species.nwk")
    ap.add_argument("--out-map", default="data/species_tree/gene_species.map")
    ap.add_argument("--update-db", action="store_true", help="force re-download of the NCBI taxdump")
    args = ap.parse_args()

    from ete3 import NCBITaxa
    ncbi = NCBITaxa()
    if args.update_db:
        ncbi.update_taxonomy_database()

    rows = list(csv.DictReader(open(args.tip_map), delimiter="\t"))
    taxids = sorted({int(r["taxid"]) for r in rows if r.get("taxid")})
    sys.stderr.write(f"[10] {len(rows)} genes across {len(taxids)} taxids\n")

    # NCBI may have merged/renamed some taxids; translate and note any it can't place
    valid = set(ncbi.get_taxid_translator(taxids).keys())
    missing = [t for t in taxids if t not in valid]
    if missing:
        sys.stderr.write(f"[10] WARN: {len(missing)} taxids not found in NCBI taxonomy "
                         f"(genes on these will be dropped from reconciliation): {missing[:10]}\n")

    tree = ncbi.get_topology(sorted(valid), intermediate_nodes=False)
    # GeneRax wants a rooted BINARY tree. Two cleanups on the raw NCBI topology:
    #  (1) suppress unifurcations -- single-child internal nodes appear when a rank (e.g. a class)
    #      has only one sampled representative; node.delete() reconnects the child to the grandparent.
    #  (2) resolve remaining polytomies (many genera in a family) into arbitrary bifurcations.
    changed = True
    while changed:
        changed = False
        for node in tree.traverse():
            if not node.is_leaf() and len(node.children) == 1:
                node.delete()
                changed = True
                break
    while len(tree.children) == 1:          # collapse a unifurcating root
        tree = tree.children[0]
        tree.up = None
    tree.resolve_polytomy(recursive=True)
    leaf_taxids = {int(l.name) for l in tree.get_leaves()}
    sys.stderr.write(f"[10] species tree: {len(leaf_taxids)} leaves "
                     f"({len(tree.get_descendants())} nodes after binary resolution)\n")

    os.makedirs(os.path.dirname(args.out_tree), exist_ok=True)
    tree.write(outfile=args.out_tree, format=9)   # leaf names only; topology is what DL needs

    # gene -> species map, dropping genes whose species isn't in the tree
    n_map = n_drop = 0
    with open(args.out_map, "w") as out:
        for r in rows:
            tx = int(r["taxid"]) if r.get("taxid") else None
            if tx in leaf_taxids:
                out.write(f"{r['tip_label']} {tx}\n")
                n_map += 1
            else:
                n_drop += 1

    sys.stderr.write(f"[10] wrote species tree -> {args.out_tree}\n")
    sys.stderr.write(f"[10] wrote {n_map} gene->species mappings -> {args.out_map}"
                     + (f" ({n_drop} genes dropped)\n" if n_drop else "\n"))
    sys.stderr.write("[10] NOTE: this is the NCBI-taxonomy FAST pass; swap in the published "
                     "phylogenomic backbone later for the rigorous answer.\n")


if __name__ == "__main__":
    main()
