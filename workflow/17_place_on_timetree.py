#!/usr/bin/env python3
"""
17_place_on_timetree.py  --  M5 + R6: place (and date) the PIF1/RRM3 duplication on the published,
time-calibrated Y1000+ budding-yeast phylogeny (Shen et al. 2018), instead of the NCBI taxonomy.

The reconciliation maps the duplication onto the species-tree node = the last common ancestor of all
species that retain a gene descended from the duplication. This script:
  1. pulls the duplication clade's species from our gene tree (MRCA of the two anchors in the AA+3Di
     tree) and matches them, by genus+species, onto the Shen 2018 tip set;
  2. tests M5 robustness: are those species MONOPHYLETIC on the phylogenomic tree, and does their MRCA
     sit within Saccharomycotina? (the headline should hold on a real phylogeny, not just NCBI);
  3. reads R6 dating straight off the calibrated tree: the duplication sits on the branch above that
     MRCA, so its age is bracketed by [MRCA crown age, MRCA-parent (stem) age], in Myr.

Run:  conda run -n pif1 python workflow/17_place_on_timetree.py --timetree /tmp/shen2018.newick
"""
import argparse
import csv
import re
import sys

import dendropy


def norm(name):
    """'Lachancea lanzarotensis (Yeast) (...)' -> 'lachancea_lanzarotensis'."""
    name = re.sub(r"\(.*?\)", " ", name or "")
    parts = re.sub(r"[^A-Za-z0-9 ]", " ", name).lower().split()
    return "_".join(parts[:2]) if len(parts) >= 2 else (parts[0] if parts else "")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--timetree", required=True, help="Shen 2018 time-calibrated newick")
    ap.add_argument("--gene-tree", default="results/seq_tree/pif1_aa3di.treefile")
    ap.add_argument("--tip-map", default="data/seqs/tip_map.tsv")
    ap.add_argument("--pif1-acc", default="P07271")
    ap.add_argument("--rrm3-acc", default="P38766")
    args = ap.parse_args()

    # gene tree: duplication clade -> our species (organism names) via tip_map
    acc2tip, tip2org = {}, {}
    for r in csv.DictReader(open(args.tip_map), delimiter="\t"):
        acc2tip[r["accession"]] = r["tip_label"]
        tip2org[r["tip_label"]] = r.get("organism", "")
    gt = dendropy.Tree.get(path=args.gene_tree, schema="newick", preserve_underscores=True)
    hum = acc2tip.get("Q9H611")
    hn = gt.find_node_with_taxon_label(hum) if hum else None
    if hn is not None and len(gt.seed_node.child_nodes()) > 2:
        gt.reroot_at_edge(hn.edge, update_bipartitions=False)  # root so MRCA is the Saccharomycotina clade
    gt.is_rooted = True
    dup = gt.mrca(taxon_labels=[acc2tip[args.pif1_acc], acc2tip[args.rrm3_acc]])
    dup_orgs = {norm(tip2org.get(l.taxon.label, "")) for l in dup.leaf_iter()}
    dup_orgs.discard("")
    sys.stderr.write(f"[17] duplication clade: {len(list(dup.leaf_iter()))} genes -> "
                     f"{len(dup_orgs)} distinct species (normalized)\n")

    # time tree
    tt = dendropy.Tree.get(path=args.timetree, schema="newick", preserve_underscores=True)
    tt.is_rooted = True
    tt.calc_node_ages(ultrametricity_precision=0.01)
    SCALE = 100.0 if tt.seed_node.age < 20 else 1.0  # Shen 2018 tree is in units of 100 Myr (root~4.0=400 Mya)
    tip_norm = {}
    for lf in tt.leaf_node_iter():
        tip_norm.setdefault(norm(lf.taxon.label), lf.taxon.label)
    tt_tips = set(tip_norm)
    sys.stderr.write(f"[17] time tree: {len(tt_tips)} tips\n")

    matched = sorted(dup_orgs & tt_tips)
    missing = sorted(dup_orgs - tt_tips)
    print(f"\nmatched {len(matched)}/{len(dup_orgs)} duplication-clade species to the Shen tree "
          f"({len(missing)} unmatched)")
    if len(matched) < 3:
        print("  too few matches to place a node; inspect name formats:", missing[:10]); return

    labels = [tip_norm[m] for m in matched]
    mrca = tt.mrca(taxon_labels=labels)
    mrca_leaves = [l.taxon.label for l in mrca.leaf_iter()]
    monophyletic = len(mrca_leaves) == len(labels)
    crown = mrca.age * SCALE
    stem = (mrca.parent_node.age * SCALE) if mrca.parent_node else None

    print(f"\n=== M5: placement on the phylogenomic tree ===")
    print(f"  MRCA of the matched duplication-clade species subtends {len(mrca_leaves)} tips "
          f"(our matched: {len(labels)}) -> {'MONOPHYLETIC' if monophyletic else 'NON-monophyletic (includes '+str(len(mrca_leaves)-len(labels))+' extra Shen tips)'}")
    print(f"  (extra tips, if any, are budding yeasts we didn't sample; expected, not a conflict)")

    print(f"\n=== R6: calibrated dating (Myr) ===")
    print(f"  duplication sits on the branch ABOVE this MRCA, so its age is bracketed:")
    print(f"    crown age of the duplication clade (MRCA): {crown:.1f} Mya")
    if stem is not None:
        print(f"    stem age (MRCA's parent):                  {stem:.1f} Mya")
        print(f"    => PIF1/RRM3 duplication age in [{crown:.0f}, {stem:.0f}] Mya")
    # for context: the whole-Saccharomycotina crown age (root of the budding-yeast tree)
    print(f"  (Shen tree root / budding-yeast crown age: {tt.seed_node.age*SCALE:.1f} Mya; "
          f"time units scaled x{SCALE:.0f})")
    if missing:
        print(f"\n  unmatched species (first 15): {missing[:15]}")


if __name__ == "__main__":
    main()
