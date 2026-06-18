#!/usr/bin/env python3
"""
12_reconcile.py  --  Stage 6: species-overlap reconciliation -> the duplication landscape.

GeneRax is unusable here (segfaults under Rosetta), so we reconcile with the species-overlap
algorithm (Huerta-Cepas et al.), the standard, citable method built into ete3 -- it classifies
every internal gene-tree node as a Duplication or Speciation from whether its two daughter clades
share any species, independent of branch lengths. For each DUPLICATION we report the NCBI clade its
species set maps to. Instead of trusting the single MRCA(Pif1,Rrm3) node (which a few rogue tips can
drag deep), this gives the whole landscape: a real, ancestral PIF1/RRM3 duplication should show up as
duplications shared across MANY species all mapping to the SAME clade.

Reports, per tree:
  - total duplication / speciation events,
  - duplications mapped to each NCBI clade (sorted), restricted to those shared by >= MIN_SP species
    (a duplication seen in 1 species is a recent/in-paralog event, not the ancestral one we want),
  - the clade the ScPif1/ScRrm3 anchor node itself maps to.

Run:  conda run -n pif1 python workflow/12_reconcile.py --tree results/seq_tree/pif1.treefile
      conda run -n pif1 python workflow/12_reconcile.py --tree results/struct_tree/fident_foldtree.rooted.nwk
"""
import argparse
import csv
import sys
from collections import Counter

from ete3 import PhyloTree, NCBITaxa


def ncbi_lca(ncbi, taxids):
    lin = {t: ncbi.get_lineage(t) for t in set(taxids) if t and ncbi.get_lineage(t)}
    if not lin:
        return "?", "?", 0
    common = set.intersection(*[set(v) for v in lin.values()])
    ref = next(iter(lin.values()))
    mrca = max(common, key=lambda x: ref.index(x)) if common else None
    name = ncbi.get_taxid_translator([mrca])[mrca] if mrca else "?"
    rank = ncbi.get_rank([mrca]).get(mrca, "?") if mrca else "?"
    return name, rank, len(lin)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tree", required=True)
    ap.add_argument("--tip-map", default="data/seqs/tip_map.tsv")
    ap.add_argument("--pif1-acc", default="P07271")
    ap.add_argument("--rrm3-acc", default="P38766")
    ap.add_argument("--root-on", default="Q9H611")
    ap.add_argument("--min-sp", type=int, default=3,
                    help="only tally duplications shared by >= this many species")
    args = ap.parse_args()

    ncbi = NCBITaxa()
    tip2tax, acc2tip = {}, {}
    for r in csv.DictReader(open(args.tip_map), delimiter="\t"):
        tip2tax[r["tip_label"]] = r["taxid"]
        acc2tip[r["accession"]] = r["tip_label"]

    import re
    nwk = re.sub(r"\[&[RU]\]", "", open(args.tree).read()).strip()  # strip dendropy's rooting tag
    gt = PhyloTree(nwk, format=1)  # internal labels are "SH-aLRT/UFBoot" strings, not floats
    og = acc2tip.get(args.root_on)
    nodes = gt.search_nodes(name=og) if og else []
    if nodes and len(gt.children) > 2:
        gt.set_outgroup(nodes[0])
    gt.set_species_naming_function(lambda n: tip2tax.get(n, "0"))

    events = gt.get_descendant_evol_events()
    nD = sum(1 for e in events if e.etype == "D")
    nS = sum(1 for e in events if e.etype == "S")

    clade_dups = Counter()       # clade -> # duplications (shared by >= min_sp species)
    clade_example_sizes = {}
    for e in events:
        if e.etype != "D":
            continue
        seqs = set(e.in_seqs) | set(e.out_seqs)
        taxids = [int(tip2tax[s]) for s in seqs if tip2tax.get(s, "").isdigit()]
        name, rank, nsp = ncbi_lca(ncbi, taxids)
        if nsp >= args.min_sp:
            clade_dups[(name, rank)] += 1
            clade_example_sizes.setdefault((name, rank), []).append(nsp)

    print(f"\n=== {args.tree} ===")
    print(f"events: {nD} duplications, {nS} speciations")
    print(f"duplications shared by >= {args.min_sp} species, by the NCBI clade they map to:")
    for (name, rank), c in clade_dups.most_common(12):
        sizes = clade_example_sizes[(name, rank)]
        print(f"  {c:4d}  {name} ({rank})   [species-spans: "
              f"min {min(sizes)}, max {max(sizes)}]")

    # the anchor node itself
    p, r = acc2tip[args.pif1_acc], acc2tip[args.rrm3_acc]
    anc = gt.get_common_ancestor([p, r])
    taxids = [int(tip2tax[l.name]) for l in anc.get_leaves() if tip2tax.get(l.name, "").isdigit()]
    name, rank, nsp = ncbi_lca(ncbi, taxids)
    print(f"anchor MRCA(Pif1,Rrm3): {len(anc.get_leaves())} genes / {nsp} species -> maps to {name} ({rank})")


if __name__ == "__main__":
    main()
