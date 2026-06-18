#!/usr/bin/env python3
"""Annotate a GeneRax run's per-species event counts with NCBI clade names.

Reads <prefix>/species_trees/inferred_species_tree.newick (GeneRax labels internal nodes
node_N) and <prefix>/reconciliations/pif1_speciesEventCounts.txt, maps every species-tree
node to the NCBI clade of its descendant leaves, prints column sums (to identify which column
is duplications vs losses against pif1_eventCounts.txt), the top nodes by each column, and the
event counts on the Saccharomycotina node specifically.

Usage: python workflow/_aps_parse_events.py results/reconciliation/aa3di/spr
"""
import sys, os
from ete3 import NCBITaxa
import dendropy
os.chdir(os.path.expanduser("~/pif1-foldtree"))

pre = sys.argv[1]
ncbi = NCBITaxa()

# --- map species-tree node label -> set of descendant leaf taxids ---
t = dendropy.Tree.get(path=f"{pre}/species_trees/inferred_species_tree.newick",
                      schema="newick", preserve_underscores=True)
def lab_of(nd):
    return nd.taxon.label if nd.is_leaf() else nd.label
node_leaves = {}
for nd in t.postorder_node_iter():
    if nd.is_leaf():
        node_leaves[lab_of(nd)] = {nd.taxon.label}
    else:
        s = set()
        for c in nd.child_nodes():
            s |= node_leaves.get(lab_of(c), set())
        node_leaves[lab_of(nd)] = s

_clade_cache = {}
def clade(taxids):
    key = frozenset(taxids)
    if key in _clade_cache: return _clade_cache[key]
    ids = [int(x) for x in taxids if str(x).isdigit()]
    lin = {i: ncbi.get_lineage(i) for i in set(ids) if ncbi.get_lineage(i)}
    if not lin:
        res = ("?", "?", 0)
    else:
        common = set.intersection(*[set(v) for v in lin.values()])
        ref = next(iter(lin.values()))
        mrca = max(common, key=lambda x: ref.index(x)) if common else None
        name = ncbi.get_taxid_translator([mrca])[mrca] if mrca else "?"
        rank = ncbi.get_rank([mrca]).get(mrca, "?") if mrca else "?"
        res = (name, rank, len(ids))
    _clade_cache[key] = res
    return res

# --- parse event counts ---
rows = []
colsum = None
with open(f"{pre}/reconciliations/pif1_speciesEventCounts.txt") as fh:
    for line in fh:
        p = line.split()
        if not p: continue
        lab = p[0]; nums = list(map(int, p[1:]))
        rows.append((lab, nums))
        if colsum is None: colsum = [0]*len(nums)
        for i, n in enumerate(nums): colsum[i] += n

print(f"=== {pre} ===")
print(f"rows={len(rows)} ncols={len(colsum)} colsums={colsum}")
glob = open(f"{pre}/reconciliations/pif1_eventCounts.txt").read().split()
print(f"global eventCounts: {glob}")
print(f"sample rows: {rows[:4]}")

# top nodes by each column (with clade)
for ci in range(len(colsum)):
    ranked = sorted(rows, key=lambda r: r[1][ci], reverse=True)[:8]
    print(f"\n-- top nodes by column {ci} (colsum={colsum[ci]}) --")
    for lab, nums in ranked:
        if nums[ci] == 0: continue
        name, rank, nsp = clade(node_leaves.get(lab, set()))
        print(f"   {lab:>10}  col{ci}={nums[ci]:<4} -> {name} ({rank}; {len(node_leaves.get(lab,set()))} desc-leaves)")

# Saccharomycotina node(s)
print("\n-- nodes whose clade == Saccharomycotina --")
for lab, nums in rows:
    name, rank, _ = clade(node_leaves.get(lab, set()))
    if name == "Saccharomycotina":
        print(f"   {lab}: counts={nums}  desc-leaves={len(node_leaves.get(lab,set()))}")
print("PARSE DONE")
