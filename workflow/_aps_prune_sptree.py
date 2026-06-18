#!/usr/bin/env python3
"""Prune species-tree leaves not covered by any gene in the GeneRax mapping.

GeneRax 2.1.3 segfaults in Scenario::savePerSpeciesEventsCounts when the species tree
contains a species covered by zero gene families. Our ncbi_species.nwk has 720 leaves but
the 941-gene PIF1 family covers only 719 species, leaving exactly one geneless leaf
(taxid 1041607). Removing geneless leaves is biologically inert for the reconciliation
(a species with no gene only ever contributes inferred losses) and lets GeneRax write the
per-species event counts. Writes ncbi_species.covered.nwk next to the original.
"""
import os, dendropy
os.chdir(os.path.expanduser("~/pif1-foldtree"))

SP_IN  = "data/species_tree/ncbi_species.nwk"
SP_OUT = "data/species_tree/ncbi_species.covered.nwk"
MAP    = "data/species_tree/gene_species.map"

covered = set()
with open(MAP) as fh:
    for line in fh:
        p = line.split()
        if len(p) >= 2:
            covered.add(p[1])

t = dendropy.Tree.get(path=SP_IN, schema="newick", preserve_underscores=True)
before = [l.taxon.label for l in t.leaf_node_iter()]
drop = [x for x in before if x not in covered]
print(f"species-tree leaves: {len(before)}; covered by genes: {len(covered)}; dropping geneless: {drop}")

t.retain_taxa_with_labels([x for x in before if x in covered])
# suppress_unifurcations cleans up the degree-2 node left by removing a leaf
t.suppress_unifurcations()
t.write(path=SP_OUT, schema="newick", unquoted_underscores=True, suppress_rooting=True)

after = [l.taxon.label for l in dendropy.Tree.get(path=SP_OUT, schema='newick', preserve_underscores=True).leaf_node_iter()]
missing = covered - set(after)
print(f"wrote {SP_OUT}: {len(after)} leaves; all covered species present: {len(missing)==0} (missing={list(missing)[:5]})")
print("PRUNE DONE")
