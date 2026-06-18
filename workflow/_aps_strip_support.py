#!/usr/bin/env python3
"""Strip SH-aLRT from the IQ-TREE 'SHaLRT/UFBoot' internal labels -> UFBoot-only tree (for GeneRax
--support-threshold), and report the support on the duplication node and its two paralog daughters."""
import re, os, csv
from ete3 import PhyloTree
os.chdir(os.path.expanduser("~/pif1-foldtree"))
IN="results/seq_tree/pif1_aa3di_supported.treefile"
OUT="results/seq_tree/pif1_aa3di_supported.ufboot.treefile"

raw=open(IN).read()
# )SHaLRT/UFBoot:  -> )UFBoot:   (keep the 2nd number = UFBoot)
stripped=re.sub(r"\)([\d.]+)/(\d+):", r")\2:", raw)
open(OUT,"w").write(stripped)
print("wrote",OUT)

# report support on the duplication node (read original SHaLRT/UFBoot labels)
tip2tax,acc2tip={},{}
for r in csv.DictReader(open("data/seqs/tip_map.tsv"),delimiter="\t"):
    tip2tax[r["tip_label"]]=r["taxid"]; acc2tip[r["accession"]]=r["tip_label"]
t=PhyloTree(re.sub(r"\[&[RU]\]","",raw).strip(),format=1)
t.set_outgroup(acc2tip["Q9H611"])
dup=t.get_common_ancestor([acc2tip["P07271"],acc2tip["P38766"]])
def lab(node): return node.name if node.name else "(none)"
print(f"duplication node MRCA(Pif1,Rrm3): {len(dup.get_leaf_names())} genes; SHaLRT/UFBoot = {lab(dup)}")
for c in dup.children:
    which="Pif1" if acc2tip['P07271'] in c.get_leaf_names() else ("Rrm3" if acc2tip['P38766'] in c.get_leaf_names() else "?")
    print(f"  {which} daughter clade: {len(c.get_leaf_names())} genes; SHaLRT/UFBoot = {lab(c)}")
print("STRIP DONE")
