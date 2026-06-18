#!/usr/bin/env python3
"""Walk outward from the duplication node, naming each successive sister clade (what is 'just
outside' the duplication, in nested order), and locate the earliest-diverging Saccharomycotina
(Lipomycetales) gene relative to the duplication clade."""
import os, csv, re
from ete3 import PhyloTree, NCBITaxa
os.chdir(os.path.expanduser("~/pif1-foldtree"))
ncbi=NCBITaxa()
tip2tax,acc2tip={},{}
for r in csv.DictReader(open("data/seqs/tip_map.tsv"),delimiter="\t"):
    tip2tax[r["tip_label"]]=r["taxid"]; acc2tip[r["accession"]]=r["tip_label"]
nwk=re.sub(r"\[&[RU]\]","",open("results/seq_tree/pif1_aa3di.treefile").read()).strip()
t=PhyloTree(nwk,format=1); t.set_species_naming_function(lambda n: tip2tax.get(n,"0"))
leaves=set(t.get_leaf_names())
def tip(a):
    x=acc2tip.get(a); return x if x in leaves else None
pif1_t,rrm3_t,hum_t=tip("P07271"),tip("P38766"),tip("Q9H611")
t.set_outgroup(hum_t)
dup=t.get_common_ancestor([pif1_t,rrm3_t])

def clade_name(node):
    txs=[int(tip2tax[l]) for l in node.get_leaf_names() if tip2tax.get(l,"").isdigit()]
    lin={i:ncbi.get_lineage(i) for i in set(txs) if ncbi.get_lineage(i)}
    if not lin: return "?","?",0
    common=set.intersection(*[set(v) for v in lin.values()])
    ref=next(iter(lin.values()))
    mrca=max(common,key=lambda x:ref.index(x)) if common else None
    return (ncbi.get_taxid_translator([mrca])[mrca] if mrca else "?",
            ncbi.get_rank([mrca]).get(mrca,"?") if mrca else "?", len(lin))

print("Walking outward from the duplication node (nested sisters = what is 'just outside'):")
cur=dup; step=0
while cur.up is not None and step<8:
    par=cur.up
    sibs=[c for c in par.children if c is not cur]
    for s in sibs:
        nm,rk,nsp=clade_name(s)
        tips=s.get_leaf_names()
        sample=[f"{l}({tip2tax.get(l)})" for l in tips[:4]]
        print(f"  step{step}: sister = {nm} ({rk}; {len(tips)} genes/{nsp} sp)  e.g. {sample}")
    cur=par; step+=1

# locate the Lipomycetales gene(s)
print("\nEarliest-diverging Saccharomycotina (Lipomycetales) tips & their placement:")
for l in leaves:
    tx=tip2tax.get(l,"")
    if not tx.isdigit(): continue
    lin=ncbi.get_lineage(int(tx))
    if not lin: continue
    names=set(ncbi.get_taxid_translator(lin).values())
    if "Lipomycetales" in names or "Lipomycetaceae" in names:
        inside = l in set(dup.get_leaf_names())
        # distance/relationship to dup clade
        anc=t.get_common_ancestor([l, pif1_t])
        print(f"  {l} (taxid {tx}) inside_dup_clade={inside}; MRCA(this, ScPif1) subtends {len(anc.get_leaf_names())} genes")
print("SISTERLADDER DONE")
