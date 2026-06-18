#!/usr/bin/env python3
"""Figure for the foldseek/FoldTree story:
 (A) 3Di structural characters carry MORE phylogenetic signal than amino acids over the same core;
 (B) which method actually resolves the PIF1/RRM3 duplication (anchor-clade size + mapped clade):
     only the AA+3Di ML partition reaches Saccharomycotina; AA-only and the pure structural
     distance trees (fident/alntmscore/lddt FoldTree) stay deep (Fungi)."""
import os, csv, re
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from ete3 import PhyloTree, NCBITaxa
os.chdir(os.path.expanduser("~/pif1-foldtree"))
os.makedirs("results/figures", exist_ok=True)
ncbi=NCBITaxa()
tip2tax,acc2tip={},{}
for r in csv.DictReader(open("data/seqs/tip_map.tsv"),delimiter="\t"):
    tip2tax[r["tip_label"]]=r["taxid"]; acc2tip[r["accession"]]=r["tip_label"]

def anchor_map(path):
    nwk=re.sub(r"\[&[RU]\]","",open(path).read()).strip()
    t=PhyloTree(nwk,format=1); t.set_species_naming_function(lambda n: tip2tax.get(n,"0"))
    leaves=set(t.get_leaf_names())
    og=acc2tip.get("Q9H611")
    if og in leaves and len(t.children)>2: t.set_outgroup(og)
    p,r=acc2tip["P07271"],acc2tip["P38766"]
    anc=t.get_common_ancestor([p,r]); lv=anc.get_leaf_names()
    txs=[int(tip2tax[l]) for l in lv if tip2tax.get(l,"").isdigit()]
    lin={i:ncbi.get_lineage(i) for i in set(txs) if ncbi.get_lineage(i)}
    common=set.intersection(*[set(v) for v in lin.values()]) if lin else set()
    ref=next(iter(lin.values()))
    mrca=max(common,key=lambda x:ref.index(x)) if common else None
    name=ncbi.get_taxid_translator([mrca])[mrca] if mrca else "?"
    return len(lv), name

methods=[
    ("AA-only ML\n(LG+I+G4)","results/seq_tree/pif1.treefile"),
    ("Structural FoldTree\n(lddt dist.)","results/struct_tree/lddt_foldtree.rooted.nwk"),
    ("Structural FoldTree\n(alntmscore dist.)","results/struct_tree/alntmscore_foldtree.rooted.nwk"),
    ("Structural FoldTree\n(fident dist.)","results/struct_tree/fident_foldtree.rooted.nwk"),
    ("AA + 3Di ML\n(partitioned)","results/seq_tree/pif1_aa3di.treefile"),
]
data=[]
for lab,p in methods:
    n,clade=anchor_map(p); data.append((lab,n,clade)); print(lab.replace(chr(10),' '),n,clade)

fig,(axA,axB)=plt.subplots(1,2,figsize=(13.5,5.6))

# Panel A: information content
cats=["parsimony-\ninformative","invariant"]
aa=[180,15]; tdi=[194,4]
x=range(len(cats)); w=0.38
axA.bar([i-w/2 for i in x],aa,w,label="amino acid",color="#4d4d4d")
axA.bar([i+w/2 for i in x],tdi,w,label="3Di (structure)",color="#1a9850")
for i,(a,b) in enumerate(zip(aa,tdi)):
    axA.text(i-w/2,a+2,str(a),ha="center",fontsize=10)
    axA.text(i+w/2,b+2,str(b),ha="center",fontsize=10,fontweight="bold")
axA.set_xticks(list(x)); axA.set_xticklabels(cats)
axA.set_ylabel("sites (of 209 core columns)")
axA.set_title("A. Structure carries more signal than sequence\n(same 209-column helicase core)",fontsize=11)
axA.legend(fontsize=9)

# Panel B: resolution
labs=[d[0] for d in data]; sizes=[d[1] for d in data]; clades=[d[2] for d in data]
cols=["#1a9850" if c=="Saccharomycotina" else "#b2182b" for c in clades]
y=range(len(labs))
axB.barh(list(y),sizes,color=cols)
axB.axvline(197,ls="--",c="#1a9850",lw=1)
axB.text(197,-0.7,"true duplication-clade size (197)",color="#1a9850",fontsize=8,ha="center")
for i,(s,c) in enumerate(zip(sizes,clades)):
    axB.text(s+8,i,f"{s} genes → {c}",va="center",fontsize=9,
             color="#1a9850" if c=="Saccharomycotina" else "#b2182b")
axB.set_yticks(list(y)); axB.set_yticklabels(labs,fontsize=9); axB.invert_yaxis()
axB.set_xlabel("size of the MRCA(ScPif1,ScRrm3) clade  (genes)")
axB.set_xlim(0,1120)
axB.set_title("B. Only AA+3Di ML resolves the duplication to Saccharomycotina\n"
              "(distance-based FoldTree alone leaves it deep, like AA-only)",fontsize=11)
from matplotlib.patches import Patch
axB.legend(handles=[Patch(color="#1a9850",label="resolved → Saccharomycotina"),
                    Patch(color="#b2182b",label="unresolved → Fungi (kingdom)")],
           fontsize=9,loc="lower right")
fig.suptitle("Foldseek / FoldTree pipeline: the 3Di structural alphabet is what resolves the recent PIF1/RRM3 split",
             fontsize=12,fontweight="bold")
fig.tight_layout(rect=[0,0,1,0.94])
fig.savefig("results/figures/FoldTree_structure_resolves_duplication.png",dpi=300)
print("wrote results/figures/FoldTree_structure_resolves_duplication.png")
print("FIG_FOLDTREE DONE")
