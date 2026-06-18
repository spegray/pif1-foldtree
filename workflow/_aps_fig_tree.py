#!/usr/bin/env python3
"""Annotated-tree figures:
   (1) tree_schematic.png  -- clean cladogram telling the whole story (for the explainer slide).
   (2) tree_dupclade_real.png -- the REAL 197-gene MRCA(ScPif1,ScRrm3) clade from the AA+3Di tree,
       rendered as a circular cladogram with the Pif1 and Rrm3 ortholog clades colored.
"""
import os, csv, re, math
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrow, Polygon, Patch
from ete3 import PhyloTree
os.chdir(os.path.expanduser("~/pif1-foldtree"))
os.makedirs("results/figures", exist_ok=True)

BLUE="#2166ac"; RED="#b2182b"; GREY="#9a9a9a"; GREEN="#1a9850"; GOLD="#f5a800"
tip2tax,acc2tip={},{}
for r in csv.DictReader(open("data/seqs/tip_map.tsv"),delimiter="\t"):
    tip2tax[r["tip_label"]]=r["taxid"]; acc2tip[r["accession"]]=r["tip_label"]

# ---------------- Figure 1: schematic ----------------
fig,ax=plt.subplots(figsize=(11,6.6)); ax.axis("off")
ax.set_xlim(0,10); ax.set_ylim(0,10)
def line(x0,y0,x1,y1,c="k",lw=2): ax.plot([x0,x1],[y0,y1],c=c,lw=lw,solid_capstyle="round",zorder=2)
def tri(x,y,w,h,c,label,sub):
    ax.add_patch(Polygon([[x,y],[x+w,y+h/2],[x+w,y-h/2]],closed=True,fc=c,ec="k",lw=1,alpha=.92,zorder=3))
    ax.text(x+w+0.12,y,label,va="center",ha="left",fontsize=11,fontweight="bold")
    ax.text(x+w+0.12,y-0.42,sub,va="center",ha="left",fontsize=8.5,color="#444")

# backbone: ancestral single-copy PIF1
line(0.5,5,2.4,5,"k",2.5)
ax.text(0.5,5.35,"ancestral fungal PIF1\n(single copy)",fontsize=9,color="#333")
# split: outgroup fungi (up) vs Saccharomycotina lineage (down)
line(2.4,5,2.4,7.6); line(2.4,5,2.4,2.2)
# outgroup triangle (1 copy)
line(2.4,7.6,4.6,7.6,"k",2)
tri(4.6,7.6,1.5,1.5,GREY,"All other fungi — 1 copy",
    "Pezizomycotina, Taphrinomycotina (Pfh1), Basidiomycota, early-diverging fungi")
# Saccharomycotina lineage down to the duplication
line(2.4,2.2,4.0,2.2,"k",2.5)
# duplication star
ax.scatter([4.0],[2.2],marker="*",s=900,c=GOLD,edgecolors="k",zorder=5)
ax.text(4.0,1.35,"PIF1/RRM3 duplication\n= Saccharomycotina ancestor",ha="center",fontsize=9.5,
        fontweight="bold",color="#a06b00")
# two paralog clades
line(4.0,2.2,4.0,3.4); line(4.0,2.2,4.0,1.0)
line(4.0,3.4,6.0,3.4,BLUE,2); line(4.0,1.0,6.0,1.0,RED,2)
tri(6.0,3.4,1.7,1.4,BLUE,"Pif1 orthologs",  "incl. S. cerevisiae Pif1 — across ~all Saccharomycotina")
tri(6.0,1.0,1.7,1.4,RED, "Rrm3 orthologs",  "incl. S. cerevisiae Rrm3 — across ~all Saccharomycotina")
# independent mushroom duplication (inset, top-right)
line(6.35,8.7,7.2,8.7,"k",1.6)
ax.scatter([7.2],[8.7],marker="*",s=300,c=GOLD,edgecolors="k",zorder=5)
line(7.2,8.7,8.0,9.05,GREEN,1.6); line(7.2,8.7,8.0,8.35,GREEN,1.6)
ax.text(8.1,8.7,"mushrooms (Agaricomycetes)\nINDEPENDENT duplication",va="center",fontsize=8.5,color=GREEN)
ax.text(6.3,9.4,"convergence, not shared ancestry:",fontsize=8.5,color=GREEN,fontweight="bold")
ax.set_title("PIF1/RRM3 arose once, in the Saccharomycotina (budding-yeast) ancestor\n"
             "before the radiation; an outside species' single PIF1 is a co-ortholog of both (Pif1-like by sequence)",
             fontsize=12,fontweight="bold")
fig.tight_layout()
fig.savefig("results/figures/tree_schematic.png",dpi=300); print("wrote tree_schematic.png")

# ---------------- Figure 2: real duplication clade, circular ----------------
nwk=re.sub(r"\[&[RU]\]","",open("results/seq_tree/pif1_aa3di.treefile").read()).strip()
t=PhyloTree(nwk,format=1)
leaves_all=set(t.get_leaf_names())
t.set_outgroup(acc2tip["Q9H611"])
dup=t.get_common_ancestor([acc2tip["P07271"],acc2tip["P38766"]])
kids=dup.children
def clade_with(a):
    for k in kids:
        if acc2tip.get(a) in k.get_leaf_names(): return k
cP,cR=clade_with("P07271"),clade_with("P38766")
Pset=set(cP.get_leaf_names()); Rset=set(cR.get_leaf_names())

leaves=dup.get_leaves(); n=len(leaves)
ang={lf.name:2*math.pi*i/n for i,lf in enumerate(leaves)}
# levels from dup
for node in dup.traverse("preorder"):
    node.add_feature("lvl",0 if node is dup else node.up.lvl+1)
maxlvl=max(node.lvl for node in dup.traverse())
for node in dup.traverse("postorder"):
    node.add_feature("ang", ang[node.name] if node.is_leaf() else
                     sum(c.ang for c in node.children)/len(node.children))
Rmax=1.0
def rad(node): return Rmax if node.is_leaf() else Rmax*node.lvl/maxlvl
def grp(node):
    lv=set(node.get_leaf_names())
    if lv<=Pset: return BLUE
    if lv<=Rset: return RED
    return GREY
def pol(r,a): return (r*math.cos(a), r*math.sin(a))

fig2,ax=plt.subplots(figsize=(8.4,8.4)); ax.axis("off"); ax.set_aspect("equal")
# edges
for node in dup.traverse():
    if node is dup: continue
    c=grp(node); rp=rad(node.up); rc=rad(node); a=node.ang
    x0,y0=pol(rp,a); x1,y1=pol(rc,a)
    ax.plot([x0,x1],[y0,y1],c=c,lw=0.7,zorder=2)
# arcs at each internal node's radius spanning its children angles
for node in dup.traverse():
    if node.is_leaf(): continue
    r=rad(node); cs=[c.ang for c in node.children]
    a0,a1=min(cs),max(cs); c=grp(node)
    aa=[a0+(a1-a0)*k/30 for k in range(31)]
    ax.plot([r*math.cos(a) for a in aa],[r*math.sin(a) for a in aa],c=c,lw=0.7,zorder=2)
# center duplication star
ax.scatter([0],[0],marker="*",s=700,c=GOLD,edgecolors="k",zorder=6)
# mark anchors
for acc,lab,col in [("P07271","ScPif1",BLUE),("P38766","ScRrm3",RED)]:
    lf=acc2tip[acc]; a=ang[lf]; x,y=pol(Rmax*1.04,a)
    ax.scatter(*pol(Rmax,a),s=42,c=col,edgecolors="k",zorder=7)
    ax.text(x,y,lab,fontsize=11,fontweight="bold",color=col,
            ha="left" if -math.pi/2<a<math.pi/2 else "right",
            va="center",rotation=math.degrees(a) if -math.pi/2<a<math.pi/2 else math.degrees(a)+180,
            rotation_mode="anchor")
# clade labels placed at each clade's angular centre, coloured to match its branches
def clade_label(clade,text,color):
    angs=[ang[l] for l in clade.get_leaf_names()]
    ma=math.atan2(sum(math.sin(a) for a in angs),sum(math.cos(a) for a in angs))
    x,y=pol(1.34,ma)
    ax.text(x,y,text,ha="center",va="center",color=color,fontsize=11.5,fontweight="bold")
clade_label(cP,f"Pif1 ortholog clade\n({len(cP.get_leaf_names())} genes)",BLUE)
clade_label(cR,f"Rrm3 ortholog clade\n({len(cR.get_leaf_names())} genes)",RED)
ax.text(0,-0.13,"duplication",ha="center",fontsize=8,color="#a06b00")
ax.set_xlim(-1.62,1.62); ax.set_ylim(-1.62,1.62)
ax.set_title("The real duplication clade (AA+3Di tree): 197 genes / 103 Saccharomycotina species\n"
             "two mirror-image paralog clades radiating from one duplication",
             fontsize=11.5,fontweight="bold")
fig2.tight_layout()
fig2.savefig("results/figures/tree_dupclade_real.png",dpi=300); print("wrote tree_dupclade_real.png")
print("FIG_TREE DONE")
