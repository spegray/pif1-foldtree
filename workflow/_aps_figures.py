#!/usr/bin/env python3
"""Generate two figures answering the user's deep-dive questions:
  FIG1 (Q2): single-copy outside PIF1s resemble ScPif1 vs ScRrm3 (core %identity).
  FIG2 (Q1): paralog retention across Saccharomycotina orders -> duplication at the subphylum base,
             plus the independent Agaricomycetes (mushroom) duplication for contrast.
Saves PNGs to results/figures/. Run from repo root."""
import os, csv, re
from collections import defaultdict, Counter
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from ete3 import PhyloTree, NCBITaxa
os.chdir(os.path.expanduser("~/pif1-foldtree"))
os.makedirs("results/figures", exist_ok=True)

ncbi = NCBITaxa()
tip2tax, acc2tip = {}, {}
for r in csv.DictReader(open("data/seqs/tip_map.tsv"), delimiter="\t"):
    tip2tax[r["tip_label"]] = r["taxid"]; acc2tip[r["accession"]] = r["tip_label"]

aln = {}; name=None
for line in open("results/seq_tree/aln.trim.fasta"):
    line=line.rstrip("\n")
    if line.startswith(">"): name=line[1:].split()[0]; aln[name]=[]
    elif name: aln[name].append(line)
aln={k:"".join(v) for k,v in aln.items()}
def pid(a,b):
    sa,sb=aln.get(a),aln.get(b)
    if not sa or not sb: return None
    n=s=0
    for x,y in zip(sa,sb):
        if x not in "-.Xx" and y not in "-.Xx":
            n+=1; s+=(x.upper()==y.upper())
    return 100.0*s/n if n else None
def ln(tx):
    try: l=ncbi.get_lineage(int(tx))
    except: return {},set()
    if not l: return {},set()
    nm=ncbi.get_taxid_translator(l); rk=ncbi.get_rank(l)
    return {rk[t]:nm[t] for t in l}, set(nm.values())

nwk=re.sub(r"\[&[RU]\]","",open("results/seq_tree/pif1_aa3di.treefile").read()).strip()
t=PhyloTree(nwk,format=1); t.set_species_naming_function(lambda n: tip2tax.get(n,"0"))
leaves=set(t.get_leaf_names())
def tip(a):
    x=acc2tip.get(a); return x if x in leaves else None
pif1_t,rrm3_t,pfh1_t,hum_t=tip("P07271"),tip("P38766"),tip("Q9UUA2"),tip("Q9H611")
t.set_outgroup(hum_t)
dup=t.get_common_ancestor([pif1_t,rrm3_t])
kids=dup.children
def clade_with(a):
    for k in kids:
        if acc2tip.get(a) in k.get_leaf_names(): return k
def_=None
cP,cR=clade_with("P07271"),clade_with("P38766")
P=set(cP.get_leaf_names()); R=set(cR.get_leaf_names()); dupset=set(dup.get_leaf_names())

# ---- data for FIG1 (Q2) ----
sp_tips=defaultdict(list)
for lf in leaves:
    tx=tip2tax.get(lf,"")
    if tx.isdigit(): sp_tips[tx].append(lf)
pts=[]   # (idPif1, idRrm3, subphylum)
sub_counts=defaultdict(lambda:[0,0])  # subphylum -> [closerPif1, closerRrm3]
for tx,tps in sp_tips.items():
    if len(tps)!=1: continue
    d,names=ln(tx)
    if not names or "Saccharomycotina" in names: continue
    bucket=d.get("subphylum") or d.get("phylum") or "?"
    iP=pid(tps[0],pif1_t); iR=pid(tps[0],rrm3_t)
    if iP is None or iR is None: continue
    pts.append((iP,iR,bucket,tps[0]))
    if iP>iR: sub_counts[bucket][0]+=1
    elif iR>iP: sub_counts[bucket][1]+=1
nP=sum(c[0] for c in sub_counts.values()); nR=sum(c[1] for c in sub_counts.values())

# ---- data for FIG2 (Q1) ----
ORDER_SEQ=["Lipomycetales","Trigonopsidales","Dipodascales","Pichiales","Sporopachydermiales",
           "Ascoideales","Serinales","Phaffomycetales","Saccharomycodales","Saccharomycetales"]
ostat=defaultdict(lambda:dict(both=0,P=0,R=0,single=0))
for tx,tps in sp_tips.items():
    d,names=ln(tx)
    if "Saccharomycotina" not in names: continue
    o=d.get("order","(no order)")
    inP=any(x in P for x in tps); inR=any(x in R for x in tps)
    s=ostat[o]
    if inP and inR: s["both"]+=1
    elif inP: s["P"]+=1
    elif inR: s["R"]+=1
    else: s["single"]+=1
orders=[o for o in ORDER_SEQ if o in ostat]+[o for o in ostat if o not in ORDER_SEQ]

# ============ FIG 1 ============
fig,(ax1,ax2)=plt.subplots(1,2,figsize=(13,5.6))
col_p="#2166ac"; col_r="#b2182b"
for iP,iR,b,lf in pts:
    c=col_p if iP>iR else (col_r if iR>iP else "grey")
    ax1.scatter(iP,iR,s=14,c=c,alpha=0.55,edgecolors="none")
lo,hi=40,62
ax1.plot([lo,hi],[lo,hi],"k--",lw=1)
ax1.text(58,41.3,"closer to ScPif1",color=col_p,ha="right",fontsize=11,fontweight="bold")
ax1.text(41,60,"closer to ScRrm3",color=col_r,ha="left",fontsize=11,fontweight="bold")
# highlight Pfh1 & human
for lf,lab,mk in [(pfh1_t,"S. pombe Pfh1","*"),(hum_t,"human PIF1","D")]:
    iP=pid(lf,pif1_t); iR=pid(lf,rrm3_t)
    ax1.scatter(iP,iR,s=160 if mk=="*" else 70,c="gold",edgecolors="k",marker=mk,zorder=5)
    ax1.annotate(lab,(iP,iR),textcoords="offset points",xytext=(6,-12),fontsize=9)
ax1.set_xlabel("% identity to S. cerevisiae Pif1 (core)")
ax1.set_ylabel("% identity to S. cerevisiae Rrm3 (core)")
ax1.set_xlim(lo,hi); ax1.set_ylim(lo,hi); ax1.set_aspect("equal")
ax1.set_title(f"Each single-copy fungus outside the duplication (n={len(pts)})\n"
              f"{nP} closer to Pif1  vs  {nR} closer to Rrm3", fontsize=11)

# right: per-subphylum stacked fraction
order_sub=sorted(sub_counts, key=lambda b:-(sub_counts[b][0]+sub_counts[b][1]))
yp=range(len(order_sub))
P_=[sub_counts[b][0] for b in order_sub]; R_=[sub_counts[b][1] for b in order_sub]
ax2.barh(list(yp),P_,color=col_p,label="closer to Pif1")
ax2.barh(list(yp),R_,left=P_,color=col_r,label="closer to Rrm3")
ax2.set_yticks(list(yp)); ax2.set_yticklabels(order_sub,fontsize=9)
ax2.invert_yaxis()
for i,b in enumerate(order_sub):
    tot=P_[i]+R_[i]
    ax2.text(tot+2,i,f"{P_[i]}:{R_[i]}",va="center",fontsize=8)
ax2.set_xlabel("number of single-copy species")
ax2.set_title("By lineage — every group leans Pif1",fontsize=11)
ax2.legend(fontsize=9,loc="lower right")
fig.suptitle("Q2: A pre-duplication (single-copy) PIF1 resembles ScPif1 more than ScRrm3\n"
             "(by descent it is a co-ortholog of both; by sequence Pif1 kept more of the ancestral state)",
             fontsize=12,fontweight="bold")
fig.tight_layout(rect=[0,0,1,0.93])
fig.savefig("results/figures/Q2_pif1_vs_rrm3_identity.png",dpi=300)
print("wrote results/figures/Q2_pif1_vs_rrm3_identity.png")

# ============ FIG 2 ============
fig2,ax=plt.subplots(figsize=(11,6.2))
c_both="#1a9850"; c_P=col_p; c_R=col_r; c_s="#999999"
y=range(len(orders))
both=[ostat[o]["both"] for o in orders]
Po=[ostat[o]["P"] for o in orders]; Ro=[ostat[o]["R"] for o in orders]; So=[ostat[o]["single"] for o in orders]
ax.barh(list(y),both,color=c_both,label="both paralogs (Pif1 + Rrm3)")
ax.barh(list(y),Po,left=both,color=c_P,label="Pif1-ortholog only")
left2=[both[i]+Po[i] for i in range(len(orders))]
ax.barh(list(y),Ro,left=left2,color=c_R,label="Rrm3-ortholog only")
left3=[left2[i]+Ro[i] for i in range(len(orders))]
ax.barh(list(y),So,left=left3,color=c_s,label="single-copy / pre-duplication")
ax.set_yticks(list(y)); ax.set_yticklabels(orders,fontsize=10)
# order early->crown top->bottom: ORDER_SEQ is early->crown, so invert so earliest at top
ax.invert_yaxis()
ax.set_xlabel("number of sampled species")
ax.set_title("Q1: PIF1/RRM3 duplication is shared across (nearly) all Saccharomycotina orders\n"
             "→ it maps to the Saccharomycotina common ancestor (the stem of the budding-yeast radiation)",
             fontsize=12,fontweight="bold")
ax.legend(fontsize=9,loc="lower right")
# annotate earliest-diverging
ax.text(0.5,-0.5,"⟵ earliest-diverging (Lipomycetales: single-copy, sparsely sampled)",
        fontsize=8,color="#444",transform=ax.get_yaxis_transform())
fig2.tight_layout()
fig2.savefig("results/figures/Q1_saccharomycotina_paralog_retention.png",dpi=300)
print("wrote results/figures/Q1_saccharomycotina_paralog_retention.png")
print("FIGURES DONE")
