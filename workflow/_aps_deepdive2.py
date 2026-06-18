#!/usr/bin/env python3
"""Refinement: (Q1) Saccharomycotina copy-number & paralog membership by ORDER/family, to pin
stem vs crown and see if the earliest-diverging budding yeasts are two-copy; (Q2) population-level
test over ALL single-copy non-Saccharomycotina species: does the single PIF1 resemble ScPif1 or
ScRrm3 by core %identity?  Run from repo root."""
import sys, os, csv, re
from collections import Counter, defaultdict
from ete3 import PhyloTree, NCBITaxa
os.chdir(os.path.expanduser("~/pif1-foldtree"))

TREE = sys.argv[1] if len(sys.argv) > 1 else "results/seq_tree/pif1_aa3di.treefile"
ALN  = "results/seq_tree/aln.trim.fasta"
PIF1, RRM3, PFH1, HUMAN = "P07271", "P38766", "Q9UUA2", "Q9H611"

ncbi = NCBITaxa()
tip2tax, acc2tip = {}, {}
for r in csv.DictReader(open("data/seqs/tip_map.tsv"), delimiter="\t"):
    tip2tax[r["tip_label"]] = r["taxid"]; acc2tip[r["accession"]] = r["tip_label"]

aln = {}; name=None
for line in open(ALN):
    line=line.rstrip("\n")
    if line.startswith(">"): name=line[1:].split()[0]; aln[name]=[]
    elif name: aln[name].append(line)
aln={k:"".join(v) for k,v in aln.items()}
def pid(a,b):
    sa,sb=aln.get(a),aln.get(b)
    if not sa or not sb: return None
    n=same=0
    for x,y in zip(sa,sb):
        if x not in "-.Xx" and y not in "-.Xx":
            n+=1; same+= (x.upper()==y.upper())
    return 100.0*same/n if n else None

def ln_names(tx):
    try: lin=ncbi.get_lineage(int(tx))
    except: return {}
    if not lin: return {}
    nm=ncbi.get_taxid_translator(lin); rk=ncbi.get_rank(lin)
    return {rk[t]:nm[t] for t in lin}, set(nm.values())

nwk=re.sub(r"\[&[RU]\]","",open(TREE).read()).strip()
t=PhyloTree(nwk,format=1)
t.set_species_naming_function(lambda n: tip2tax.get(n,"0"))
leaves=set(t.get_leaf_names())
def tip(a):
    tp=acc2tip.get(a); return tp if tp in leaves else None
pif1_t,rrm3_t,pfh1_t,hum_t=tip(PIF1),tip(RRM3),tip(PFH1),tip(HUMAN)
if hum_t: t.set_outgroup(hum_t)
dup=t.get_common_ancestor([pif1_t,rrm3_t])
kids=dup.children
def clade_with(a):
    for k in kids:
        if acc2tip.get(a) in k.get_leaf_names(): return k
    return None
cP,cR=clade_with(PIF1),clade_with(RRM3)
P=set(cP.get_leaf_names()); R=set(cR.get_leaf_names())
dupset=set(dup.get_leaf_names())

# ---------- Q1: Saccharomycotina by order, copy number, paralog membership ----------
print(f"TREE={TREE}")
print("\n=== Q1 refine: Saccharomycotina by ORDER (copy number + paralog membership) ===")
# group tips by species, restricted to Saccharomycotina
sacc_sp=defaultdict(list)   # taxid -> [tips]
order_of={}; family_of={}
for lf in leaves:
    tx=tip2tax.get(lf,"")
    if not tx.isdigit(): continue
    d=ln_names(tx)
    if not d: continue
    ranks,names=d
    if "Saccharomycotina" in names:
        sacc_sp[tx].append(lf)
        order_of[tx]=ranks.get("order","(no order)")
        family_of[tx]=ranks.get("family","(no family)")

order_stats=defaultdict(lambda: dict(sp=0, two=0, both=0, onlyP=0, onlyR=0, inDup=0))
for tx,tps in sacc_sp.items():
    o=order_of[tx]; s=order_stats[o]; s["sp"]+=1
    inP=sum(1 for x in tps if x in P); inR=sum(1 for x in tps if x in R)
    indup=sum(1 for x in tps if x in dupset)
    if len(tps)>=2: s["two"]+=1
    if inP and inR: s["both"]+=1
    elif inP: s["onlyP"]+=1
    elif inR: s["onlyR"]+=1
    if indup: s["inDup"]+=1
print(f"{'order':28s} {'#sp':>4} {'2copy':>5} {'both':>5} {'onlyPif1':>8} {'onlyRrm3':>8}")
for o,s in sorted(order_stats.items(), key=lambda kv:-kv[1]["sp"]):
    print(f"{o:28s} {s['sp']:>4} {s['two']:>5} {s['both']:>5} {s['onlyP']:>8} {s['onlyR']:>8}")

# earliest-diverging orders explicitly
EARLY=["Lipomycetales","Trigonopsidales","Dipodascales","Pichiales","Alaninales","Ascoideales",
       "Sporopachydermiales","Serinales","Phaffomycetales","Saccharomycodales","Saccharomycetales"]
print("\n  earliest-diverging Saccharomycotina orders present, with both-paralog status:")
for o in EARLY:
    if o in order_stats:
        s=order_stats[o]
        print(f"    {o:22s} sp={s['sp']:2d} both-paralog={s['both']:2d} onlyPif1={s['onlyP']} onlyRrm3={s['onlyR']}")

# ---------- Q2: population identity over single-copy non-Saccharomycotina species ----------
print("\n=== Q2 refine: single-copy NON-Saccharomycotina species — resemble ScPif1 or ScRrm3? ===")
# count copies per species across whole tree
sp_tips=defaultdict(list)
for lf in leaves:
    tx=tip2tax.get(lf,"")
    if tx.isdigit(): sp_tips[tx].append(lf)

groups=defaultdict(lambda: dict(n=0, closerP=0, closerR=0, idP=[], idR=[], nest_in_dup=0))
for tx,tps in sp_tips.items():
    if len(tps)!=1: continue          # single-copy only
    d=ln_names(tx)
    if not d: continue
    ranks,names=d
    if "Saccharomycotina" in names: continue   # outside the duplication only
    # subphylum bucket
    bucket = ranks.get("subphylum") or ranks.get("phylum") or "?"
    lf=tps[0]
    iP=pid(lf,pif1_t); iR=pid(lf,rrm3_t)
    if iP is None or iR is None: continue
    g=groups[bucket]; g["n"]+=1; g["idP"].append(iP); g["idR"].append(iR)
    if iP>iR: g["closerP"]+=1
    elif iR>iP: g["closerR"]+=1
    if lf in dupset: g["nest_in_dup"]+=1

def mean(x): return sum(x)/len(x) if x else float('nan')
print(f"{'group(subphylum)':26s} {'#1copy':>6} {'closerPif1':>10} {'closerRrm3':>10} {'meanIdPif1':>10} {'meanIdRrm3':>10} {'nestInDup':>9}")
tot=dict(n=0,cP=0,cR=0,idP=[],idR=[])
for b,g in sorted(groups.items(), key=lambda kv:-kv[1]["n"]):
    print(f"{b:26s} {g['n']:>6} {g['closerP']:>10} {g['closerR']:>10} {mean(g['idP']):>10.2f} {mean(g['idR']):>10.2f} {g['nest_in_dup']:>9}")
    tot["n"]+=g["n"]; tot["cP"]+=g["closerP"]; tot["cR"]+=g["closerR"]; tot["idP"]+=g["idP"]; tot["idR"]+=g["idR"]
print(f"{'TOTAL':26s} {tot['n']:>6} {tot['cP']:>10} {tot['cR']:>10} {mean(tot['idP']):>10.2f} {mean(tot['idR']):>10.2f}")

# focus: Ascomycota single-copy (Pezizomycotina+Taphrinomycotina) = the biologically nearest 'just outside'
print("\n  Pfh1 (S.pombe) specifically: %id ScPif1 =", round(pid(pfh1_t,pif1_t),2),
      " %id ScRrm3 =", round(pid(pfh1_t,rrm3_t),2))
print("DEEPDIVE2 DONE")
