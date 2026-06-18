#!/usr/bin/env python3
"""Deep dive: (Q1) where within Saccharomycotina the PIF1/RRM3 duplication maps (stem vs crown,
which sub-lineages have both paralogs), and (Q2) whether the single PIF1 of a species just
outside the duplication resembles ScPif1 or ScRrm3 -- by tree position (orthology), by sequence
identity, and controlling for evolutionary-rate asymmetry.

Usage: python workflow/_aps_deepdive.py [genetree]   (default: AA+3Di tree)
"""
import sys, os, csv, re
from collections import Counter, defaultdict
from ete3 import PhyloTree, NCBITaxa
os.chdir(os.path.expanduser("~/pif1-foldtree"))

TREE = sys.argv[1] if len(sys.argv) > 1 else "results/seq_tree/pif1_aa3di.treefile"
ALN  = "results/seq_tree/aln.trim.fasta"
TIPMAP = "data/seqs/tip_map.tsv"
MANIFEST = "manifest.csv"
PIF1, RRM3, PFH1, HUMAN = "P07271", "P38766", "Q9UUA2", "Q9H611"

ncbi = NCBITaxa()
tip2tax, acc2tip = {}, {}
for r in csv.DictReader(open(TIPMAP), delimiter="\t"):
    tip2tax[r["tip_label"]] = r["taxid"]
    acc2tip[r["accession"]] = r["tip_label"]
tip2acc = {v: k for k, v in acc2tip.items()}
sub = {}
for r in csv.DictReader(open(MANIFEST)):
    sub[r["accession"]] = r.get("subphylum", "") or "(none)"

# alignment
aln = {}
name = None
for line in open(ALN):
    line = line.rstrip("\n")
    if line.startswith(">"):
        name = line[1:].split()[0]; aln[name] = []
    elif name:
        aln[name].append(line)
aln = {k: "".join(v) for k, v in aln.items()}

def pid(a, b):
    """percent identity over columns where BOTH are non-gap."""
    sa, sb = aln.get(a), aln.get(b)
    if not sa or not sb: return None
    n = same = 0
    for x, y in zip(sa, sb):
        if x not in "-.Xx" and y not in "-.Xx":
            n += 1
            if x.upper() == y.upper(): same += 1
    return (100.0*same/n, n) if n else (None, 0)

def lineage_names(taxid):
    try:
        lin = ncbi.get_lineage(int(taxid))
    except Exception:
        return {}
    if not lin: return {}
    names = ncbi.get_taxid_translator(lin)
    ranks = ncbi.get_rank(lin)
    return {ranks[t]: names[t] for t in lin}

def sacc_subclade(taxid):
    """Return a within-Saccharomycotina label (family or order) for a species, else None."""
    ln = lineage_names(taxid)
    names = set(ln.values())
    if "Saccharomycotina" not in names:
        return None
    # prefer family, then order, then genus
    return ln.get("family") or ln.get("order") or ln.get("genus") or "Saccharomycotina(other)"

# ---- tree ----
nwk = re.sub(r"\[&[RU]\]", "", open(TREE).read()).strip()
t = PhyloTree(nwk, format=1)
t.set_species_naming_function(lambda n: tip2tax.get(n, "0"))
leaf_names = set(t.get_leaf_names())

def tip(acc):
    tp = acc2tip.get(acc)
    return tp if tp in leaf_names else None

pif1_t, rrm3_t, pfh1_t, hum_t = tip(PIF1), tip(RRM3), tip(PFH1), tip(HUMAN)
print(f"TREE = {TREE}")
print(f"tips: total={len(leaf_names)}  ScPif1={pif1_t}  ScRrm3={rrm3_t}  Pfh1={pfh1_t}  human={hum_t}")

if hum_t:
    t.set_outgroup(hum_t)

dup = t.get_common_ancestor([pif1_t, rrm3_t])
dup_leaves = dup.get_leaf_names()
dup_taxids = [int(tip2tax[l]) for l in dup_leaves if tip2tax.get(l, "").isdigit()]
print(f"\n=== Q1: duplication node = MRCA(ScPif1,ScRrm3) ===")
print(f"  subtends {len(dup_leaves)} genes / {len(set(dup_taxids))} species")

# the two daughter clades
kids = dup.children
def clade_with(acc):
    for k in kids:
        if acc2tip.get(acc) in k.get_leaf_names():
            return k
    return None
cP = clade_with(PIF1); cR = clade_with(RRM3)

def describe(clade, label):
    leaves = clade.get_leaf_names()
    taxids = [tip2tax[l] for l in leaves if tip2tax.get(l, "").isdigit()]
    saccs = Counter()
    nonsacc = Counter()
    for tx in taxids:
        sc = sacc_subclade(tx)
        if sc: saccs[sc] += 1
        else:
            ln = lineage_names(tx)
            nonsacc[ln.get("subphylum") or ln.get("phylum") or "?"] += 1
    print(f"\n  -- {label} clade: {len(leaves)} genes / {len(set(taxids))} species --")
    print(f"     Saccharomycotina sub-lineages present ({sum(saccs.values())} genes):")
    for k, v in saccs.most_common():
        print(f"        {v:3d}  {k}")
    if nonsacc:
        print(f"     non-Saccharomycotina present: {dict(nonsacc)}")
    return set(saccs)

if cP and cR:
    setP = describe(cP, "Pif1(P07271)")
    setR = describe(cR, "Rrm3(P38766)")
    both = setP & setR
    onlyP = setP - setR
    onlyR = setR - setP
    print(f"\n  Saccharomycotina sub-lineages in BOTH paralog clades ({len(both)}): {sorted(both)}")
    print(f"  only in Pif1 clade: {sorted(onlyP)}")
    print(f"  only in Rrm3 clade: {sorted(onlyR)}")
else:
    print("  !! could not cleanly split duplication node into two paralog daughters")

# sister of the duplication node (the 'just outside' lineage)
print(f"\n=== Q2: the lineage just OUTSIDE the duplication (sister of the dup node) ===")
par = dup.up
sisters = [c for c in par.children if c is not dup] if par else []
sis_reps = []
for s in sisters:
    sl = s.get_leaf_names()
    stax = [tip2tax[l] for l in sl if tip2tax.get(l, "").isdigit()]
    name, rank, _ = (lambda ids: (lambda lin: (
        (ncbi.get_taxid_translator([max(set.intersection(*[set(v) for v in lin.values()]),
          key=lambda x: next(iter(lin.values())).index(x))])[
          max(set.intersection(*[set(v) for v in lin.values()]),
          key=lambda x: next(iter(lin.values())).index(x))],
         "", len(lin)) if lin else ("?","",0)))(
        {i: ncbi.get_lineage(i) for i in set(int(x) for x in ids if str(x).isdigit()) if ncbi.get_lineage(i)}))(stax)
    print(f"  sister clade: {len(sl)} genes -> ~{name}; sample tips: {sl[:6]}")
    sis_reps += sl[:5]

# Q2 sequence identity + patristic distance for chosen 'outside' orthologs
print(f"\n=== Q2: does an outside single-copy PIF1 resemble ScPif1 or ScRrm3? ===")
print(f"  baseline: %id(ScPif1,ScRrm3) = {pid(pif1_t, rrm3_t)}")
cands = []
if pfh1_t: cands.append(("S.pombe Pfh1 (Q9UUA2, Taphrinomycotina)", pfh1_t))
if hum_t:  cands.append(("Human PIF1 (Q9H611, outgroup root)", hum_t))
# add a couple of sister-lineage reps (closest pre-duplication relatives in the tree)
for r in sis_reps[:3]:
    cands.append((f"sister-lineage rep {r} ({sub.get(tip2acc.get(r,''),'?')})", r))

for label, tp in cands:
    idP = pid(tp, pif1_t); idR = pid(tp, rrm3_t)
    try:
        dP = t.get_distance(tp, pif1_t); dR = t.get_distance(tp, rrm3_t)
    except Exception:
        dP = dR = None
    verdict = "?"
    if idP and idR:
        verdict = "ScPif1" if idP[0] > idR[0] else ("ScRrm3" if idR[0] > idP[0] else "tie")
    print(f"\n  {label}")
    print(f"     %id to ScPif1 = {idP}    %id to ScRrm3 = {idR}   -> more similar to: {verdict}")
    print(f"     patristic dist to ScPif1 = {dP}   to ScRrm3 = {dR}")

# rate-asymmetry control: branch length from dup node to each anchor, and clade mean root-to-tip
print(f"\n=== rate-asymmetry control (interprets Q2) ===")
print(f"  dist(dupnode -> ScPif1) = {dup.get_distance(pif1_t)}")
print(f"  dist(dupnode -> ScRrm3) = {dup.get_distance(rrm3_t)}")
if cP and cR:
    def mean_rt(clade):
        ds = [clade.get_distance(l) for l in clade.get_leaf_names()]
        return sum(ds)/len(ds) if ds else None
    print(f"  mean root-to-tip within Pif1 clade = {mean_rt(cP):.4f}")
    print(f"  mean root-to-tip within Rrm3 clade = {mean_rt(cR):.4f}")
print("DEEPDIVE DONE")
