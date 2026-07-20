#!/usr/bin/env python3
"""
score_ancestor_fold.py -- score the ancestral PIF1 core AF3 folds with the SAME metric as the 956 modern
models (workflow/22): G4_folded (a K+ coordinated by >=4 guanine O6 within 3.5 A), and the wedge->5'-tetrad
min distance. The wedge is residue 89 of the 206-aa ancestral core (the reconstructed Arg; see
extract_ancestor.R). Scores all 5 AF3 models of the +G4 job for a pose-robust read, and confirms the wedge
in the apo job. Compare against the modern PIF1 aggregate (median 5.5 A; <=4 A 25%, <=6 A 57%, <=8 A 83%;
ScPif1 anchor R324 -> 1.98 A; G4 folded 94%) from results/g4/CONSERVATION_SUMMARY.md.

Run from repo root:  python3 results/asr/score_ancestor_fold.py
"""
import glob
import json
import os
import numpy as np
from collections import defaultdict

AA3 = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G",
       "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P", "SER": "S",
       "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"}
WEDGE = 89  # ancestral-core position of the reconstructed Arg-wedge


def parse(cif):
    col = {}
    K, o6 = [], []
    guan = defaultdict(list)
    prot = defaultdict(lambda: (None, []))
    wedge_plddt = None
    with open(cif) as fh:
        for line in fh:
            if line.startswith("_atom_site."):
                col[line.strip()] = len(col); continue
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            f = line.split()
            comp = f[col["_atom_site.label_comp_id"]]
            atom = f[col["_atom_site.label_atom_id"]]
            try:
                xyz = (float(f[col["_atom_site.Cartn_x"]]), float(f[col["_atom_site.Cartn_y"]]),
                       float(f[col["_atom_site.Cartn_z"]]))
            except (ValueError, KeyError):
                continue
            if comp == "K":
                K.append(xyz)
            elif comp == "DG":
                try: g = int(f[col["_atom_site.label_seq_id"]])
                except ValueError: continue
                guan[g].append(xyz)
                if atom == "O6": o6.append(xyz)
            elif comp in AA3:
                try: r = int(f[col["_atom_site.label_seq_id"]])
                except ValueError: continue
                c, pts = prot[r]; pts.append(xyz); prot[r] = (comp, pts)
                if r == WEDGE and atom == "CA":
                    try: wedge_plddt = float(f[col["_atom_site.B_iso_or_equiv"]])
                    except (ValueError, KeyError): pass
    return K, o6, guan, prot, wedge_plddt


def mindist(a, b):
    if not len(a) or not len(b): return np.inf
    a = np.asarray(a); b = np.asarray(b)
    return float(np.sqrt(((a[:, None, :] - b[None, :, :]) ** 2).sum(-1)).min())


def score(cif, with_g4):
    K, o6, guan, prot, wplddt = parse(cif)
    comp, pts = prot.get(WEDGE, (None, []))
    waa = AA3.get(comp, "?")
    n = int(cif.split("_model_")[1].split(".")[0])
    sc = cif.replace("model_%d.cif" % n, "summary_confidences_%d.json" % n)
    iptm = ptm = None
    if os.path.exists(sc):
        s = json.load(open(sc)); iptm = s.get("iptm"); ptm = s.get("ptm")
    row = {"model": n, "wedge_aa": waa, "wedge_plddt": wplddt, "iptm": iptm, "ptm": ptm}
    if with_g4:
        g4 = any(sum(1 for q in o6 if np.linalg.norm(np.array(k) - np.array(q)) <= 3.5) >= 4 for k in K)
        gres = sorted(guan)
        tet = [p for g in gres[:4] for p in guan[g]]      # 5' tetrad (as workflow/22)
        row["G4_folded"] = int(g4)
        row["wedge_to_5tetrad"] = round(mindist(pts, tet), 2)
    return row


base = "results/asr/af3/folds"
print("=== +G4 (holo): all 5 models ===")
g4rows = [score(c, True) for c in sorted(glob.glob(base + "/g4/*_model_*.cif"))]
for r in g4rows:
    print("  model %d: wedge=%s pLDDT=%.1f  G4_folded=%d  wedge->5'tetrad=%.2f A  ipTM=%.2f" %
          (r["model"], r["wedge_aa"], r["wedge_plddt"] or -1, r["G4_folded"], r["wedge_to_5tetrad"], r["iptm"] or -1))
d = [r["wedge_to_5tetrad"] for r in g4rows]
folded = [r for r in g4rows if r["G4_folded"]]
print("  --> G4 folded in %d/5 models; wedge->5'tetrad median %.2f A (min %.2f); of folded, best %.2f A" %
      (len(folded), float(np.median(d)), min(d), min([r["wedge_to_5tetrad"] for r in folded]) if folded else -1))

print("=== apo (-G4): wedge present + confidence ===")
aporows = [score(c, False) for c in sorted(glob.glob(base + "/apo/*_model_*.cif"))]
for r in aporows:
    print("  model %d: wedge=%s pLDDT=%.1f  pTM=%.2f" % (r["model"], r["wedge_aa"], r["wedge_plddt"] or -1, r["ptm"] or -1))

import csv
with open("results/asr/ancestor_fold_scores.tsv", "w", newline="") as fh:
    cols = ["job", "model", "wedge_aa", "wedge_plddt", "G4_folded", "wedge_to_5tetrad", "iptm", "ptm"]
    w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t"); w.writeheader()
    for r in g4rows: w.writerow({**{"job": "g4"}, **r})
    for r in aporows: w.writerow({**{"job": "apo"}, **r})
print("\nwrote results/asr/ancestor_fold_scores.tsv")
