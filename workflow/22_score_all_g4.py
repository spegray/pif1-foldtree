#!/usr/bin/env python3
"""
22_score_all_g4.py  --  Phase 3 at scale: score the G4-engaging interface in every AF3 model.

Walks the collated AF3 structure set (one top model, model_0, per helicase) and, for each, measures the
three things the Phase-0 gate said AF3 reports reliably (contact-level, not atomic pose):
  * G4_folded        -- is a quadruplex channel formed? (a K+ coordinated by >=4 guanine O6 within 3.5 A)
  * motif_to_tetrad  -- distance from the G4-engaging motif to the engaged G-tetrad
  * confidence       -- ipTM, ranking_score, protein<->DNA PAE_min (from summary_confidences_0.json)

Family-specific motif (per the two experimental anchors, 8XAK / 9I22):
  PIF1  -- the Arg-WEDGE. A single residue (aligns to ScPif1 R324; per-protein # from pif1_wedge_map.tsv).
           Measured to the 5' G-tetrad (T7-AT11 loads a 5' tail; wedge caps the 5' face).
  RecQ  -- the RQC beta-WING. The tip is a divergent loop that sequence alignment cannot pin (established
           earlier), so we let the STRUCTURE speak: within the RQC envelope (recq_tip_map core_from..core_to)
           we find the residue closest to the 3' G-tetrad (AT11-T7 loads a 3' tail; beta-wing caps the 3'
           face) and report its identity -- testing "a conserved AROMATIC engages the 3' tetrad" without
           trusting a fragile per-residue alignment.

Output: results/g4/interface_all.tsv (one row per helicase) -> feeds the conservation + cross-family summary.

Run:  conda run -n pif1 python workflow/22_score_all_g4.py \
        --struct-dir ~/Desktop/PIF1_RecQ_Structures
"""
import argparse
import csv
import glob
import json
import os
import re
import sys
from collections import defaultdict

import numpy as np

AA3 = {"ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E","GLY":"G","HIS":"H",
       "ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F","PRO":"P","SER":"S","THR":"T","TRP":"W",
       "TYR":"Y","VAL":"V"}
DN = {"DA","DG","DC","DT"}


def parse_model(cif, keep_prot):
    """Fast targeted CIF scan. keep_prot = set of protein residue numbers to retain (AA atoms only)."""
    col = {}
    K, o6 = [], []
    guan = defaultdict(list)                 # DG residue seqid -> [xyz]
    prot = defaultdict(lambda: (None, []))   # residue seqid -> (comp, [xyz])
    chain_comps = defaultdict(set)
    with open(cif) as fh:
        for line in fh:
            if line.startswith("_atom_site."):
                col[line.strip()] = len(col); continue
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            f = line.split()
            comp = f[col["_atom_site.label_comp_id"]]; asym = f[col["_atom_site.label_asym_id"]]
            chain_comps[asym].add(comp)
            try:
                xyz = (float(f[col["_atom_site.Cartn_x"]]), float(f[col["_atom_site.Cartn_y"]]),
                       float(f[col["_atom_site.Cartn_z"]]))
            except (ValueError, KeyError):
                continue
            atom = f[col["_atom_site.label_atom_id"]]
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
                if r in keep_prot:
                    c, pts = prot[r]; pts.append(xyz); prot[r] = (comp, pts)
    return K, o6, guan, prot, chain_comps


def mindist(a, b):
    if not len(a) or not len(b): return np.inf
    a = np.asarray(a); b = np.asarray(b)
    return float(np.sqrt(((a[:, None, :] - b[None, :, :]) ** 2).sum(-1)).min())


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--struct-dir", default=os.path.expanduser("~/Desktop/PIF1_RecQ_Structures"))
    ap.add_argument("--pif1-map", default="data/g4/pif1_wedge_map.tsv")
    ap.add_argument("--recq-map", default="data/recq/recq_tip_map.tsv")
    ap.add_argument("--tetrad-n", type=int, default=4)
    ap.add_argument("--contact", type=float, default=4.0)
    ap.add_argument("--out", default="results/g4/interface_all.tsv")
    args = ap.parse_args()

    pif1 = {r["accession"]: r for r in csv.DictReader(open(args.pif1_map), delimiter="\t")}
    recq = {r["accession"]: r for r in csv.DictReader(open(args.recq_map), delimiter="\t")}

    rows = []; n = 0; fail = 0
    for cif in glob.glob(os.path.join(args.struct_dir, "**", "*_model_0.cif"), recursive=True):
        m = re.search(r"fold_(pif1|recq)_([a-z0-9]+)_model_0\.cif$", os.path.basename(cif))
        if not m: continue
        fam, acc_l = m.group(1), m.group(2)
        n += 1
        try:
            if fam == "pif1":
                rec = next((v for k, v in pif1.items() if k.lower() == acc_l), None)
                if not rec or not rec["wedge_res"]: continue
                wres = int(rec["wedge_res"]); keep = {wres}; group = rec.get("group", "")
            else:
                rec = next((v for k, v in recq.items() if k.lower() == acc_l), None)
                if not rec: continue
                lo, hi = int(rec["core_from"]), int(rec["core_to"]); keep = set(range(lo, hi + 1))
                group = rec.get("group", "")
            K, o6, guan, prot, chain_comps = parse_model(cif, keep)
            if not guan:
                fail += 1; continue
            gres = sorted(guan)
            tetrad = gres[:args.tetrad_n] if fam == "pif1" else gres[-args.tetrad_n:]  # 5' vs 3'
            tet_atoms = [p for g in tetrad for p in guan[g]]
            g4 = any(sum(1 for q in o6 if (np.array(k) - np.array(q)).dot(np.array(k) - np.array(q)) ** .5 <= 3.5) >= 4 for k in K)

            if fam == "pif1":
                comp, pts = prot.get(wres, (None, []))
                d = mindist(pts, tet_atoms); mres, maa = wres, AA3.get(comp, "?")
                cls = maa in ("R", "K")
            else:
                best = (np.inf, None, "?")
                for r, (comp, pts) in prot.items():
                    dd = mindist(pts, tet_atoms)
                    if dd < best[0]: best = (dd, r, AA3.get(comp, "?"))
                d, mres, maa = best; cls = maa in ("F", "Y", "W", "H")

            conf = {}
            sc = cif.replace("model_0.cif", "summary_confidences_0.json")
            if os.path.exists(sc):
                s = json.load(open(sc)); chains = sorted(chain_comps)
                prot_ch = next((c for c in chains if chain_comps[c] & set(AA3)), None)
                dna_ch = next((c for c in chains if chain_comps[c] & DN), None)
                conf["iptm"] = s.get("iptm"); conf["ranking"] = s.get("ranking_score")
                if prot_ch and dna_ch and s.get("chain_pair_pae_min"):
                    pi, di = chains.index(prot_ch), chains.index(dna_ch)
                    cpp = s["chain_pair_pae_min"]
                    vals = [cpp[pi][di], cpp[di][pi]]
                    conf["pae"] = min(v for v in vals if v is not None) if any(v is not None for v in vals) else None
            rows.append({"family": fam, "accession": acc_l.upper(), "group": group,
                         "motif_res": mres, "motif_aa": maa, "class_ok": int(cls),
                         "motif_to_tetrad": round(d, 2), "G4_folded": int(g4),
                         "iptm": conf.get("iptm"), "ranking": conf.get("ranking"),
                         "protDNA_pae": conf.get("pae"), "contact": int(d <= args.contact)})
        except Exception as e:
            fail += 1; sys.stderr.write(f"  [warn] {acc_l}: {e}\n")
        if n % 250 == 0: sys.stderr.write(f"  ...{n} scored\n")

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    cols = ["family","accession","group","motif_res","motif_aa","class_ok","motif_to_tetrad",
            "G4_folded","iptm","ranking","protDNA_pae","contact"]
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t"); w.writeheader(); w.writerows(rows)
    sys.stderr.write(f"[22] scored {len(rows)}/{n} models ({fail} skipped/failed) -> {args.out}\n")


if __name__ == "__main__":
    main()
