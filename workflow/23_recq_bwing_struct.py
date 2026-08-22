#!/usr/bin/env python3
"""
23_recq_bwing_struct.py  --  RecQ β-wing conservation by STRUCTURAL homology (the rigorous test).

Sequence alignment cannot pin the RecQ β-wing tip, and the min-distance-to-tetrad metric is backbone-biased
(it returned RECQL1 T562, not the real β-wing Y564). So we let STRUCTURE define homology: extract every RecQ
RQC domain from its AF3 model, Foldseek-superpose all of them onto the RECQL1 RQC, read the residue that
aligns to RECQL1's β-wing Y564, and tabulate ITS identity across the family. This asks "is the residue at the
β-wing tip position structurally conserved as an aromatic?" -- independent of AF3's coarse G4 docking.

Out: results/g4/recq_bwing_struct.tsv (accession, group, bwing_res, bwing_aa, aromatic, bwing_to_3tetrad, TMscore).

Run inside the pif1 env (needs foldseek):
    conda run -n pif1 python workflow/23_recq_bwing_struct.py --struct-dir ~/Desktop/PIF1_RecQ_Structures
"""
import argparse, csv, glob, os, re, subprocess, sys
from collections import defaultdict
import numpy as np

AA3 = {"ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E","GLY":"G","HIS":"H","ILE":"I",
       "LEU":"L","LYS":"K","MET":"M","PHE":"F","PRO":"P","SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V"}
AROM = set("FYWH")
REF_ACC, REF_RES = "p46063", 564   # RECQL1 β-wing Y564


def mindist(a, b):
    if not len(a) or not len(b): return None
    a = np.asarray(a); b = np.asarray(b)
    return float(np.sqrt(((a[:, None, :] - b[None, :, :]) ** 2).sum(-1)).min())


def process(cif, lo, hi):
    """Extract the RQC (protein chain, residues lo..hi) as a minimal CIF; return text + ordered residue meta
    + each residue's distance to the 3' G-tetrad."""
    header, col, data = [], {}, []
    res = defaultdict(lambda: [None, []]); dg = defaultdict(list); pchain = None
    with open(cif) as fh:
        for line in fh:
            if line.startswith("_atom_site."):
                col[line.strip()] = len(header); header.append(line.strip()); continue
            if not (line.startswith("ATOM") or line.startswith("HETATM")): continue
            f = line.split()
            comp = f[col["_atom_site.label_comp_id"]]; asym = f[col["_atom_site.label_asym_id"]]
            try:
                xyz = (float(f[col["_atom_site.Cartn_x"]]), float(f[col["_atom_site.Cartn_y"]]), float(f[col["_atom_site.Cartn_z"]]))
            except (ValueError, KeyError):
                continue
            if comp == "DG":
                try: dg[int(f[col["_atom_site.label_seq_id"]])].append(xyz)
                except ValueError: pass
            elif comp in AA3:
                try: r = int(f[col["_atom_site.label_seq_id"]])
                except ValueError: continue
                if lo <= r <= hi:
                    if pchain is None: pchain = asym
                    if asym == pchain:
                        res[r][0] = comp; res[r][1].append(xyz); data.append(line.rstrip("\n"))
    if not res: return None
    txt = "data_x\n#\nloop_\n" + "\n".join(header) + "\n" + "\n".join(data) + "\n#\n"
    ordered = sorted(res)
    tet = [p for g in sorted(dg)[-4:] for p in dg[g]] if dg else []
    dist = {r: mindist(res[r][1], tet) for r in ordered}
    meta = [(r, AA3[res[r][0]]) for r in ordered]
    return txt, meta, dist


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--struct-dir", default=os.path.expanduser("~/Desktop/PIF1_RecQ_Structures"))
    ap.add_argument("--recq-map", default="data/recq/recq_tip_map.tsv")
    ap.add_argument("--rqc-dir", default="/tmp/recq_rqc")
    ap.add_argument("--foldseek", default=os.path.join(os.environ.get("CONDA_PREFIX", ""), "bin", "foldseek"))
    ap.add_argument("--out", default="results/g4/recq_bwing_struct.tsv")
    args = ap.parse_args()

    recq = {r["accession"].lower(): r for r in csv.DictReader(open(args.recq_map), delimiter="\t")}
    os.makedirs(args.rqc_dir, exist_ok=True)
    meta_all, dist_all, group = {}, {}, {}
    n = 0
    for cif in glob.glob(os.path.join(args.struct_dir, "**", "fold_recq_*_model_0.cif"), recursive=True):
        acc = re.search(r"fold_recq_([a-z0-9]+)_model_0", os.path.basename(cif)).group(1)
        if acc not in recq or acc in meta_all: continue
        r = recq[acc]; out = process(cif, int(r["core_from"]), int(r["core_to"]))
        if not out: continue
        txt, meta, dist = out
        open(os.path.join(args.rqc_dir, f"{acc}.cif"), "w").write(txt)
        meta_all[acc] = meta; dist_all[acc] = dist; group[acc] = r.get("group", ""); n += 1
        if n % 250 == 0: sys.stderr.write(f"  extracted {n}\n")
    sys.stderr.write(f"[23] extracted {n} RQC domains -> {args.rqc_dir}\n")

    # RECQL1 Y564 -> its 1-based position in the ordered RQC residue list
    ref = meta_all.get(REF_ACC)
    if not ref: sys.exit(f"[23] reference {REF_ACC} not extracted")
    ref_pos = next((i + 1 for i, (rr, aa) in enumerate(ref) if rr == REF_RES), None)
    sys.stderr.write(f"[23] RECQL1 Y{REF_RES} is ordered-position {ref_pos} ({ref[ref_pos-1][1] if ref_pos else '?'})\n")

    # Foldseek: RECQL1 RQC vs all
    aln = "/tmp/recq_bwing_aln.tsv"
    cmd = [args.foldseek, "easy-search", os.path.join(args.rqc_dir, f"{REF_ACC}.cif"), args.rqc_dir, aln,
           "/tmp/fs_bwing", "--format-output", "target,qstart,tstart,qaln,taln,alntmscore",
           "--exhaustive-search", "1", "-e", "1000", "--max-seqs", "5000"]
    sys.stderr.write("[23] " + " ".join(cmd) + "\n")
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    # map Y564 -> each target's β-wing residue by walking the structural alignment
    rows = {}
    for line in open(aln):
        p = line.rstrip("\n").split("\t")
        if len(p) < 6: continue
        tgt, qs, ts, qaln, taln, tm = p[0], int(p[1]), int(p[2]), p[3], p[4], float(p[5])
        acc = re.match(r"([a-z0-9]+)", os.path.basename(tgt)).group(1)
        if acc not in meta_all or acc in rows: continue
        qp, tp, hit = qs - 1, ts - 1, None
        for qa, ta in zip(qaln, taln):
            if qa != "-": qp += 1
            if ta != "-": tp += 1
            if qa != "-" and qp == ref_pos and ta != "-":
                if tp - 1 < len(meta_all[acc]): hit = meta_all[acc][tp - 1]
                break
        if hit:
            rr, aa = hit
            rows[acc] = {"accession": acc.upper(), "group": group[acc], "bwing_res": rr, "bwing_aa": aa,
                         "aromatic": int(aa in AROM),
                         "bwing_to_3tetrad": round(dist_all[acc].get(rr), 2) if dist_all[acc].get(rr) is not None else "",
                         "TMscore": round(tm, 2)}

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    cols = ["accession", "group", "bwing_res", "bwing_aa", "aromatic", "bwing_to_3tetrad", "TMscore"]
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t"); w.writeheader()
        for acc in sorted(rows): w.writerow(rows[acc])
    sys.stderr.write(f"[23] mapped β-wing residue for {len(rows)} RecQ proteins -> {args.out}\n")


if __name__ == "__main__":
    main()
