#!/usr/bin/env python3
"""
20_g4_interface.py  --  Phase 3 / Phase 0: measure the wedge-(or beta-wing)->G4 interface in AF3 models.

For each AlphaFold-Server model of a helicase + T7-AT11 G4 + ATP + 2 K+ complex, quantify whether the
G4-engaging residue actually contacts the quadruplex, and how confident the model is:

  * wedge->G4 distance   : min distance from the engaging residue's functional atoms (Arg guanidinium
                           NH1/NH2/NE/CZ by default) to any DNA guanine; also to the 5'-most guanine
                           (the tetrad the ScPif1 wedge stacks on in 8XAK).
  * K+-in-channel        : for each K+ ion, how many guanine O6 lie within 3.5 A -- a folded G4 channel
                           sandwiches each inter-tetrad K+ between 8 O6 (>=4 per ion => in-channel).
  * interface confidence : ipTM, pTM, ranking_score, and the protein<->DNA chain-pair ipTM / PAE_min
                           read from the matching summary_confidences JSON.

Chains are auto-detected (protein = amino-acid chain; DNA = DA/DG/DC/DT chain; ions = K; ligand = ATP),
so the same script runs on every PIF1 and RecQ job. The engaging residue is given per protein
(--wedge-resid for a single residue, or --wedge-range LO-HI for the RecQ beta-wing loop, which the
structural pilot showed is too mobile to pin to one residue).

Phase 0 (validation gate): run on the ScPif1 job; GO if the top model folds a G4 (K+ in-channel) AND the
wedge contacts the 5' guanines (<= --contact A), matching the experimental 8XAK contact.

Run:  conda run -n pif1 python workflow/20_g4_interface.py \
        --job-dir data/g4_models/pif1_P07271 --wedge-resid 324
"""
import argparse
import glob
import json
import os
import re
import sys

import numpy as np
from Bio.PDB.MMCIF2Dict import MMCIF2Dict

AA = {"ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE","LEU","LYS",
      "MET","PHE","PRO","SER","THR","TRP","TYR","VAL"}
DNT = {"DA","DG","DC","DT"}


def load_atoms(cif):
    """Return a dict of parallel numpy arrays for the atom_site loop."""
    d = MMCIF2Dict(cif)
    asym = np.array(d["_atom_site.label_asym_id"])
    comp = np.array(d["_atom_site.label_comp_id"])
    atom = np.array(d["_atom_site.label_atom_id"])
    seq = np.array(d["_atom_site.label_seq_id"])
    xyz = np.array([d["_atom_site.Cartn_x"], d["_atom_site.Cartn_y"],
                    d["_atom_site.Cartn_z"]], dtype=float).T
    return asym, comp, atom, seq, xyz


def detect_chains(asym, comp):
    prot = dna = None
    for ch in sorted(set(asym)):
        comps = set(comp[asym == ch])
        if comps & AA and prot is None:
            prot = ch
        if comps & DNT and dna is None:
            dna = ch
    ions = sorted(set(asym[comp == "K"]))
    return prot, dna, ions


def mindist(a, b):
    if len(a) == 0 or len(b) == 0:
        return np.inf
    return float(np.sqrt(((a[:, None, :] - b[None, :, :]) ** 2).sum(-1)).min())


def analyze(cif, summary, wedge_resid, wedge_range, wedge_atoms, contact, tetrad_n):
    asym, comp, atom, seq, xyz = load_atoms(cif)
    prot, dna, ions = detect_chains(asym, comp)

    # engaging-residue atoms (single residue, or a residue range for the RecQ beta-wing loop)
    if wedge_range:
        lo, hi = wedge_range
        seqn = np.array([int(s) if s.lstrip("-").isdigit() else -1 for s in seq])
        wmask = (asym == prot) & (seqn >= lo) & (seqn <= hi)
    else:
        wmask = (asym == prot) & (seq == str(wedge_resid))
        if wedge_atoms != "all":
            wmask &= np.isin(atom, wedge_atoms.split(","))
    w = xyz[wmask]

    # DNA guanines; the 5'-most tetrad ~ the tetrad_n lowest-numbered guanines (the face the wedge caps)
    gmask = (asym == dna) & (comp == "DG")
    g = xyz[gmask]
    seqn = np.array([int(s) if s.lstrip("-").isdigit() else -1 for s in seq])
    gres = sorted(set(seqn[gmask].tolist()))
    tetrad5 = gres[:tetrad_n]
    g5 = xyz[gmask & np.isin(seqn, tetrad5)]
    o6 = xyz[gmask & (atom == "O6")]

    # which guanine the wedge actually contacts
    nearest_g, nearest_d = None, np.inf
    for gr in gres:
        dm = mindist(w, xyz[gmask & (seqn == gr)])
        if dm < nearest_d:
            nearest_d, nearest_g = dm, gr

    res = {"model": os.path.basename(cif),
           "wedge_to_G4": round(mindist(w, g), 2), "nearest_G": nearest_g,
           "wedge_to_5ptetrad": round(mindist(w, g5), 2),
           "tetrad5_G": tetrad5, "n_guanines": len(gres)}

    # K+ coordination by guanine O6 (G4 channel check)
    kcoord = []
    for ki in ions:
        k = xyz[asym == ki]
        if len(k) and len(o6):
            d = np.sqrt(((k[0] - o6) ** 2).sum(-1))
            kcoord.append(int((d <= 3.5).sum()))
        else:
            kcoord.append(0)
    res["K_O6_within3.5"] = kcoord
    res["G4_folded"] = any(c >= 4 for c in kcoord)
    res["wedge_contact"] = res["wedge_to_5ptetrad"] <= contact

    # confidence
    if summary and os.path.exists(summary):
        s = json.load(open(summary))
        chains = sorted(set(asym))
        pi, di = chains.index(prot), chains.index(dna)
        cpi = s.get("chain_pair_iptm"); cpp = s.get("chain_pair_pae_min")
        res.update({
            "iptm": s.get("iptm"), "ptm": s.get("ptm"),
            "ranking_score": s.get("ranking_score"), "has_clash": s.get("has_clash"),
            "protDNA_iptm": cpi[pi][di] if cpi else None,
            "protDNA_pae_min": min(x for x in (cpp[pi][di], cpp[di][pi]) if x is not None) if cpp else None,
        })
    return res


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--job-dir", help="a job's extracted dir (processes all model_*.cif)")
    ap.add_argument("--cif", help="a single model CIF instead of --job-dir")
    ap.add_argument("--wedge-resid", type=int, help="engaging residue # in the protein chain (e.g. ScPif1 324)")
    ap.add_argument("--wedge-range", help="LO-HI residue range (RecQ beta-wing loop) instead of a single residue")
    ap.add_argument("--wedge-atoms", default="NH1,NH2,NE,CZ",
                    help="functional atoms of the engaging residue ('all' for whole residue)")
    ap.add_argument("--contact", type=float, default=4.0, help="contact cutoff to the 5' G-tetrad (A)")
    ap.add_argument("--tetrad-n", type=int, default=4, help="# lowest-numbered guanines treated as the 5' tetrad")
    ap.add_argument("--out", help="optional TSV path for the per-model table")
    args = ap.parse_args()

    wrange = tuple(int(x) for x in args.wedge_range.split("-")) if args.wedge_range else None
    if args.cif:
        cifs = [args.cif]
    else:
        cifs = sorted(glob.glob(os.path.join(args.job_dir, "*_model_*.cif")))
    if not cifs:
        sys.exit("[20] no model CIFs found")

    rows = []
    for cif in cifs:
        summary = re.sub(r"_model_(\d+)\.cif$", r"_summary_confidences_\1.json", cif)
        rows.append(analyze(cif, summary, args.wedge_resid, wrange, args.wedge_atoms,
                            args.contact, args.tetrad_n))

    rows.sort(key=lambda r: -(r.get("ranking_score") or 0))
    cols = ["model", "ranking_score", "iptm", "protDNA_iptm", "protDNA_pae_min",
            "wedge_to_5ptetrad", "wedge_to_G4", "nearest_G", "K_O6_within3.5", "G4_folded", "wedge_contact"]
    print("  ".join(c for c in cols))
    for r in rows:
        print("  ".join(str(r.get(c, "")) for c in cols))
    top = rows[0]
    print(f"\n[20] TOP model: G4_folded={top['G4_folded']}, wedge_contact={top['wedge_contact']} "
          f"(wedge->5'tetrad {top['wedge_to_5ptetrad']} A via G{top['nearest_G']}), "
          f"iptm={top.get('iptm')}, prot-DNA PAE_min={top.get('protDNA_pae_min')}")
    if args.out:
        import csv
        allk = list({k for r in rows for k in r})
        with open(args.out, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=cols + [k for k in allk if k not in cols])
            w.writeheader(); w.writerows(rows)
        sys.stderr.write(f"[20] table -> {args.out}\n")


if __name__ == "__main__":
    main()
