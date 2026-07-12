#!/usr/bin/env python3
"""
24_recql4_g4.py  --  R6: do the RQC-less RecQ helicases (RECQL4/5-like) share a way of engaging a G4?

RECQL4/5-like helicases LACK the RQC/β-wing that RECQL1 uses, so there is no known motif column to key off
(scripts 22 and 23 both need one). We therefore scan UNBIASED. For the top AF3 model of each of the 82 jobs:
  * G4_folded     -- did a quadruplex channel form? (a K+ coordinated by >=4 guanine O6 within 3.5 A; same
                     test as script 22, the Phase-0-validated readout)
  * engaging res  -- the protein residue whose SIDE CHAIN comes closest to the 3' G-tetrad (AT11-T7 loads a
                     3' tail, so the 3' face is the loading face). Side-chain (not backbone) distance is used
                     on purpose: the backbone-biased "closest atom" metric is what returned RECQL1 T562 instead
                     of the real Y564 in the earlier RecQ pass, so here we ask specifically which residue
                     *reaches in* to the tetrad.
  * footprint     -- every protein residue with a side-chain atom within --contact of the 3' tetrad, and the
                     chemical make-up of that set (basic / aromatic / acidic / polar / hydrophobic).
  * pos_frac      -- where the engaging residue sits in the modeled chain (0=N-term, 1=C-term); flags cases
                     where a disordered N-terminal arm, not the helicase core, is what grabbed the DNA.
  * confidence    -- iptm, ranking_score, protein<->DNA PAE_min (summary_confidences_0.json)

The aggregate then answers R6 plainly: do RQC-less helicases converge on a chemistry at the 3' face
(a PIF1-like basic wedge? a lingering aromatic?), or is engagement diffuse / not reliably modeled?

Out: results/g4/recql4_interface.tsv (one row per helicase) + a printed summary.

Run:  conda run -n pif1 python workflow/24_recql4_g4.py --struct-dir ~/Desktop/PIF1_RecQ_Structures
"""
import argparse, csv, glob, json, os, re, sys
from collections import defaultdict
import numpy as np

AA3 = {"ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E","GLY":"G","HIS":"H","ILE":"I",
       "LEU":"L","LYS":"K","MET":"M","PHE":"F","PRO":"P","SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V"}
DN = {"DA","DG","DC","DT"}
BACKBONE = {"N", "CA", "C", "O", "OXT"}          # exclude to isolate side-chain engagement
AROM = set("FYWH"); BASIC = set("RK"); ACIDIC = set("DE")
POLAR = set("STNQCYH"); HYDRO = set("AVLIMFWPG")  # rough Venn overlap on Y/H is fine for a coarse footprint tally


def classify(aa):
    if aa in BASIC:  return "basic"
    if aa in AROM:   return "aromatic"
    if aa in ACIDIC: return "acidic"
    if aa in "STNQC": return "polar"
    return "hydrophobic"


def parse_model(cif):
    """Fast targeted CIF scan. Returns K+ xyz, guanine-O6 xyz, {DG seqid: [xyz]},
    {prot seqid: (aa1, [sidechain xyz], [all xyz])}, and per-chain component sets."""
    col = {}
    K, o6 = [], []
    guan = defaultdict(list)
    prot = {}
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
                e = prot.get(r)
                if e is None:
                    e = [AA3[comp], [], []]; prot[r] = e
                e[2].append(xyz)
                if atom not in BACKBONE:
                    e[1].append(xyz)
    return K, o6, guan, prot, chain_comps


def mindist(a, b):
    if not len(a) or not len(b): return np.inf
    a = np.asarray(a); b = np.asarray(b)
    return float(np.sqrt(((a[:, None, :] - b[None, :, :]) ** 2).sum(-1)).min())


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--struct-dir", default=os.path.expanduser("~/Desktop/PIF1_RecQ_Structures"))
    ap.add_argument("--tip-map", default="data/recq/recql4_tip_map.tsv")
    ap.add_argument("--tetrad-n", type=int, default=4)
    ap.add_argument("--contact", type=float, default=4.0)
    ap.add_argument("--out", default="results/g4/recql4_interface.tsv")
    args = ap.parse_args()

    # subfamily (RECQL4-like / RECQL5-like) + lineage labels, keyed by accession (lower)
    tip = {}
    if os.path.exists(args.tip_map):
        for r in csv.DictReader(open(args.tip_map), delimiter="\t"):
            tip[r["accession"].lower()] = (r.get("group", ""), r.get("subphylum", ""), r.get("organism", ""))

    rows = []; n = 0; fail = 0
    for cif in glob.glob(os.path.join(args.struct_dir, "**", "fold_recql4_*_model_0.cif"), recursive=True):
        m = re.search(r"fold_recql4_([a-z0-9]+)_model_0\.cif$", os.path.basename(cif))
        if not m: continue
        acc = m.group(1); n += 1
        try:
            K, o6, guan, prot, chain_comps = parse_model(cif)
            if not guan or not prot:
                fail += 1; continue
            gres = sorted(guan)
            tetrad = gres[-args.tetrad_n:]                         # 3' tetrad = loading face
            tet_atoms = [p for g in tetrad for p in guan[g]]
            g4 = any(sum(1 for q in o6 if (np.array(k) - np.array(q)).dot(np.array(k) - np.array(q)) ** .5 <= 3.5) >= 4
                     for k in K)

            # closest residue BY SIDE CHAIN (Gly falls back to its all-atom set, having no side chain)
            best = (np.inf, None, "?")
            footprint = []
            for r, (aa, sc, allx) in prot.items():
                probe = sc if sc else allx
                dd = mindist(probe, tet_atoms)
                if dd < best[0]: best = (dd, r, aa)
                if dd <= args.contact: footprint.append(aa)
            d, mres, maa = best
            fp_counts = defaultdict(int)
            for aa in footprint: fp_counts[classify(aa)] += 1
            resnums = sorted(prot)
            lo, hi = resnums[0], resnums[-1]
            pos_frac = round((mres - lo) / (hi - lo), 3) if hi > lo and mres is not None else ""

            conf = {}
            sc = cif.replace("model_0.cif", "summary_confidences_0.json")
            if os.path.exists(sc):
                s = json.load(open(sc)); chains = sorted(chain_comps)
                prot_ch = next((c for c in chains if chain_comps[c] & set(AA3)), None)
                dna_ch = next((c for c in chains if chain_comps[c] & DN), None)
                conf["iptm"] = s.get("iptm"); conf["ranking"] = s.get("ranking_score")
                if prot_ch and dna_ch and s.get("chain_pair_pae_min"):
                    pi, di = chains.index(prot_ch), chains.index(dna_ch)
                    cpp = s["chain_pair_pae_min"]; vals = [cpp[pi][di], cpp[di][pi]]
                    conf["pae"] = min(v for v in vals if v is not None) if any(v is not None for v in vals) else None

            grp, sub, org = tip.get(acc, ("", "", ""))
            rows.append({"accession": acc.upper(), "group": grp, "subphylum": sub,
                         "engage_res": mres, "engage_aa": maa,
                         "engage_class": classify(maa) if maa in AA3.values() else "?",
                         "sc_to_3tetrad": round(d, 2) if d != np.inf else "",
                         "contact": int(d <= args.contact), "G4_folded": int(g4),
                         "pos_frac": pos_frac, "n_footprint": len(footprint),
                         "fp_basic": fp_counts["basic"], "fp_aromatic": fp_counts["aromatic"],
                         "fp_acidic": fp_counts["acidic"], "fp_polar": fp_counts["polar"],
                         "fp_hydrophobic": fp_counts["hydrophobic"],
                         "iptm": conf.get("iptm"), "ranking": conf.get("ranking"),
                         "protDNA_pae": conf.get("pae")})
        except Exception as e:
            fail += 1; sys.stderr.write(f"  [warn] {acc}: {e}\n")

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    cols = ["accession","group","subphylum","engage_res","engage_aa","engage_class","sc_to_3tetrad",
            "contact","G4_folded","pos_frac","n_footprint","fp_basic","fp_aromatic","fp_acidic","fp_polar",
            "fp_hydrophobic","iptm","ranking","protDNA_pae"]
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t"); w.writeheader(); w.writerows(rows)

    # ---- printed summary ----
    def frac(pred, sub=None):
        pool = [r for r in rows if (sub is None or sub(r))]
        if not pool: return "n/a"
        k = sum(1 for r in pool if pred(r))
        return f"{k}/{len(pool)} ({100*k/len(pool):.0f}%)"

    folded = [r for r in rows if r["G4_folded"]]
    contacts = [r for r in rows if r["contact"]]
    sys.stderr.write(f"\n[24] scored {len(rows)}/{n} models ({fail} skipped/failed) -> {args.out}\n")
    sys.stderr.write(f"  G4 folded:            {frac(lambda r: r['G4_folded'])}\n")
    sys.stderr.write(f"  side-chain <=4A@3':   {frac(lambda r: r['contact'])}\n")
    sys.stderr.write(f"  ...among folded:      {frac(lambda r: r['contact'], sub=lambda r: r['G4_folded'])}\n")
    def chem_split(pool, label):
        if not pool:
            sys.stderr.write(f"  {label}: (none)\n"); return
        cls = defaultdict(int)
        for r in pool: cls[r["engage_class"]] += 1
        parts = "  ".join(f"{c} {k} ({100*k/len(pool):.0f}%)"
                          for c, k in sorted(cls.items(), key=lambda x: -x[1]))
        sys.stderr.write(f"  {label} (n={len(pool)}):  {parts}\n")

    if contacts:
        chem_split(contacts, "ALL engaging-residue chemistry")
        aa = defaultdict(int)
        for r in contacts: aa[r["engage_aa"]] += 1
        top = " ".join(f"{a}:{k}" for a, k in sorted(aa.items(), key=lambda x: -x[1])[:8])
        sys.stderr.write(f"      top residues: {top}\n")
        nterm = sum(1 for r in contacts if isinstance(r["pos_frac"], float) and r["pos_frac"] < 0.15)
        sys.stderr.write(f"      engaged via N-terminal 15% of chain: {nterm}/{len(contacts)}\n\n")

        # --- by subfamily (RECQL4-like vs RECQL5-like) : the key R6 comparison ---
        for g in ("RECQL4-like", "RECQL5-like"):
            sub = [r for r in rows if r["group"] == g]
            subc = [r for r in sub if r["contact"]]
            sys.stderr.write(f"  [{g}] n={len(sub)}  G4-folded {sum(r['G4_folded'] for r in sub)}/{len(sub)}"
                             f"  contact {len(subc)}/{len(sub)}\n")
            chem_split(subc, f"    {g} chemistry")
        # anchors
        for acc, name in (("O94761", "RECQL4/human"), ("O94762", "RECQL5/human")):
            r = next((x for x in rows if x["accession"] == acc), None)
            if r: sys.stderr.write(f"    anchor {name}: {r['engage_aa']}{r['engage_res']} "
                                   f"({r['engage_class']}) {r['sc_to_3tetrad']}A folded={r['G4_folded']} iptm={r['iptm']}\n")
        # robustness: folded-only, and confident (iptm>=0.70 & pae<=6) subsets
        sys.stderr.write("\n")
        chem_split([r for r in contacts if r["G4_folded"]], "FOLDED-only chemistry")
        conf = [r for r in contacts if (r["iptm"] or 0) >= 0.70 and (r["protDNA_pae"] or 99) <= 6.0]
        chem_split(conf, "CONFIDENT (iptm>=.70,pae<=6) chemistry")
        # lineage: eukaryote vs bacteria (RECQL4/5-like are near-exclusively eukaryotic; verify)
        euk = [r for r in contacts if r["subphylum"]]
        sys.stderr.write(f"  rows with a subphylum label (eukaryotic lineage): {len(euk)}/{len(contacts)}\n")


if __name__ == "__main__":
    main()
