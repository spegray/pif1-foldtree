#!/usr/bin/env python3
"""
19_make_af3_jobs.py  --  Phase 2: emit AlphaFold Server JSON batches for protein-G4 complexes.

Each job is the standardized complex (locked with user):
  - 1 monomer helicase (PIF1- or RecQ-family)
  - the T7-AT11 G4+ssDNA construct (parallel G4 + 5' ssDNA tail; matches ScPif1 co-crystal 8XAK)
  - 1 ATP ligand
  - 2 K+ ions (the two inter-tetrad cations of a 3-tetrad parallel G4)
  - one fixed model seed (run-to-run reproducibility)

Everything except the protein is identical across jobs, so differences reflect the helicase.

AlphaFold SERVER JSON format (top-level array; entities proteinChain/dnaSequence/ligand/ion; NO per-entity
id; ligand by CCD code WITH the CCD_ prefix e.g. "CCD_ATP"; ion by bare code e.g. "K"). VALIDATE ONE JOB in the server UI before
uploading batches in bulk -- a malformed batch costs a day's quota.

Construct: full-length protein unless the WHOLE job (protein + DNA + ATP + 2 K+) would exceed
--max-tokens (5000); only then is it windowed to the helicase domain (corecut +/- pad from tip_map),
keeping the wedge. With the fixed substrate (~68 tokens) the protein budget is ~4932 aa, so nothing
in either family truncates.

Ordering: anchors first, then round-robin across clades, so early batches span the tree of life
(front-loads the conservation signal for the ~30 jobs/day server limit).

Run:  conda run -n pif1 python workflow/19_make_af3_jobs.py --family pif1 \
        --faa data/seqs/selected.faa --tip-map data/seqs/tip_map.tsv
"""
import argparse
import csv
import json
import os
import sys
from collections import defaultdict

from Bio import SeqIO

T7_AT11 = "TTTTTTTTGGTGGTGGTTGTTGTGGTGGTGGTGGT"  # Hu et al. 2024 (8XAK): 8-nt 5' ssDNA + parallel G4
ANCHORS = ["P07271", "P38766", "Q9UUA2", "Q9H611"]  # ScPif1, ScRrm3, SpPfh1, human PIF1

# AF3 server accepts ONLY the 20 standard one-letter codes; a single bad character rejects the whole
# upload batch. UniProt (esp. unreviewed TrEMBL) carries 'X' (unsequenced/uncertain) and occasionally the
# ambiguity codes B/Z/U/O/J. Map ambiguity codes to their nearest standard residue, and X -> G (a glycine
# placeholder approximating AlphaFold's internal 'UNK' handling: the stretch folds as a low-confidence loop
# rather than fabricating a specific side chain). Proteins that are MOSTLY X are excluded instead (--max-x-frac).
STD_AA = set("ACDEFGHIKLMNPQRSTVWY")
SUBMAP = {"B": "D", "Z": "E", "U": "C", "O": "K", "J": "L", "X": "G"}


def sanitize(seq):
    """Return (cleaned sequence, n_substituted). Non-standard -> SUBMAP, else 'G'."""
    out, nsub = [], 0
    for c in seq:
        if c in STD_AA:
            out.append(c)
        else:
            out.append(SUBMAP.get(c, "G")); nsub += 1
    return "".join(out), nsub


def acc_of(header):
    parts = header.split("|")
    return parts[1] if len(parts) >= 2 else header.split()[0]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--family", default="pif1", help="label prefix for job names + files")
    ap.add_argument("--faa", default="data/seqs/selected.faa", help="full-length protein sequences")
    ap.add_argument("--tip-map", default="data/seqs/tip_map.tsv")
    ap.add_argument("--dna", default=T7_AT11)
    ap.add_argument("--max-tokens", type=int, default=5000,
                    help="AF3 server token cap; protein stays FULL-LENGTH unless the job exceeds it, then window")
    ap.add_argument("--atp-atoms", type=int, default=31, help="ATP heavy-atom count (its token cost as a ligand)")
    ap.add_argument("--pad", type=int, default=120, help="residues padded each side of the corecut window if trimming")
    ap.add_argument("--max-x-frac", type=float, default=0.05,
                    help="exclude a protein if >this fraction of its construct is non-standard (else substitute X->G)")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--batch-size", type=int, default=30, help="jobs per JSON file (~1 server day)")
    ap.add_argument("--outdir", default="data/g4/af3_jobs")
    ap.add_argument("--anchors", default=",".join(ANCHORS),
                    help="comma-sep accessions to lead the first batch (family anchors)")
    ap.add_argument("--only-in-tipmap", action="store_true",
                    help="skip FAA sequences absent from --tip-map (RecQ: model only the RQC-bearing set)")
    args = ap.parse_args()
    anchors = [a for a in args.anchors.split(",") if a]

    meta = {}
    for r in csv.DictReader(open(args.tip_map), delimiter="\t"):
        meta[r["accession"]] = r
    seqs = {acc_of(rec.id): str(rec.seq) for rec in SeqIO.parse(args.faa, "fasta")}

    # token budget: full-length unless the WHOLE job (protein + DNA + ATP + 2 K+) exceeds --max-tokens.
    # AF3 tokenizes 1/residue, 1/nucleotide, 1/heavy-atom for ligands, 1/ion.
    overhead = len(args.dna) + args.atp_atoms + 2     # DNA nt + ATP heavy atoms + 2 K+ ions
    budget = args.max_tokens - overhead               # residues of protein that fit
    jobs, manifest = {}, []
    n_full = n_win = n_excl = 0
    subst = []  # (accession, n_residues_substituted) for proteins that needed scrubbing
    for acc, full in seqs.items():
        if args.only_in_tipmap and acc not in meta:
            continue
        m = meta.get(acc, {})
        L = len(full)
        if L <= budget:
            sub, ctype, rng = full, "full", f"1-{L}"
            n_full += 1
        else:
            cf = int(m.get("core_from", 1)); ct = int(m.get("core_to", L))
            lo = max(0, cf - 1 - args.pad); hi = min(L, ct + args.pad)
            if hi - lo > budget:                      # still over -> center the window on the helicase core
                mid = (cf + ct) // 2
                lo = max(0, mid - budget // 2); hi = min(L, lo + budget)
            sub, ctype, rng = full[lo:hi], "domain", f"{lo+1}-{hi}"
            n_win += 1
        # scrub non-standard residues (else AF3 rejects the whole batch); exclude if mostly non-standard
        nonstd = sum(1 for c in sub if c not in STD_AA)
        if nonstd and nonstd / len(sub) > args.max_x_frac:
            n_excl += 1
            if ctype == "full":
                n_full -= 1
            else:
                n_win -= 1
            continue
        sub, nsub = sanitize(sub)
        if nsub:
            subst.append((acc, nsub))
        name = f"{args.family}_{acc}"
        jobs[acc] = {
            "name": name,
            "modelSeeds": [args.seed],
            "sequences": [
                {"proteinChain": {"sequence": sub, "count": 1}},
                {"dnaSequence": {"sequence": args.dna, "count": 1}},
                {"ligand": {"ligand": "CCD_ATP", "count": 1}},  # server requires the CCD_ prefix
                {"ion": {"ion": "K", "count": 2}},
            ],
        }
        manifest.append({"job_name": name, "accession": acc, "construct": ctype,
                         "range": rng, "full_len": L, "model_len": len(sub), "n_subst": nsub,
                         "group": m.get("group", ""), "organism": m.get("organism", "")})

    # ---- ordering: anchors first, then round-robin across clades ----
    by_group = defaultdict(list)
    anchor_set = set(anchors)
    for acc in jobs:
        if acc not in anchor_set:
            by_group[meta.get(acc, {}).get("group", "?")].append(acc)
    for g in by_group:
        by_group[g].sort()
    order = [a for a in anchors if a in jobs]
    groups = sorted(by_group)
    while any(by_group[g] for g in groups):
        for g in groups:
            if by_group[g]:
                order.append(by_group[g].pop(0))

    # ---- write batches + manifest ----
    os.makedirs(args.outdir, exist_ok=True)
    nb = 0
    for i in range(0, len(order), args.batch_size):
        nb += 1
        batch = [jobs[acc] for acc in order[i:i + args.batch_size]]
        with open(os.path.join(args.outdir, f"{args.family}_batch_{nb:03d}.json"), "w") as fh:
            json.dump(batch, fh, indent=2)
    order_idx = {acc: k for k, acc in enumerate(order)}
    manifest.sort(key=lambda r: order_idx.get(r["accession"], 1e9))
    with open(os.path.join(args.outdir, f"{args.family}_jobs_manifest.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["batch", "order", "job_name", "accession", "construct",
                                           "range", "full_len", "model_len", "n_subst", "group", "organism"])
        w.writeheader()
        for k, r in enumerate(manifest):
            r["order"] = k
            r["batch"] = k // args.batch_size + 1
            w.writerow(r)

    sys.stderr.write(f"[19] token budget: {args.max_tokens} - {overhead} overhead = {budget} aa protein cap\n")
    sys.stderr.write(f"[19] {len(jobs)} jobs ({n_full} full-length, {n_win} domain-windowed) "
                     f"-> {nb} batches of {args.batch_size} in {args.outdir}\n")
    sys.stderr.write(f"[19] substrate: {len(args.dna)} nt DNA + ATP + 2 K+, seed {args.seed}\n")
    if subst or n_excl:
        tot = sum(n for _, n in subst)
        sys.stderr.write(f"[19] scrubbed {len(subst)} proteins ({tot} non-standard residues -> X->G etc.); "
                         f"excluded {n_excl} (>{args.max_x_frac:.0%} non-standard). "
                         f"flagged in manifest column n_subst.\n")
        for a, n in sorted(subst, key=lambda x: -x[1])[:8]:
            sys.stderr.write(f"      {a}: {n} substituted\n")
    sys.stderr.write(f"[19] VALIDATE {args.family}_batch_001.json (or one job) in the server UI before bulk upload.\n")


if __name__ == "__main__":
    main()
