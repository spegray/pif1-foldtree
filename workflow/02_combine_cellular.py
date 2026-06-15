#!/usr/bin/env python3
"""
02_combine_cellular.py  --  Stage 0/1 curation (keep-all variant).

After switching the family definition from Pfam PF05970 to InterPro IPR048293
("PIF1_RRM3_pfh1" = the CELLULAR PIF1/RRM3/Pfh1 group, which excludes the Helitron
rolling-circle transposon helicases that also carry a bare PF05970 domain), the whole
cellular PIF1 set is small enough (~960 proteins) that we keep ALL of it rather than
subsampling. (The older stratified-subsample path lives in 02_select_representatives.py
and is no longer needed at this size.)

This script concatenates the per-taxon cellular tables produced by workflow/01
(data/seqs/cellular/*.tsv), tags each protein with its broad GROUP and ROLE, and writes
the auditable master list -> data/seqs/selected.tsv. It then prints a breakdown,
including how many species carry >=2 cellular paralogs (the duplication signal).

Groups / roles
--------------
  ascomycota       ingroup   -- contains the PIF1/RRM3 duplication we are dating
  basidiomycota    outgroup  -- sister phylum; helps root + reveals independent duplications
  mucoromycota     outgroup  -- early-diverging fungi
  chytridiomycota  outgroup  -- early-diverging fungi
  human            outgroup  -- deep metazoan anchor (single ancestral PIF1, Q9H611)

Usage
-----
    python workflow/02_combine_cellular.py --indir data/seqs/cellular --out data/seqs/selected.tsv

Uses only the standard library.
"""
import argparse
import csv
import glob
import os
import sys
from collections import defaultdict

# group (= basename of each per-taxon TSV) -> role in the analysis
GROUP_ROLE = {
    "ascomycota": "ingroup",
    "basidiomycota": "outgroup",
    "mucoromycota": "outgroup",
    "chytridiomycota": "outgroup",
    "human": "outgroup",
}


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--indir", default="data/seqs/cellular")
    ap.add_argument("--out", default="data/seqs/selected.tsv")
    args = ap.parse_args()

    files = sorted(glob.glob(os.path.join(args.indir, "*.tsv")))
    if not files:
        sys.exit(f"[02b] no *.tsv in {args.indir} -- run workflow/01 first")

    header = None
    rows = []
    for f in files:
        group = os.path.splitext(os.path.basename(f))[0]
        role = GROUP_ROLE.get(group, "outgroup")
        with open(f) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for r in reader:
                r["group"] = group
                r["role"] = role
                rows.append(r)
            if header is None:
                header = list(reader.fieldnames)

    out_cols = header + ["group", "role"]
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=out_cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)

    # ---- summary ----
    sys.stderr.write(f"\n[02b] combined {len(rows)} proteins -> {args.out}\n\n")
    sys.stderr.write("[02b] proteins / species per group:\n")
    grp_prot = defaultdict(int)
    grp_sp = defaultdict(set)
    for r in rows:
        grp_prot[r["group"]] += 1
        grp_sp[r["group"]].add(r["Organism (ID)"])
    for g in sorted(grp_prot, key=lambda k: -grp_prot[k]):
        sys.stderr.write(f"    {grp_prot[g]:5d} prot  {len(grp_sp[g]):4d} sp  "
                         f"{g}  ({GROUP_ROLE.get(g,'outgroup')})\n")

    # duplication signal: within each subphylum, how many species carry >=2 paralogs?
    sys.stderr.write("\n[02b] cellular-PIF1 copy number by subphylum "
                     "(species with 1 vs >=2 paralogs):\n")
    per_sp = defaultdict(int)          # taxid -> n paralogs
    sp_sub = {}                        # taxid -> subphylum
    for r in rows:
        tid = r["Organism (ID)"]
        per_sp[tid] += 1
        sp_sub[tid] = r.get("subphylum") or r.get("class") or "(unparsed)"
    sub_one = defaultdict(int)
    sub_multi = defaultdict(int)
    for tid, n in per_sp.items():
        if n >= 2:
            sub_multi[sp_sub[tid]] += 1
        else:
            sub_one[sp_sub[tid]] += 1
    subs = sorted(set(sub_one) | set(sub_multi),
                  key=lambda s: -(sub_one[s] + sub_multi[s]))
    sys.stderr.write(f"    {'subphylum/class':28s} {'1-copy':>7s} {'>=2-copy':>9s}\n")
    for s in subs:
        sys.stderr.write(f"    {s:28s} {sub_one[s]:7d} {sub_multi[s]:9d}\n")
    sys.stderr.write("\n[02b] NB: >=2-copy species are where a PIF1-family duplication is "
                     "retained; the gene tree + reconciliation decide which split is the "
                     "true PIF1/RRM3 one.\n")


if __name__ == "__main__":
    main()
