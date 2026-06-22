#!/usr/bin/env python3
"""
18_wedge_conservation.py  --  Phase 1: sequence conservation of a G4-engaging anchor residue.

Maps an anchor residue (default: ScPif1 R324, the 1A-domain wedge Arg that cation-π stacks on the
5'-most G-tetrad; Hu et al. 2024, PDB 8XAK) onto the existing family alignment and tabulates the
residue at that column across all homologs -- the sequence half of "how conserved is the wedge".

Uses the UNTRIMMED core alignment (results/seq_tree/aln.fasta) so no column is lost to trimAl. The
anchor residue is given in FULL-sequence numbering; we convert via the corecut coordinates in
tip_map.tsv (core_from), then find the alignment column of that core residue and read it across the family.

If the anchor residue falls OUTSIDE the corecut envelope, the script says so (then we extend the
alignment N-terminally to cover the wedge) rather than silently mis-mapping.

Run:  conda run -n pif1 python workflow/18_wedge_conservation.py            # PIF1 / R324
      conda run -n pif1 python workflow/18_wedge_conservation.py --anchor-acc <acc> --anchor-res <n> ...
"""
import argparse
import csv
import sys
from collections import Counter

from Bio import SeqIO


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--anchor-acc", default="P07271")
    ap.add_argument("--anchor-res", type=int, default=324, help="1-based residue # in the FULL sequence")
    ap.add_argument("--anchor-aa", default="R", help="expected residue identity (sanity check)")
    ap.add_argument("--full-faa", default="data/seqs/selected.faa")
    ap.add_argument("--aln", default="results/seq_tree/aln.fasta")
    ap.add_argument("--tip-map", default="data/seqs/tip_map.tsv")
    args = ap.parse_args()

    # tip_map: accession -> tip_label, core_from, core_to, group
    acc2 = {}
    for r in csv.DictReader(open(args.tip_map), delimiter="\t"):
        acc2[r["accession"]] = r
    if args.anchor_acc not in acc2:
        sys.exit(f"[18] anchor {args.anchor_acc} not in tip_map")
    a = acc2[args.anchor_acc]
    tip, cf, ct = a["tip_label"], int(a["core_from"]), int(a["core_to"])

    # full anchor sequence (verify the residue identity + numbering)
    full = {acc_of(rec.id): str(rec.seq) for rec in SeqIO.parse(args.full_faa, "fasta")}
    fseq = full.get(args.anchor_acc)
    if fseq:
        obs = fseq[args.anchor_res - 1] if args.anchor_res <= len(fseq) else "?"
        sys.stderr.write(f"[18] {args.anchor_acc} full len {len(fseq)}; residue {args.anchor_res} = "
                         f"{obs} (expected {args.anchor_aa})\n")
    sys.stderr.write(f"[18] corecut envelope for {args.anchor_acc}: {cf}-{ct}\n")

    if not (cf <= args.anchor_res <= ct):
        print(f"\n[18] *** anchor residue {args.anchor_res} is OUTSIDE the corecut ({cf}-{ct}). ***")
        print("     The wedge is not covered by the current alignment -> extend the corecut "
              "N-terminally (re-run 05 with padding / a wider envelope) before mapping.")
        return

    core_idx = args.anchor_res - cf            # 0-based index into the (ungapped) core sequence

    # aligned anchor core -> find the alignment column of core residue `core_idx`
    aln = {rec.id: str(rec.seq) for rec in SeqIO.parse(args.aln, "fasta")}
    if tip not in aln:
        sys.exit(f"[18] anchor tip {tip} not in alignment")
    ng, col = -1, None
    for ci, c in enumerate(aln[tip]):
        if c != "-":
            ng += 1
            if ng == core_idx:
                col = ci
                break
    if col is None:
        sys.exit("[18] could not locate the anchor residue's column")
    sys.stderr.write(f"[18] wedge maps to alignment column {col} (anchor residue there: {aln[tip][col]})\n")

    # tabulate that column across the family
    column = [seq[col] for seq in aln.values()]
    n = len(column)
    cnt = Counter(column)
    gaps = cnt.get("-", 0)
    nongap = n - gaps
    print(f"\n=== wedge-column conservation (anchor {args.anchor_acc} {args.anchor_aa}{args.anchor_res}) ===")
    print(f"  {n} sequences; {gaps} gap, {nongap} residue")
    rfrac = cnt.get(args.anchor_aa, 0) / nongap if nongap else 0
    print(f"  fraction = {args.anchor_aa}: {cnt.get(args.anchor_aa,0)}/{nongap} = {rfrac:.1%} (of non-gap)")
    print("  residue distribution: " + ", ".join(f"{aa}:{c}" for aa, c in cnt.most_common(8)))

    # by broad clade (tip_map 'group')
    tip2grp = {acc2[acc]["tip_label"]: acc2[acc].get("group", "") for acc in acc2}
    by = {}
    for tipid, seq in aln.items():
        g = tip2grp.get(tipid, "?")
        by.setdefault(g, Counter())[seq[col]] += 1
    print("  by clade (fraction = anchor residue, of non-gap):")
    for g, c in sorted(by.items()):
        ngc = sum(v for k, v in c.items() if k != "-")
        fr = c.get(args.anchor_aa, 0) / ngc if ngc else 0
        print(f"    {g:16s} {args.anchor_aa}={c.get(args.anchor_aa,0)}/{ngc} ({fr:.0%})  "
              f"top: {', '.join(f'{a}:{n2}' for a,n2 in c.most_common(3))}")


def acc_of(header):
    parts = header.split("|")
    return parts[1] if len(parts) >= 2 else header.split()[0]


if __name__ == "__main__":
    main()
