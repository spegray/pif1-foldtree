#!/usr/bin/env python3
"""
13_build_aa3di.py  --  build the AA+3Di concatenated alignment + IQ-TREE partition (FoldTree method).

The deep PIF1/RRM3 nodes are unresolved because the 209-col AA core has saturated. The FoldTree
paper's remedy: add the 3Di structural alphabet, which stays informative in the twilight zone, as a
SECOND partition with the SAME column/gap structure as the AA alignment, then run IQ-TREE with a
partition model (LG for AA + the 3DiPhy 3Di substitution matrix for 3Di).

Method (matches Moi et al. + the fold_tree pipeline):
  - 3Di per structure comes from `foldseek createdb` (data/3di/coresdb).
  - We DON'T re-align 3Di independently; we map each residue's 3Di state onto the EXISTING AA
    alignment columns (validated: trimAl -colnumbering reproduces aln.trim.fasta exactly), so AA and
    3Di share an identical gap pattern -- exactly the paper's requirement.
  - Concatenate -> alnAA_3di.fasta (209 AA + 209 3Di = 418 cols); write a partition NEXUS.

Outputs (results/seq_tree/):
  3di.trim.fasta        the column-matched 3Di alignment (209 cols, same gaps as aln.trim.fasta)
  alnAA_3di.fasta       concatenated 957 x 418
  aa3di_partition.nex   IQ-TREE partition: LG+G : AA(1-209),  3di_substmat+G : 3Di(210-418)

Tips whose structure-core length doesn't match the sequence-core length get an all-gap 3Di row
(they stay in the tree via the AA partition; flagged in the report).

Run:  conda run -n pif1 python workflow/13_build_aa3di.py
"""
import argparse
import os
import re
import subprocess
import sys

from Bio import SeqIO


def read_3di(dbprefix):
    ss = [l.rstrip("\n").replace("\x00", "") for l in open(dbprefix + "_ss")]
    ids = [l.split("\t")[1].strip() for l in open(dbprefix + ".lookup")]
    ids = [i[:-4] if i.endswith(".pdb") else i for i in ids]
    return {i: s for i, s in zip(ids, ss) if s}


def kept_columns(trimal, aln_full, expect):
    import tempfile
    tmp = tempfile.NamedTemporaryFile(suffix=".fa", delete=False).name
    out = subprocess.run([trimal, "-in", aln_full, "-automated1", "-out", tmp, "-colnumbering"],
                         capture_output=True, text=True)   # -out sends the aln to file; stdout = colmap only
    if os.path.exists(tmp):
        os.unlink(tmp)
    nums = [int(x) for x in re.findall(r"\d+", out.stdout)]
    if len(nums) != expect:
        sys.stderr.write(f"[13] WARNING: trimAl colmap has {len(nums)} cols, expected {expect}\n")
    return nums


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--aln-full", default="results/seq_tree/aln.fasta")
    ap.add_argument("--aln-trim", default="results/seq_tree/aln.trim.fasta")
    ap.add_argument("--db", default="data/3di/coresdb")
    ap.add_argument("--matrix", default="data/3di/3di_substmat.txt")
    ap.add_argument("--trimal", default="trimal")
    ap.add_argument("--out-3di", default="results/seq_tree/3di.trim.fasta")
    ap.add_argument("--out-concat", default="results/seq_tree/alnAA_3di.fasta")
    ap.add_argument("--out-part", default="results/seq_tree/aa3di_partition.nex")
    ap.add_argument("--aa-model", default="LG+G")
    args = ap.parse_args()

    full = {r.id: str(r.seq) for r in SeqIO.parse(args.aln_full, "fasta")}
    order = [r.id for r in SeqIO.parse(args.aln_trim, "fasta")]
    trim = {r.id: str(r.seq) for r in SeqIO.parse(args.aln_trim, "fasta")}
    ncol = len(trim[order[0]])
    threedi = read_3di(args.db)
    kept = kept_columns(args.trimal, args.aln_full, ncol)
    sys.stderr.write(f"[13] {len(order)} tips, {ncol} AA cols, {len(threedi)} 3Di seqs\n")

    tdi_trim = {}
    n_gapfill = 0
    bad = []
    for tip in order:
        aa = full[tip]
        col2res, ri = {}, 0
        for ci, ch in enumerate(aa):
            if ch != "-":
                col2res[ci] = ri
                ri += 1
        tdi = threedi.get(tip)
        if tdi is None or len(tdi) != ri:
            tdi_trim[tip] = "-" * ncol            # keep tip; 3Di partition contributes nothing
            n_gapfill += 1
            if tdi is not None:
                bad.append(f"{tip}(res{ri}/3di{len(tdi)})")
            continue
        row = ["-" if aa[ci] == "-" else tdi[col2res[ci]] for ci in kept]
        s = "".join(row)
        # gap pattern MUST match the AA trimmed row
        assert [c == "-" for c in s] == [c == "-" for c in trim[tip]], f"gap mismatch {tip}"
        tdi_trim[tip] = s

    with open(args.out_3di, "w") as fh:
        for tip in order:
            fh.write(f">{tip}\n{tdi_trim[tip]}\n")
    with open(args.out_concat, "w") as fh:
        for tip in order:
            fh.write(f">{tip}\n{trim[tip]}{tdi_trim[tip]}\n")

    matrix_abs = os.path.abspath(args.matrix)
    with open(args.out_part, "w") as fh:
        fh.write("#nexus\nbegin sets;\n")
        fh.write(f"  charset AA = 1-{ncol};\n")
        fh.write(f"  charset TDi = {ncol+1}-{2*ncol};\n")
        fh.write(f"  charpartition pif1 = {args.aa_model}:AA, {matrix_abs}+G:TDi;\n")
        fh.write("end;\n")

    sys.stderr.write(f"[13] wrote {args.out_concat} (957 x {2*ncol}) + partition {args.out_part}\n")
    sys.stderr.write(f"[13] gap-filled 3Di for {n_gapfill} tips"
                     + (f" (length mismatches: {bad[:8]})" if bad else " (all length-matched)") + "\n")


if __name__ == "__main__":
    main()
