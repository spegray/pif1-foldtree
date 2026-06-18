#!/usr/bin/env python3
"""
07_structcorecut.py  --  Stage 3 (structures): slice every structure to the helicase core.

The sequence tree (Stage 4a) was built from cores trimmed to the PF05970 HMM envelope
(workflow/05 -> data/seqs/tip_map.tsv, columns core_from/core_to, 1-based inclusive, in
FULL-LENGTH UniProt numbering). To make the sequence-vs-structure concordance check fair
(QC #3 in the README), the STRUCTURAL tree must describe the SAME residues. So here we slice
each structure to the very same core_from..core_to span and write one .pdb per protein, named
by the SAME tree-safe tip_label, so the FoldTree tips line up 1:1 with the IQ-TREE tips.

Inputs
------
  data/seqs/tip_map.tsv     tip_label, accession, ... core_from, core_to   (defines the tip set)
  manifest.csv              accession -> structure_path (.cif AFDB / .pdb ColabFold), structure_note

Output
------
  data/structures/cores/<tip_label>.pdb     core-only structures for foldseek (Stage 4b)

Special case -- the 3 core-only giants (manifest structure_note ~ "core_only"): their ColabFold
.pdb is ALREADY just the helicase core, renumbered from residue 1, so core_from/core_to (which are
full-length coords like 166-369) do NOT exist in the file. For these we keep the WHOLE structure.

Foldseek reads .cif and .pdb; we normalize everything to .pdb cores for a uniform input folder.

Run in the pif1 env (needs biopython):
    conda run -n pif1 python workflow/07_structcorecut.py
    conda run -n pif1 python workflow/07_structcorecut.py --limit 3   # quick test
"""
import argparse
import csv
import os
import sys
import warnings

from Bio.PDB import MMCIFParser, PDBParser, PDBIO, Select
from Bio import BiopythonWarning

warnings.simplefilter("ignore", BiopythonWarning)


class CoreSelect(Select):
    """Keep standard polymer residues whose author seq number is within [lo, hi]."""
    def __init__(self, lo, hi):
        self.lo, self.hi = lo, hi

    def accept_residue(self, residue):
        hetflag, resseq, icode = residue.id
        return 1 if (hetflag == " " and self.lo <= resseq <= self.hi) else 0


def load_structure(path):
    ext = os.path.splitext(path)[1].lower()
    parser = MMCIFParser(QUIET=True) if ext in (".cif", ".mmcif") else PDBParser(QUIET=True)
    return parser.get_structure("s", path)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tip-map", default="data/seqs/tip_map.tsv")
    ap.add_argument("--manifest", default="manifest.csv")
    ap.add_argument("--outdir", default="data/structures/cores")
    ap.add_argument("--limit", type=int, default=0)
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    # accession -> structure_path, structure_note
    struct = {}
    with open(args.manifest) as fh:
        for r in csv.DictReader(fh):
            struct[r["accession"]] = (r.get("structure_path", ""),
                                      (r.get("structure_note", "") or ""))

    rows = list(csv.DictReader(open(args.tip_map), delimiter="\t"))
    if args.limit:
        rows = rows[:args.limit]
    sys.stderr.write(f"[07] slicing {len(rows)} structures to their helicase core\n")

    io = PDBIO()
    n_ok = n_core_only = 0
    missing, failed = [], []
    for i, r in enumerate(rows, 1):
        tip, acc = r["tip_label"], r["accession"]
        path, note = struct.get(acc, ("", ""))
        if not path or not os.path.exists(path):
            missing.append(acc)
            continue
        try:
            s = load_structure(path)
            if "core_only" in note:
                lo, hi = -10**9, 10**9          # already the core -> keep whole structure
                n_core_only += 1
            else:
                lo, hi = int(r["core_from"]), int(r["core_to"])
            io.set_structure(s)
            io.save(os.path.join(args.outdir, f"{tip}.pdb"), CoreSelect(lo, hi))
            n_ok += 1
        except Exception as e:                  # noqa: BLE001 -- log & continue, don't abort batch
            failed.append(f"{acc}({type(e).__name__})")
        if i % 200 == 0:
            sys.stderr.write(f"  [{i}/{len(rows)}] ...\n")

    sys.stderr.write(f"\n[07] wrote {n_ok} core .pdb -> {args.outdir} "
                     f"({n_core_only} kept whole = core-only giants)\n")
    if missing:
        sys.stderr.write(f"[07] {len(missing)} had no structure on disk: {missing[:10]}\n")
    if failed:
        sys.stderr.write(f"[07] {len(failed)} failed to parse/slice: {failed[:10]}\n")


if __name__ == "__main__":
    main()
