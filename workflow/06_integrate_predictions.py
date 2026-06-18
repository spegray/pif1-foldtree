#!/usr/bin/env python3
"""
06_integrate_predictions.py  --  Stage 2b: fold the ColabFold gap predictions into the manifest.

After Stage 2b (localcolabfold on the RTX 4090) the 129 AFDB-absent proteins exist as
data/structures/af3/<ACCESSION>.pdb. This script promotes their manifest.csv rows from
"MISSING_predict_with_AF3" to real structures so the structural pipeline (corecut -> FoldTree,
Stage 4b) treats them exactly like the 828 AFDB entries.

For each gap row it:
  - locates data/structures/af3/<accession>.pdb,
  - reads the per-residue pLDDT from the B-factor column (mean over CA atoms),
  - sets structure_source = AF2_ColabFold, structure_path, mean_plddt,
  - flags the 3 core-only giants (from af3/_cores_note.txt) in a new `structure_note`
    column WITH their original-protein envelope coords -- this is critical for Stage 3:
    those PDBs are already the helicase core renumbered from 1, so the structural corecut
    must NOT re-slice them by full-length coordinates (it would slice the wrong residues).

Stdlib only (matches 03/04); re-runnable. Reads + rewrites manifest.csv in place
(it is git-tracked, so the previous version is always recoverable via git).

Usage:
    python3 workflow/06_integrate_predictions.py
    python3 workflow/06_integrate_predictions.py --manifest manifest.csv --af3dir data/structures/af3
"""
import argparse
import csv
import os
import sys


def mean_ca_plddt(pdb_path):
    """Mean of the B-factor (= AF pLDDT) over CA atoms = mean per-residue pLDDT."""
    vals = []
    with open(pdb_path) as fh:
        for line in fh:
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                try:
                    vals.append(float(line[60:66]))
                except ValueError:
                    pass
    return (sum(vals) / len(vals)) if vals else None


def load_core_only(note_path):
    """af3/_cores_note.txt -> {accession: 'core_only ... env a-b (full N aa)'} for the giants."""
    out = {}
    if not os.path.exists(note_path):
        return out
    with open(note_path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            acc = (r.get("accession") or "").strip()
            if acc:
                out[acc] = (f"core_only_AF2 (env {r.get('env_from')}-{r.get('env_to')} of "
                            f"{r.get('full_len')}aa; AF3 full-length pending)")
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--manifest", default="manifest.csv")
    ap.add_argument("--af3dir", default="data/structures/af3")
    ap.add_argument("--source-tag", default="AF2_ColabFold")
    args = ap.parse_args()

    core_only = load_core_only(os.path.join(args.af3dir, "_cores_note.txt"))

    with open(args.manifest) as fh:
        reader = csv.DictReader(fh)
        fields = list(reader.fieldnames)
        rows = list(reader)
    if "structure_note" not in fields:
        fields.append("structure_note")

    n_updated = n_core = 0
    still_missing = []
    plddts = []
    for row in rows:
        is_gap = ("MISSING" in (row.get("structure_source") or "")) or not (row.get("structure_path") or "")
        if not is_gap:
            continue
        acc = row["accession"]
        pdb = os.path.join(args.af3dir, f"{acc}.pdb")
        if not os.path.exists(pdb):
            still_missing.append(acc)
            continue
        plddt = mean_ca_plddt(pdb)
        row["structure_source"] = args.source_tag
        row["structure_path"] = pdb
        row["afdb_version"] = ""
        row["mean_plddt"] = f"{plddt:.2f}" if plddt is not None else ""
        if acc in core_only:
            row["structure_note"] = core_only[acc]
            n_core += 1
        else:
            row.setdefault("structure_note", "")
        n_updated += 1
        if plddt is not None:
            plddts.append(plddt)

    # ensure every row has the new column key so DictWriter doesn't choke
    for row in rows:
        row.setdefault("structure_note", "")

    with open(args.manifest, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)

    # ---- report ----
    sys.stderr.write(f"\n[06] integrated {n_updated} ColabFold structures into {args.manifest}\n")
    sys.stderr.write(f"[06]   of which {n_core} are core-only giants (flagged in structure_note)\n")
    if plddts:
        plddts.sort()
        med = plddts[len(plddts) // 2]
        sys.stderr.write(f"[06]   pLDDT over the new set: min {min(plddts):.1f} / "
                         f"median {med:.1f} / max {max(plddts):.1f}\n")
    if still_missing:
        sys.stderr.write(f"[06] WARNING: {len(still_missing)} gap rows still have NO pdb: "
                         f"{still_missing[:15]}\n")
    else:
        sys.stderr.write("[06] all gaps resolved -- 0 rows remain MISSING.\n")

    # final manifest sanity: source distribution
    from collections import Counter
    dist = Counter(r.get("structure_source", "") for r in rows)
    sys.stderr.write("[06] manifest structure_source now: "
                     + ", ".join(f"{k}={v}" for k, v in sorted(dist.items())) + "\n")


if __name__ == "__main__":
    main()
