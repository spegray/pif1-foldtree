#!/usr/bin/env python3
"""
03_fetch_structures.py  --  Stage 2 of the PIF1/RRM3 dating pipeline.

Download AlphaFold DB structures for the selected proteins and build the project
manifest. Hits the AFDB *prediction API* per accession rather than guessing a file
URL, because the AFDB version string changes over time (verified 2026-06-15: current
version is v6, NOT the v4 quoted in older docs). The API also returns provenance we
record for auditability: model version, mean pLDDT, organism, taxid, and whether the
entry is a reference proteome.

Outputs
-------
- data/structures/afdb/AF-<acc>-F1-model_v<ver>.<fmt>   (one file per protein)
- manifest.csv  (one row per protein: accession, gene, organism, taxid, subphylum,
  class, length, structure source, AFDB version, mean pLDDT, file path, paralog label)

Proteins with NO AFDB entry are written to manifest.csv with source="MISSING_predict_with_AF3"
so you know exactly which handful need AlphaFold 3 (ColabFold / AF3 server) -- see README Stage 2.

Usage
-----
    python workflow/03_fetch_structures.py --in data/seqs/selected.tsv \
        --manifest manifest.csv --fmt cif
    python workflow/03_fetch_structures.py --in data/seqs/selected.tsv --limit 6   # quick test

Uses only the standard library.
"""
import argparse
import csv
import json
import os
import sys
import time
import urllib.error
import urllib.request

AFDB_API = "https://alphafold.ebi.ac.uk/api/prediction/"
USER_AGENT = "pif1-foldtree/1.0 (PIF1-RRM3 phylogenetics; contact: spencer.j.gray@gmail.com)"

MANIFEST_COLS = [
    "accession", "gene", "organism", "taxid", "subphylum", "class", "order",
    "family", "genus", "length", "reviewed", "structure_source", "afdb_version",
    "mean_plddt", "is_reference_proteome", "structure_path", "paralog_label",
]


def get_json(url, timeout=60):
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return json.load(resp)


def download(url, dest, timeout=120):
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(req, timeout=timeout) as resp, open(dest, "wb") as out:
        out.write(resp.read())


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--in", dest="infile", default="data/seqs/selected.tsv")
    ap.add_argument("--manifest", default="manifest.csv")
    ap.add_argument("--outdir", default="data/structures/afdb")
    ap.add_argument("--fmt", choices=["cif", "pdb"], default="cif",
                    help="structure format to download (cif recommended; foldseek reads both)")
    ap.add_argument("--limit", type=int, default=0, help="process only first N (testing)")
    ap.add_argument("--sleep", type=float, default=0.3, help="seconds between requests (politeness)")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    with open(args.infile) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    if args.limit:
        rows = rows[:args.limit]
    sys.stderr.write(f"[03] processing {len(rows)} proteins\n")

    manifest = []
    n_ok = n_missing = 0
    for i, r in enumerate(rows, 1):
        acc = r["Entry"]
        gene = (r.get("Gene Names", "") or "").split()[0] if r.get("Gene Names") else ""
        base = {
            "accession": acc, "gene": gene, "organism": r.get("Organism", ""),
            "taxid": r.get("Organism (ID)", ""), "subphylum": r.get("subphylum", ""),
            "class": r.get("class", ""), "order": r.get("order", ""),
            "family": r.get("family", ""), "genus": r.get("genus", ""),
            "length": r.get("Length", ""), "reviewed": r.get("Reviewed", ""),
            "paralog_label": "",  # filled later from the gene tree, NOT assumed here
        }
        try:
            data = get_json(AFDB_API + acc)
            entry = data[0] if isinstance(data, list) else data
            ver = entry.get("latestVersion", "")
            url = entry.get("cifUrl") if args.fmt == "cif" else entry.get("pdbUrl")
            fname = os.path.basename(url)
            dest = os.path.join(args.outdir, fname)
            if not os.path.exists(dest):
                download(url, dest)
            base.update({
                "structure_source": "AFDB", "afdb_version": ver,
                "mean_plddt": entry.get("globalMetricValue", ""),
                "is_reference_proteome": entry.get("isUniProtReferenceProteome", ""),
                "structure_path": dest,
            })
            n_ok += 1
            sys.stderr.write(f"  [{i}/{len(rows)}] {acc}  v{ver}  pLDDT={entry.get('globalMetricValue')}  OK\n")
        except urllib.error.HTTPError as e:
            if e.code == 404:
                base.update({"structure_source": "MISSING_predict_with_AF3",
                             "afdb_version": "", "mean_plddt": "",
                             "is_reference_proteome": "", "structure_path": ""})
                n_missing += 1
                sys.stderr.write(f"  [{i}/{len(rows)}] {acc}  NO AFDB ENTRY -> needs AF3\n")
            else:
                raise
        manifest.append(base)
        time.sleep(args.sleep)

    with open(args.manifest, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=MANIFEST_COLS)
        w.writeheader()
        w.writerows(manifest)

    sys.stderr.write(f"\n[03] {n_ok} structures downloaded, {n_missing} missing (need AF3) "
                     f"-> {args.manifest}\n")


if __name__ == "__main__":
    main()
