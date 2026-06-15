#!/usr/bin/env python3
"""
04_fetch_fasta.py  --  fetch protein sequences for the selected set.

Reads the master list (data/seqs/selected.tsv), pulls each protein's amino-acid
sequence from UniProt in batches via the REST 'stream' endpoint, and writes a single
FASTA (data/seqs/selected.faa). This FASTA feeds Stage 3 (corecut to the helicase
domain) and Stage 4a (MAFFT alignment + IQ-TREE). Verifies every requested accession
came back. Stdlib only -- runs on base Python.

Tree-tip naming is handled later (at the alignment step): UniProt FASTA headers look
like '>sp|P07271|PIF1_YEAST ...' / '>tr|<acc>|...'; we keep them verbatim here and
sanitize to safe tip labels (accession + organism) during corecut/alignment.

Usage
-----
    python workflow/04_fetch_fasta.py --in data/seqs/selected.tsv --out data/seqs/selected.faa
"""
import argparse
import csv
import sys
import time
import urllib.parse
import urllib.request

STREAM = "https://rest.uniprot.org/uniprotkb/stream"
USER_AGENT = "pif1-foldtree/1.0 (PIF1-RRM3 phylogenetics; contact: spencer.j.gray@gmail.com)"


def chunks(lst, n):
    for i in range(0, len(lst), n):
        yield lst[i:i + n]


def fetch_fasta(accs, timeout=120):
    query = "accession:(" + " OR ".join(accs) + ")"
    url = STREAM + "?" + urllib.parse.urlencode({"query": query, "format": "fasta"})
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.read().decode("utf-8")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--in", dest="infile", default="data/seqs/selected.tsv")
    ap.add_argument("--out", default="data/seqs/selected.faa")
    ap.add_argument("--batch", type=int, default=100, help="accessions per request")
    ap.add_argument("--sleep", type=float, default=0.3)
    args = ap.parse_args()

    with open(args.infile) as fh:
        accs = [r["Entry"] for r in csv.DictReader(fh, delimiter="\t") if r.get("Entry")]
    sys.stderr.write(f"[04] requesting {len(accs)} sequences in batches of {args.batch}\n")

    parts = []
    for i, batch in enumerate(chunks(accs, args.batch), 1):
        fa = fetch_fasta(batch)
        parts.append(fa)
        sys.stderr.write(f"  batch {i}: got {fa.count('>')} (asked {len(batch)})\n")
        time.sleep(args.sleep)

    text = "".join(parts)
    with open(args.out, "w") as out:
        out.write(text)

    got = text.count(">")
    sys.stderr.write(f"[04] wrote {got} sequences -> {args.out} (requested {len(accs)})\n")
    if got != len(accs):
        # identify which accessions are missing so they can be chased by hand
        present = set()
        for line in text.splitlines():
            if line.startswith(">"):
                # header like >sp|P07271|... or >tr|A0A...|...
                f = line.split("|")
                if len(f) >= 2:
                    present.add(f[1])
        missing = [a for a in accs if a not in present]
        sys.stderr.write(f"[04] WARNING: {len(missing)} missing: {missing[:20]}"
                         f"{' ...' if len(missing) > 20 else ''}\n")


if __name__ == "__main__":
    main()
