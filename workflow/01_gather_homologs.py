#!/usr/bin/env python3
"""
01_gather_homologs.py  --  Stage 1 of the PIF1/RRM3 dating pipeline.

Gather every PIF1-family helicase annotated in Ascomycota directly from UniProt,
WITHOUT downloading proteomes or running HMMER locally. UniProt + InterPro have
already annotated these proteins (Pfam PF05970, "PIF1-like helicase"), so a single
paginated API query returns the whole landscape.

What it does
------------
- Queries the UniProt REST API for proteins carrying the PIF1 Pfam domain within
  Ascomycota (NCBI taxon 4890), restricted to reference proteomes by default
  (one canonical proteome per species -> guarantees AlphaFold DB coverage).
- Follows cursor pagination to retrieve ALL hits.
- Parses the taxonomic lineage into subphylum / class / order / family / genus so
  downstream selection can sample evenly across the tree.
- Writes a tab-separated candidate table and prints a breakdown by subphylum/class.

Uses only the Python standard library, so it runs on the base anaconda Python with
zero installs. No API key required; UniProt asks only that you be polite (we page
in chunks and identify ourselves via User-Agent).

Usage
-----
    python workflow/01_gather_homologs.py                 # default: all Ascomycota ref-proteome PF05970
    python workflow/01_gather_homologs.py --out data/seqs/candidates.tsv
    python workflow/01_gather_homologs.py --no-refproteome # include non-reference proteomes too

Seeds for reference: S. cerevisiae PIF1 = P07271, RRM3 = P38766.
Pfam family verified 2026-06-15: PF05970 = "PIF1-like helicase".
"""
import argparse
import csv
import re
import sys
import time
import urllib.parse
import urllib.request

API = "https://rest.uniprot.org/uniprotkb/search"
USER_AGENT = "pif1-foldtree/1.0 (PIF1-RRM3 phylogenetics; contact: spencer.j.gray@gmail.com)"

# UniProt field IDs requested (TSV column order follows this list).
FIELDS = [
    "accession", "id", "protein_name", "gene_names",
    "organism_name", "organism_id", "length",
    "reviewed", "annotation_score", "xref_alphafolddb", "lineage",
]

# Output columns: the requested fields plus the parsed rank columns we add.
RANKS = ["subphylum", "class", "order", "family", "genus"]


def build_query(taxon, xref_term, ref_only):
    # xref_term is e.g. "pfam-PF05970" (Pfam) or "interpro-IPR048293" (InterPro).
    parts = [f"(xref:{xref_term})", f"(taxonomy_id:{taxon})"]
    if ref_only:
        parts.append("(keyword:KW-1185)")  # KW-1185 = Reference proteome
    return " AND ".join(parts)


def fetch_all(query, fields, page_size=500, sleep=0.2):
    """Yield TSV rows (as dicts) across all pages, following the Link: rel=next cursor."""
    params = {
        "query": query,
        "fields": ",".join(fields),
        "format": "tsv",
        "size": str(page_size),
    }
    url = API + "?" + urllib.parse.urlencode(params)
    header = None
    total = None
    n = 0
    while url:
        req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
        with urllib.request.urlopen(req, timeout=60) as resp:
            if total is None:
                total = resp.headers.get("x-total-results")
                if total:
                    sys.stderr.write(f"[01] UniProt reports {total} total hits\n")
            link = resp.headers.get("Link", "")
            body = resp.read().decode("utf-8")
        lines = body.splitlines()
        if header is None:
            header = lines[0].split("\t")
            data_lines = lines[1:]
        else:
            data_lines = lines[1:]  # every page repeats the header row
        for line in data_lines:
            cols = line.split("\t")
            n += 1
            yield dict(zip(header, cols))
        # find next-page cursor
        m = re.search(r'<([^>]+)>;\s*rel="next"', link)
        url = m.group(1) if m else None
        if url:
            time.sleep(sleep)
    sys.stderr.write(f"[01] retrieved {n} rows\n")


def parse_lineage(lineage):
    """Extract named ranks from a UniProt 'Taxonomic lineage' string.

    Tokens look like 'Saccharomycotina (subphylum)'. Returns {rank: name}.
    """
    out = {r: "" for r in RANKS}
    for token in lineage.split(","):
        token = token.strip()
        m = re.match(r"^(.*)\s+\(([^)]+)\)$", token)
        if not m:
            continue
        name, rank = m.group(1).strip(), m.group(2).strip()
        if rank in out and not out[rank]:
            out[rank] = name
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--taxon", default="4890", help="NCBI taxon id (default 4890 = Ascomycota)")
    ap.add_argument("--pfam", default="PF05970", help="Pfam family (default PF05970 = PIF1)")
    ap.add_argument("--interpro", default=None,
                    help="InterPro id to use INSTEAD of --pfam, e.g. IPR048293 "
                         "('PIF1_RRM3_pfh1' = cellular PIF1/RRM3/Pfh1; excludes Helitron "
                         "transposon helicases that also carry the bare PF05970 domain)")
    ap.add_argument("--out", default="data/seqs/candidates.tsv", help="output TSV path")
    ap.add_argument("--no-refproteome", dest="ref_only", action="store_false",
                    help="do not restrict to reference proteomes")
    ap.add_argument("--page-size", type=int, default=500)
    args = ap.parse_args()

    xref_term = f"interpro-{args.interpro}" if args.interpro else f"pfam-{args.pfam}"
    query = build_query(args.taxon, xref_term, args.ref_only)
    sys.stderr.write(f"[01] query: {query}\n")

    out_cols = FIELDS + RANKS
    rows = []
    by_subphylum, by_class = {}, {}
    species = set()
    for rec in fetch_all(query, FIELDS, page_size=args.page_size):
        ranks = parse_lineage(rec.get("Taxonomic lineage", ""))
        # UniProt TSV header names differ from field IDs; map back by position-safe keys
        merged = dict(rec)
        merged.update(ranks)
        rows.append(merged)
        sp = rec.get("Organism (ID)", "")
        species.add(sp)
        by_subphylum[ranks["subphylum"]] = by_subphylum.get(ranks["subphylum"], 0) + 1
        by_class[ranks["class"]] = by_class.get(ranks["class"], 0) + 1

    # The TSV header names from UniProt -> our normalized column keys.
    # We write the UniProt header verbatim for the requested fields, then the rank cols.
    if rows:
        uniprot_headers = [h for h in rows[0].keys() if h not in RANKS]
        header_row = uniprot_headers + RANKS
        with open(args.out, "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(header_row)
            for r in rows:
                w.writerow([r.get(h, "") for h in header_row])

    sys.stderr.write(f"\n[01] wrote {len(rows)} proteins from {len(species)} species -> {args.out}\n")
    sys.stderr.write("\n[01] proteins per subphylum:\n")
    for k, v in sorted(by_subphylum.items(), key=lambda kv: -kv[1]):
        sys.stderr.write(f"    {v:5d}  {k or '(unparsed)'}\n")
    sys.stderr.write("\n[01] proteins per class (top 20):\n")
    for k, v in sorted(by_class.items(), key=lambda kv: -kv[1])[:20]:
        sys.stderr.write(f"    {v:5d}  {k or '(unparsed)'}\n")


if __name__ == "__main__":
    main()
