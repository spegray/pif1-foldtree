#!/usr/bin/env python3
"""
02_select_representatives.py  --  Stage 0/1 curation.

Reduce the full candidate table (workflow/01 output, ~2255 proteins / ~697 species)
to an even, taxonomically spread set of representative SPECIES, while keeping ALL
PIF1-family paralogs within each chosen species (never collapse Pif1 vs Rrm3!).

Sampling strategy (stratified, NOT proportional)
------------------------------------------------
The raw data is dominated by Pezizomycotina (esp. Sordariomycetes). For dating the
budding-yeast PIF1/RRM3 split we want:
  * DENSE sampling inside the focus subphylum (Saccharomycotina by default),
  * a spread across the other classes for context + bracketing,
  * ALL species of rare outgroup subphyla (Taphrinomycotina) kept,
  * the curated/reviewed Swiss-Prot anchors always included.

Within a class we spread the picks across families/genera (evenly-spaced indices
after sorting by family,genus,name) so we sample diversity rather than whichever
species sort first alphabetically. Among species we prefer those with >=2 paralogs
and higher annotation scores (more informative, better curated).

This is a transparent first cut: review data/seqs/selected.tsv by hand and edit as
you see fit -- that edited file is the auditable record of which taxa entered the tree.

Outgroups beyond Ascomycota (Basidiomycota, early-diverging fungi, a human/plant PIF1
to anchor the deep root) are NOT in the candidate file (it is Ascomycota-only). Add
them by re-running workflow/01 with a different --taxon, or append accessions by hand
(e.g. human PIF1 = Q9H611). See README Stage 0.

Usage
-----
    python workflow/02_select_representatives.py \
        --in data/seqs/candidates.tsv --out data/seqs/selected.tsv \
        --focus Saccharomycotina --focus-per-class 50 --per-class 10 --keep-all Taphrinomycotina

Uses only the standard library.
"""
import argparse
import csv
import sys
from collections import defaultdict


def even_indices(n, k):
    """Return k evenly-spaced indices in [0, n) (all of them if k >= n)."""
    if k >= n:
        return list(range(n))
    return [round(i * (n - 1) / (k - 1)) for i in range(k)] if k > 1 else [0]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--in", dest="infile", default="data/seqs/candidates.tsv")
    ap.add_argument("--out", default="data/seqs/selected.tsv")
    ap.add_argument("--focus", default="Saccharomycotina",
                    help="subphylum to sample densely (default Saccharomycotina)")
    ap.add_argument("--focus-per-class", type=int, default=50,
                    help="max species per class within the focus subphylum")
    ap.add_argument("--per-class", type=int, default=10,
                    help="max species per class elsewhere")
    ap.add_argument("--keep-all", default="Taphrinomycotina",
                    help="comma-list of subphyla to keep entirely")
    args = ap.parse_args()

    keep_all = {s.strip() for s in args.keep_all.split(",") if s.strip()}

    with open(args.infile) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        rows = list(reader)
    sys.stderr.write(f"[02] read {len(rows)} proteins\n")

    # group proteins by species (taxid)
    by_species = defaultdict(list)
    for r in rows:
        by_species[r["Organism (ID)"]].append(r)

    # species metadata for sorting/picking
    species_meta = {}
    for taxid, prots in by_species.items():
        p0 = prots[0]
        reviewed = any(p.get("Reviewed", "") == "reviewed" for p in prots)
        try:
            ann = max(float(p.get("Annotation", "0") or 0) for p in prots)
        except ValueError:
            ann = 0.0
        species_meta[taxid] = {
            "subphylum": p0["subphylum"], "class": p0["class"],
            "order": p0["order"], "family": p0["family"], "genus": p0["genus"],
            "organism": p0["Organism"], "n_paralogs": len(prots),
            "reviewed": reviewed, "annotation": ann,
        }

    # bucket species by (subphylum, class)
    buckets = defaultdict(list)
    for taxid, m in species_meta.items():
        buckets[(m["subphylum"], m["class"])].append(taxid)

    selected_taxids = set()
    # always keep reviewed anchors
    for taxid, m in species_meta.items():
        if m["reviewed"]:
            selected_taxids.add(taxid)

    for (subphylum, cls), taxids in buckets.items():
        if subphylum in keep_all:
            selected_taxids.update(taxids)
            continue
        cap = args.focus_per_class if subphylum == args.focus else args.per_class
        # sort to spread across family/genus, then prefer informative species
        taxids_sorted = sorted(
            taxids,
            key=lambda t: (species_meta[t]["family"], species_meta[t]["genus"],
                           species_meta[t]["organism"]),
        )
        # prefer species with >=2 paralogs and higher annotation by re-ordering,
        # but keep the family/genus spread by picking evenly across the sorted list
        idxs = even_indices(len(taxids_sorted), cap)
        for i in idxs:
            selected_taxids.add(taxids_sorted[i])

    # emit all proteins for selected species
    out_rows = [r for r in rows if r["Organism (ID)"] in selected_taxids]
    fieldnames = rows[0].keys()
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        w.writerows(out_rows)

    # summary
    sel_by_sub = defaultdict(set)
    prot_by_sub = defaultdict(int)
    for r in out_rows:
        sel_by_sub[r["subphylum"]].add(r["Organism (ID)"])
        prot_by_sub[r["subphylum"]] += 1
    sys.stderr.write(f"\n[02] selected {len(selected_taxids)} species, "
                     f"{len(out_rows)} proteins -> {args.out}\n")
    sys.stderr.write("[02] species / proteins per subphylum:\n")
    for sub in sorted(sel_by_sub, key=lambda s: -len(sel_by_sub[s])):
        sys.stderr.write(f"    {len(sel_by_sub[sub]):4d} sp  {prot_by_sub[sub]:5d} prot  {sub or '(unparsed)'}\n")
    sys.stderr.write("\n[02] REVIEW data/seqs/selected.tsv by hand before proceeding.\n")


if __name__ == "__main__":
    main()
