#!/usr/bin/env python3
"""
21_recq_select.py  --  RecQ Phase-1 step: subsample the tree-of-life RecQ landscape to a
phylogenetically-even, anchor-inclusive representative set.

WHY. Unlike the fungi-bounded PIF1 family (957 proteins), RecQ is in nearly every cellular
organism, so gathering IPR004589 across reference proteomes returns ~37.7k proteins
(euk 14026, bac 23500, arc 132) whose counts track SEQUENCING DENSITY (vertebrates, plants,
Actinobacteria), not phylogenetic breadth. Using all of them would (a) bias the conservation
estimate toward over-sequenced lineages and (b) be intractable for AF3 (30 jobs/day). So we take
ONE best representative per taxonomic family across all three domains -> an even, de-biased set on
the PIF1 scale, with the experimentally-anchored members forced in.

WHAT.
  - Reads the three gather TSVs (recq_euk/arc/bac.tsv), tags each row with its domain.
  - Groups by --rank (default: family; falls back order->class if a row lacks family) and keeps the
    best --max-per-group reps per group, ranked by quality: reviewed > annotation score > has-AFDB
    > length (prefer curated, well-annotated, structurally-covered, reasonably complete).
  - ALWAYS keeps the anchors (RECQL1/4/5, BLM, WRN, Sgs1, E. coli RecQ) regardless of grouping.
  - Writes data/recq/recq_selected.tsv (gather columns + 'domain' + 'anchor') and prints the
    domain / class breakdown so the spread is auditable.

Next: 04_fetch_fasta.py on recq_selected.tsv -> 05_corecut.py (RQC/helicase) -> align -> b-wing conservation.

Run:  python workflow/21_recq_select.py            # 1 per family, all domains, anchors forced
      python workflow/21_recq_select.py --rank order   # tighter ~660-tip backbone
"""
import argparse
import csv
from collections import Counter, defaultdict

ANCHORS = {  # accession -> label (experimentally characterised RecQ + the b-wing structures)
    "P46063": "RECQL1", "O94761": "RECQL4", "O94762": "RECQL5", "P54132": "BLM",
    "Q14191": "WRN", "P35187": "Sgs1", "P15043": "EcRecQ",
}
FILES = {"euk": "data/recq/recq_euk.tsv", "arc": "data/recq/recq_arc.tsv", "bac": "data/recq/recq_bac.tsv"}


def quality(r):
    """Higher = better representative for a group."""
    reviewed = 1 if r.get("Reviewed", "").lower() == "reviewed" else 0
    try:
        ann = float(r.get("Annotation", "") or 0)
    except ValueError:
        ann = 0.0
    afdb = 1 if r.get("AlphaFoldDB", "").strip() else 0
    try:
        length = int(r.get("Length", "") or 0)
    except ValueError:
        length = 0
    return (reviewed, ann, afdb, length)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--rank", default="family", choices=["family", "order", "class", "genus"],
                    help="taxonomic level to even out across (default: family)")
    ap.add_argument("--max-per-group", type=int, default=1, help="reps kept per group (default 1)")
    ap.add_argument("--out", default="data/recq/recq_selected.tsv")
    args = ap.parse_args()

    rows = []
    header = None
    for dom, f in FILES.items():
        with open(f) as fh:
            rdr = csv.DictReader(fh, delimiter="\t")
            header = rdr.fieldnames
            for r in rdr:
                r["domain"] = dom
                rows.append(r)
    acc_col = header[0]  # 'Entry'

    # group key: chosen rank, falling back to coarser ranks so nothing is silently dropped
    def gkey(r):
        for rk in (args.rank, "order", "class", "genus"):
            v = r.get(rk, "").strip()
            if v:
                return f"{r['domain']}:{rk}:{v}"
        return f"{r['domain']}:NA"

    groups = defaultdict(list)
    anchors_rows = {}
    for r in rows:
        if r.get(acc_col, "") in ANCHORS:
            anchors_rows[r[acc_col]] = r
        groups[gkey(r)].append(r)

    selected, seen = [], set()
    for r in anchors_rows.values():          # anchors first, always kept
        r["anchor"] = ANCHORS[r[acc_col]]
        selected.append(r); seen.add(r[acc_col])
    for k in sorted(groups):                 # one (or N) best per group
        reps = sorted(groups[k], key=quality, reverse=True)
        kept = 0
        for r in reps:
            if kept >= args.max_per_group:
                break
            if r[acc_col] in seen:
                continue
            r.setdefault("anchor", "")
            selected.append(r); seen.add(r[acc_col]); kept += 1

    for r in selected:
        r["group"] = r["domain"]   # 'group' is the clade key the downstream 05/18 scripts read
    out_cols = header + ["domain", "group", "anchor"]
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=out_cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in selected:
            w.writerow(r)

    by_dom = Counter(r["domain"] for r in selected)
    by_class = Counter(r.get("class", "") or "(unparsed)" for r in selected)
    print(f"[21] selected {len(selected)} RecQ reps (rank={args.rank}, <= {args.max_per_group}/group) -> {args.out}")
    print(f"[21] by domain: " + ", ".join(f"{d}:{n}" for d, n in by_dom.most_common()))
    print(f"[21] anchors kept: {len(anchors_rows)}/{len(ANCHORS)} (" +
          ", ".join(ANCHORS[a] for a in anchors_rows) + ")")
    print("[21] top classes:")
    for c, n in by_class.most_common(15):
        print(f"     {n:4d}  {c}")


if __name__ == "__main__":
    main()
