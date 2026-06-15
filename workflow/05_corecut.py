#!/usr/bin/env python3
"""
05_corecut.py  --  Stage 3: trim every protein to the conserved PIF1 helicase core.

PIF1-family proteins carry long, disordered, non-homologous N/C-terminal extensions
(full-length AFDB mean pLDDT here is only ~60-65; the folded helicase core is much
higher). Aligning full-length sequences lets those termini mislead both the aligner
and the tree. Following the FoldTree paper's "corecut" idea (Fig. 4), we trim each
protein to the region that matches the PIF1 helicase HMM (Pfam PF05970), so MAFFT and
IQ-TREE see only the homologous core. The SAME residue ranges are written to tip_map.tsv
so Stage 4b can slice the *structures* to the matching residues (sequence and structure
trees then describe the same core -> a fair concordance check).

Inputs
------
  data/seqs/selected.faa     (957 sequences, from 04)
  data/seqs/selected.tsv     (metadata: organism, taxid, group/role, gene)
  data/hmm/PF05970.hmm       (PIF1 helicase HMM, fetched from InterPro)

Outputs
-------
  data/seqs/cores.faa              core sequences with tree-safe tip labels
  data/seqs/tip_map.tsv            tip_label, accession, taxid, organism, group, role,
                                   gene, core_from, core_to, full_len, core_len
  results/seq_tree/corecut.domtbl  raw hmmsearch domain table (audit)

Tree-safe tip label = "<accession>_<Genus>" (sanitized; accession guarantees uniqueness).
The taxid in tip_map is the gene->species link GeneRax/NOTUNG need at Stage 6.

Run inside the pif1 env (needs hmmsearch + biopython):
    conda run -n pif1 python workflow/05_corecut.py
"""
import argparse
import csv
import os
import re
import subprocess
import sys
from collections import defaultdict

from Bio import SeqIO


def acc_of(header_or_name):
    """'tr|A0A123|NAME' or 'sp|P07271|PIF1_YEAST' -> 'A0A123' / 'P07271'."""
    parts = header_or_name.split("|")
    return parts[1] if len(parts) >= 2 else header_or_name.split()[0]


def run_hmmsearch(hmm, faa, domtbl, evalue):
    cmd = ["hmmsearch", "--domtblout", domtbl, "-E", str(evalue),
           "--domE", str(evalue), "--noali", hmm, faa]
    sys.stderr.write("[05] " + " ".join(cmd) + "\n")
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL)


def parse_domtbl(path):
    """target-name -> (min env_from, max env_to), 1-based inclusive, over all passing domains."""
    spans = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split()
            if len(f) < 21:
                continue
            target = f[0]
            env_from, env_to = int(f[19]), int(f[20])  # domtblout cols 20,21 (0-based 19,20)
            lo, hi = spans.get(target, (10**9, 0))
            spans[target] = (min(lo, env_from), max(hi, env_to))
    return spans


def safe_label(acc, genus, taxid):
    g = re.sub(r"[^A-Za-z0-9]", "", genus or "")
    if not g:
        g = f"tx{taxid}" if taxid else "sp"
    return f"{acc}_{g}"


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--faa", default="data/seqs/selected.faa")
    ap.add_argument("--tsv", default="data/seqs/selected.tsv")
    ap.add_argument("--hmm", default="data/hmm/PF05970.hmm")
    ap.add_argument("--out-faa", default="data/seqs/cores.faa")
    ap.add_argument("--out-map", default="data/seqs/tip_map.tsv")
    ap.add_argument("--domtbl", default="results/seq_tree/corecut.domtbl")
    ap.add_argument("--evalue", type=float, default=1e-5)
    ap.add_argument("--pad", type=int, default=0, help="residues to pad each side of the envelope")
    args = ap.parse_args()

    os.makedirs(os.path.dirname(args.domtbl), exist_ok=True)

    # metadata by accession
    meta = {}
    with open(args.tsv) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            meta[r["Entry"]] = {
                "taxid": r.get("Organism (ID)", ""), "organism": r.get("Organism", ""),
                "genus": r.get("genus", ""), "group": r.get("group", ""),
                "role": r.get("role", ""),
                "gene": (r.get("Gene Names", "") or "").split()[0] if r.get("Gene Names") else "",
            }

    run_hmmsearch(args.hmm, args.faa, args.domtbl, args.evalue)
    spans = parse_domtbl(args.domtbl)

    n_in = n_out = n_nohit = 0
    seen_labels = set()
    no_hit = []
    with open(args.out_faa, "w") as out_faa, open(args.out_map, "w", newline="") as out_map:
        w = csv.writer(out_map, delimiter="\t")
        w.writerow(["tip_label", "accession", "taxid", "organism", "group", "role",
                    "gene", "core_from", "core_to", "full_len", "core_len"])
        for rec in SeqIO.parse(args.faa, "fasta"):
            n_in += 1
            acc = acc_of(rec.id)
            span = spans.get(rec.id) or spans.get(acc)
            if not span:
                n_nohit += 1
                no_hit.append(acc)
                continue
            lo = max(1, span[0] - args.pad)
            hi = min(len(rec.seq), span[1] + args.pad)
            core = str(rec.seq)[lo - 1:hi]
            m = meta.get(acc, {})
            label = safe_label(acc, m.get("genus", ""), m.get("taxid", ""))
            # guarantee uniqueness (paralogs in the same genus would collide on genus alone)
            base = label
            k = 1
            while label in seen_labels:
                k += 1
                label = f"{base}{k}"
            seen_labels.add(label)
            out_faa.write(f">{label}\n{core}\n")
            w.writerow([label, acc, m.get("taxid", ""), m.get("organism", ""),
                        m.get("group", ""), m.get("role", ""), m.get("gene", ""),
                        lo, hi, len(rec.seq), len(core)])
            n_out += 1

    sys.stderr.write(f"\n[05] {n_in} in -> {n_out} cores -> {args.out_faa}\n")
    sys.stderr.write(f"[05] tip map -> {args.out_map}\n")
    if n_nohit:
        sys.stderr.write(f"[05] WARNING: {n_nohit} sequences had no PF05970 hit "
                         f"(E<{args.evalue}); dropped. First few: {no_hit[:15]}\n")
        sys.stderr.write("[05] (these are cellular-PIF1 by IPR048293 but missed the Pfam "
                         "HMM cutoff -- inspect; may be fragments or divergent.)\n")


if __name__ == "__main__":
    main()
