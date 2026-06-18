#!/usr/bin/env python3
"""
14_prep_generax.py  --  build ready-to-run GeneRax inputs for a given gene tree.

GeneRax (the gold-standard ML duplication-loss reconciliation) segfaults under Rosetta on the Mac,
so this just PREPARES the inputs; the actual `generax` run happens on native x86 (WSL2 on the
Windows desktop -- see WINDOWS_HANDOFF.md). GeneRax needs every gene tip mapped to a species-tree
leaf, so we prune the gene tree + alignment to the 941 genes that map to NCBI species (16 genes on
8 taxa NCBI doesn't recognize are dropped), and write the GeneRax families file.

For each gene tree it writes (into <outdir>/<label>/):
  genetree.pruned.nwk   gene tree pruned to mapped taxa
  aln.pruned.fasta      aln.trim.fasta pruned to the same taxa (GeneRax likelihood term)
  families.txt          GeneRax families file (UndatedDL, LG+G)

Run:  conda run -n pif1 python workflow/14_prep_generax.py \
          --genetree results/seq_tree/pif1_aa3di.treefile --label aa3di
"""
import argparse
import csv
import os
import sys

import dendropy
from Bio import SeqIO


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genetree", required=True)
    ap.add_argument("--label", required=True, help="subfolder name, e.g. aa3di / lgig")
    ap.add_argument("--aln", default="results/seq_tree/aln.trim.fasta")
    ap.add_argument("--map", default="data/species_tree/gene_species.map")
    ap.add_argument("--species-tree", default="data/species_tree/ncbi_species.nwk")
    ap.add_argument("--outdir", default="results/reconciliation")
    args = ap.parse_args()

    mapped = {}
    for line in open(args.map):
        p = line.split()
        if len(p) >= 2:
            mapped[p[0]] = p[1]
    keep = set(mapped)

    outdir = os.path.join(args.outdir, args.label)
    os.makedirs(outdir, exist_ok=True)

    # prune gene tree to mapped taxa
    t = dendropy.Tree.get(path=args.genetree, schema="newick", preserve_underscores=True)
    present = {l.taxon.label for l in t.leaf_node_iter()}
    t.retain_taxa_with_labels([x for x in present if x in keep])
    tree_path = os.path.join(outdir, "genetree.pruned.nwk")
    t.write(path=tree_path, schema="newick", unquoted_underscores=True,
            suppress_rooting=True)
    kept_tips = {l.taxon.label for l in t.leaf_node_iter()}

    # prune alignment to the same taxa
    aln_path = os.path.join(outdir, "aln.pruned.fasta")
    n_aln = 0
    with open(aln_path, "w") as out:
        for r in SeqIO.parse(args.aln, "fasta"):
            if r.id in kept_tips:
                out.write(f">{r.id}\n{r.seq}\n")
                n_aln += 1

    # GeneRax families file. mapping = the project gene_species.map (gene<space>species, 941 lines)
    fam_path = os.path.join(outdir, "families.txt")
    with open(fam_path, "w") as fh:
        # RELATIVE paths (portable across machines) -- run generax from the repo root
        fh.write("[FAMILIES]\n- pif1\n")
        fh.write(f"starting_gene_tree = {tree_path}\n")
        fh.write(f"alignment = {aln_path}\n")
        fh.write(f"mapping = {args.map}\n")
        fh.write("subst_model = LG+G\n")

    sys.stderr.write(f"[14] {args.label}: pruned tree {len(kept_tips)} tips, aln {n_aln} seqs "
                     f"(dropped {len(present)-len(kept_tips)} unmapped) -> {outdir}\n")
    sys.stderr.write(f"[14] GeneRax (on WSL2):\n"
                     f"     generax -f {fam_path} -s {args.species_tree} \\\n"
                     f"        -r UndatedDL --unrooted-gene-tree --strategy SPR -p {outdir}/run\n")


if __name__ == "__main__":
    main()
