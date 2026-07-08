#!/usr/bin/env python3
"""
16_graft_species_tree.py  --  Graft the Shen 2018 budding-yeast topology into the NCBI species-tree
backbone, producing a phylogenomic species tree for the M5/handoff-R4 GeneRax reconciliation.

WHY. Reconciliation in this project used an NCBI-taxonomy species tree. The headline placement is
already shown robust to a real phylogeny for the DATING (workflow/17, on Shen 2018). This script
provides the remaining piece: a species tree whose Saccharomycotina portion follows the Shen 2018
phylogenomic topology, so GeneRax can be re-run and the duplication/loss pattern confirmed on a
phylogeny rather than on NCBI's polytomy resolution.

KEY CONSTRAINT. The species-tree leaves are NCBI **taxids** (gene_species.map keys on taxids, so
GeneRax's gene->species mapping does too). Shen tips are `Genus_species`. The graft therefore relabels
Shen tips back to the corresponding taxids and PRESERVES the exact taxid leaf set (verified at the end).

WHAT IT DOES.
  1. taxid -> normalized `genus_species` via data/seqs/tip_map.tsv (same norm() as workflow/17).
  2. matched = NCBI taxid-leaves whose species is present in the Shen tree (the budding yeasts).
  3. NCBI MRCA(matched) = the clade to replace. Prune Shen to the matched species, relabel its tips
     Genus_species -> taxid, re-parse it into the NCBI taxon namespace, and splice it in.
  4. Unmatched Saccharomycotina taxids under that MRCA (budding yeasts absent from Shen) are reattached
     as sisters to a congeneric matched taxid when one exists, else at the grafted clade's base
     ("leave the rest near their NCBI position"). UndatedDL ignores branch lengths, so none are kept.
  5. Write data/species_tree/grafted_species.nwk with the SAME taxid leaf set as the input (verified).

RUN (Mac/any; dendropy only -- no x86 needed):
  python workflow/16_graft_species_tree.py
NEXT (on the x86/WSL2 box):
  python workflow/14_prep_generax.py --genetree results/seq_tree/pif1_aa3di.treefile --label aa3di_phylo
  # point the families species tree at data/species_tree/grafted_species.nwk, then:
  # GeneRax 2.0.4 (env gx204): generax --strategy SPR -r UndatedDL ...   (2.1.3 segfaults -- see WINDOWS_RESULTS.md)
"""
import argparse, csv, re, sys
import dendropy


def norm(name):
    """'Lachancea lanzarotensis (Yeast) (...)' -> 'lachancea_lanzarotensis' (matches workflow/17)."""
    name = re.sub(r"\(.*?\)", " ", name or "")
    parts = re.sub(r"[^A-Za-z0-9 ]", " ", name).lower().split()
    return "_".join(parts[:2]) if len(parts) >= 2 else (parts[0] if parts else "")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ncbi", default="data/species_tree/ncbi_species.covered.nwk")
    ap.add_argument("--shen", default="data/species_tree/shen2018_timetree.newick")
    ap.add_argument("--tip-map", default="data/seqs/tip_map.tsv")
    ap.add_argument("--out", default="data/species_tree/grafted_species.nwk")
    args = ap.parse_args()

    # taxid -> normalized species
    taxid2norm = {}
    for r in csv.DictReader(open(args.tip_map), delimiter="\t"):
        t = str(r.get("taxid", "")).strip()
        if t:
            taxid2norm.setdefault(t, norm(r.get("organism", "")))

    ncbi = dendropy.Tree.get(path=args.ncbi, schema="newick", preserve_underscores=True)
    shen = dendropy.Tree.get(path=args.shen, schema="newick", preserve_underscores=True)
    ncbi.is_rooted = True
    shen.is_rooted = True

    ncbi_leaf_labels = [lf.taxon.label for lf in ncbi.leaf_node_iter()]
    shen_norm = {}                                  # normalized species -> Shen tip label
    for lf in shen.leaf_node_iter():
        shen_norm.setdefault(norm(lf.taxon.label), lf.taxon.label)

    named = sum(1 for tx in ncbi_leaf_labels if taxid2norm.get(str(tx)))
    matched, claimed = {}, set()                    # taxid -> Shen tip (ONE representative taxid per Shen species)
    for tx in ncbi_leaf_labels:
        nm = taxid2norm.get(str(tx))
        if nm and nm in shen_norm:
            sl = shen_norm[nm]
            if sl not in claimed:                   # first strain represents the species on the Shen tip
                claimed.add(sl)
                matched[str(tx)] = sl
            # additional same-species strains fall through to 'unmatched' and are reattached as sisters
    sys.stderr.write(f"[16] NCBI leaves: {len(ncbi_leaf_labels)} | with a species name: {named} | "
                     f"matched to Shen: {len(matched)}\n")
    if len(matched) < 3:
        sys.exit("[16] too few matches -- check that species-tree leaves are the taxids in tip_map.tsv.")

    # clade to replace = NCBI MRCA of matched taxids
    mrca = ncbi.mrca(taxon_labels=list(matched.keys()))
    clade_leaves = [l.taxon.label for l in mrca.leaf_iter()]
    unmatched = [str(tx) for tx in clade_leaves if str(tx) not in matched]
    n_nonsacch = len(ncbi_leaf_labels) - len(clade_leaves)
    sys.stderr.write(f"[16] Saccharomycotina clade (NCBI MRCA of matched): {len(clade_leaves)} leaves = "
                     f"{len(matched)} matched + {len(unmatched)} unmatched; non-Saccharomycotina: {n_nonsacch}\n")

    # Shen subtree pruned to matched, tips relabeled Genus_species -> taxid, re-parsed into NCBI namespace
    shen_sub = shen.extract_tree_with_taxa_labels(labels=list(matched.values()))
    shen2taxid = {v: k for k, v in matched.items()}
    for lf in shen_sub.leaf_node_iter():
        lf.taxon.label = shen2taxid[lf.taxon.label]
    sub_newick = shen_sub.as_string(schema="newick", suppress_rooting=True, unquoted_underscores=True)
    sub = dendropy.Tree.get(data=sub_newick, schema="newick", taxon_namespace=ncbi.taxon_namespace,
                            preserve_underscores=True)

    # reattach unmatched as sister to a congeneric matched taxid, else at the grafted clade base
    genus_of = {tx: (taxid2norm.get(str(tx), "").split("_")[0]) for tx in clade_leaves}
    matched_by_genus = {}
    for tx in matched:
        matched_by_genus.setdefault(genus_of.get(tx, ""), tx)
    reatt_congener = 0
    for tx in unmatched:
        tax = ncbi.taxon_namespace.get_taxon(label=str(tx))
        host_taxid = matched_by_genus.get(genus_of.get(tx, ""))
        host = sub.find_node_with_taxon_label(str(host_taxid)) if host_taxid else None
        if host is not None and host.parent_node is not None:
            par = host.parent_node
            par.remove_child(host)
            newint = par.new_child()
            newint.add_child(host)
            newint.new_child(taxon=tax)
            reatt_congener += 1
        else:
            sub.seed_node.new_child(taxon=tax)

    # splice the grafted clade in place of the NCBI MRCA
    par = mrca.parent_node
    if par is None:
        grafted = sub
    else:
        par.remove_child(mrca)
        par.add_child(sub.seed_node)
        grafted = ncbi

    # verify the taxid leaf set is preserved exactly
    out_leaves = set(l.taxon.label for l in grafted.leaf_node_iter())
    in_leaves = set(str(x) for x in ncbi_leaf_labels)
    missing, extra = in_leaves - out_leaves, out_leaves - in_leaves
    sys.stderr.write(f"[16] reattached unmatched: {reatt_congener} to a congener, "
                     f"{len(unmatched) - reatt_congener} at clade base\n")
    sys.stderr.write(f"[16] leaf-set check: out {len(out_leaves)} vs in {len(in_leaves)} | "
                     f"missing {len(missing)} | extra {len(extra)}\n")
    if missing:
        sys.stderr.write(f"  WARNING missing taxids (first 10): {sorted(missing)[:10]}\n")
    if extra:
        sys.stderr.write(f"  WARNING extra taxids (first 10): {sorted(extra)[:10]}\n")

    grafted.write(path=args.out, schema="newick", suppress_rooting=True, unquoted_underscores=True)
    ok = (not missing and not extra)
    print(f"[16] wrote {args.out}  ({len(out_leaves)} taxid leaves) "
          f"{'OK -- leaf set preserved' if ok else 'CHECK WARNINGS ABOVE'}")
    print(f"     {len(matched)} Shen-placed + {len(unmatched)} reattached + {n_nonsacch} non-Saccharomycotina")


if __name__ == "__main__":
    main()
