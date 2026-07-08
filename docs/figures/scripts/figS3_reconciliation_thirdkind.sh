#!/usr/bin/env bash
# figS3_reconciliation_thirdkind.sh  --  Figure S3: literal recPhyloXML reconciliation (thirdkind).
#
# Cross-check of Figure 4: renders the RAW GeneRax 2.0.4 (UndatedDL) reconciliation - the species tree
# with the 941-gene family embedded, duplications and losses drawn in place - straight from the
# recPhyloXML, using thirdkind. Non-focal clades are collapsed to triangles so the two duplication
# clades (Saccharomycotina node_448, mushroom node_578) stay legible. This is intentionally the LITERAL
# output; a 197-gene reconciliation is dense by nature (thirdkind is the QC/cross-check, ggtree Fig 4 is
# the communicative figure). Confirms the Saccharomycotina placement against the actual recPhyloXML.
#
# Requires: thirdkind (cargo install thirdkind), R+ape, rsvg-convert, ghostscript, pdffonts.
# Run from repo root:  bash docs/figures/scripts/figS3_reconciliation_thirdkind.sh   ·  Date: 2026-07-08.
set -euo pipefail

RUN=results/reconciliation/aa3di_grafted/run
XML="$RUN/reconciliations/pif1_reconciliated.xml"
ST="$RUN/species_trees/starting_species_tree.newick"
OUT=docs/figures/figS3_reconciliation_thirdkind

# collapse set = internal nodes OFF the path root -> {Saccharomycotina node_448, mushroom node_578}
COLL=$(/usr/local/bin/Rscript - "$ST" <<'RS'
suppressMessages(library(ape))
st <- read.tree(commandArgs(TRUE)[1]); root <- Ntip(st) + 1
n448 <- Ntip(st) + which(st$node.label == "node_448")
n578 <- Ntip(st) + which(st$node.label == "node_578")
spine <- unique(c(nodepath(st, root, n448), nodepath(st, root, n578)))
coll <- character(0)
for (nd in spine) for (cc in st$edge[st$edge[, 1] == nd, 2])
  if (!(cc %in% spine) && cc > Ntip(st)) coll <- c(coll, st$node.label[cc - Ntip(st)])
cat(paste(coll, collapse = ","))
RS
)
echo "collapsing $(echo "$COLL" | tr ',' '\n' | wc -l | tr -d ' ') non-focal clades: $COLL"

~/.cargo/bin/thirdkind -f "$XML" -Y "$COLL" -o "$OUT.svg"

# QC: editable thirdkind SVG -> outlined submission PDF (0 font objects)
rsvg-convert -f pdf -o /tmp/figS4.pdf "$OUT.svg"
gs -q -o "${OUT}_outlined.pdf" -dNoOutputFonts -sDEVICE=pdfwrite /tmp/figS4.pdf
rm -f /tmp/figS4.pdf
echo "wrote $OUT.svg + ${OUT}_outlined.pdf  (font objects: $(pdffonts "${OUT}_outlined.pdf" | tail -n +3 | wc -l | tr -d ' '))"
