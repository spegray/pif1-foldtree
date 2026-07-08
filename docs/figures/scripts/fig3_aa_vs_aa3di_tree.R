#!/usr/bin/env Rscript
# fig3_aa_vs_aa3di_tree.R  --  Figure 3A (HEADLINE): the duplication node moves when structure is added.
#
# Two fan gene trees of the SAME 957 cellular PIF1-family proteins, human-rooted:
#   (a) amino acids alone (IQ-TREE LG+I+G4, pif1.treefile)
#   (b) amino acids + 3Di structure (partitioned, pif1_aa3di.treefile)
# Coral node = MRCA(ScPif1 P07271, ScRrm3 P38766) = inferred Pif1/Rrm3 duplication. Tips are coloured by
# paralog identity defined on the AA+3Di tree: the Pif1 daughter clade (node 971, 98 tips, teal) and the
# Rrm3 daughter clade (node 1731, 99 tips, purple). Amino acids scatter the two paralog groups (their MRCA
# = base of Fungi, 950 tips); adding 3Di makes them adjacent sister clades (MRCA = Saccharomycotina, 197).
# Verified on-disk (human-rooted): AA MRCA -> 950 tips; AA+3Di MRCA -> 197 (all Saccharomycotina),
# splitting into a Pif1 clade (98) and an Rrm3 clade (99).
#
# Toolchain: R 4.4.3 arm64 + ggtree/treeio + ape + patchwork + svglite. Vector-first; gs-outline for PDF.
# Run from repo root:  /usr/local/bin/Rscript docs/figures/scripts/fig3_aa_vs_aa3di_tree.R
# Date: 2026-07-08.

suppressMessages({
  library(ape); library(ggtree); library(ggplot2); library(patchwork); library(svglite)
})

man <- read.csv("manifest.csv", stringsAsFactors = FALSE)
tm  <- read.delim("data/seqs/tip_map.tsv", stringsAsFactors = FALSE)
tip2acc <- setNames(tm$accession, tm$tip_label)
acc2sub <- setNames(man$subphylum, man$accession)
subof   <- function(tips) { s <- acc2sub[tip2acc[tips]]; s[is.na(s)] <- "?"; s }

teal <- "#009E73"; purple <- "#CC79A7"; coral <- "#D55E00"   # Okabe-Ito, CVD-safe

prep <- function(path) {
  tr <- read.tree(path)
  ht <- tr$tip.label[subof(tr$tip.label) == "Craniata"]      # human PIF1 outgroup
  root(tr, outgroup = ht, resolve.root = TRUE)
}

# --- define paralog identity ONCE, from the AA+3Di tree's two daughter clades ---
a3d <- prep("results/seq_tree/pif1_aa3di.treefile")
dup0 <- getMRCA(a3d, which(grepl("P07271", a3d$tip.label) | grepl("P38766", a3d$tip.label)))
kids <- a3d$edge[a3d$edge[, 1] == dup0, 2]
tipset <- function(k) if (k <= Ntip(a3d)) a3d$tip.label[k] else extract.clade(a3d, k)$tip.label
k_by <- sapply(kids, function(k) any(grepl("P07271", tipset(k))))   # TRUE = the Pif1-containing child
PIF1_TIPS <- tipset(kids[which(k_by)])
RRM3_TIPS <- tipset(kids[which(!k_by)])
cat(sprintf("Pif1 clade %d tips | Rrm3 clade %d tips\n", length(PIF1_TIPS), length(RRM3_TIPS)))

panel <- function(tr, title, highlight) {
  dat <- data.frame(label   = tr$tip.label,
                    paralog = ifelse(tr$tip.label %in% PIF1_TIPS, "Pif1 clade",
                              ifelse(tr$tip.label %in% RRM3_TIPS, "Rrm3 clade", NA)),
                    anchor  = ifelse(grepl("P07271", tr$tip.label), "ScPif1",
                              ifelse(grepl("P38766", tr$tip.label), "ScRrm3", NA)),
                    stringsAsFactors = FALSE)
  dup  <- getMRCA(tr, which(grepl("P07271", tr$tip.label) | grepl("P38766", tr$tip.label)))
  nd   <- length(extract.clade(tr, dup)$tip.label)
  frac <- round(mean(subof(extract.clade(tr, dup)$tip.label) == "Saccharomycotina") * 100)

  p <- ggtree(tr, layout = "fan", branch.length = "none", open.angle = 16, linewidth = 0.13,
              colour = "grey55") %<+% dat
  if (highlight) p <- p + geom_hilight(node = dup, fill = coral, alpha = 0.12, to.bottom = TRUE)
  p <- p +
    geom_tippoint(aes(subset = !is.na(paralog), colour = paralog), size = 0.6, stroke = 0) +
    geom_point2(aes(subset = (node == dup), colour = "Pif1/Rrm3 duplication"), size = 2.6) +
    geom_tiplab2(aes(subset = !is.na(anchor), label = anchor), colour = "grey15", size = 2.2,
                 offset = 1, fontface = "bold") +
    scale_colour_manual(name = NULL,
      breaks = c("Pif1 clade", "Rrm3 clade", "Pif1/Rrm3 duplication"),
      values = c("Pif1 clade" = teal, "Rrm3 clade" = purple, "Pif1/Rrm3 duplication" = coral),
      guide  = guide_legend(override.aes = list(size = c(2.4, 2.4, 2.6)))) +
    labs(title = title,
         subtitle = sprintf("Pif1/Rrm3 MRCA -> %d of 957 tips (%d%% Saccharomycotina)", nd, frac)) +
    theme(plot.title    = element_text(hjust = 0.5, size = 9,   face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size = 6.6, colour = "grey30"),
          plot.margin   = margin(1, 1, 1, 1))
  cat(sprintf("[%s] dup %d -> %d tips (%d%% Sacch)\n", title, dup, nd, frac))
  p
}

aa <- prep("results/seq_tree/pif1.treefile")

pA <- panel(aa,  "amino acids alone",           highlight = FALSE)
pB <- panel(a3d, "amino acids + 3Di structure", highlight = TRUE)

fig <- (pA | pB) +
  plot_layout(guides = "collect") +
  plot_annotation(
    caption = paste("Human-rooted cladograms of 957 PIF1-family proteins; tips coloured by paralog",
                    "(Pif1 vs Rrm3 daughter clade). Adding the 3Di structural alphabet collapses the",
                    "Pif1/Rrm3 duplication from the base of Fungi to Saccharomycotina."),
    theme = theme(plot.caption = element_text(size = 6, colour = "grey40", hjust = 0.5))) &
  theme(legend.position = "bottom", legend.text = element_text(size = 7),
        legend.key.size = unit(3, "mm"))

ggsave("docs/figures/fig3_aa_vs_aa3di_tree.svg", fig, width = 174, height = 100, units = "mm",
       device = svglite::svglite)
cat("wrote docs/figures/fig3_aa_vs_aa3di_tree.svg\n")
