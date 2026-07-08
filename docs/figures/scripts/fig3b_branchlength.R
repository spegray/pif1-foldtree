#!/usr/bin/env Rscript
# fig3b_branchlength.R  --  Figure 3B: branch lengths within the Saccharomycotina duplication clade.
#
# The 197-tip Saccharomycotina Pif1/Rrm3 clade, extracted from aa_bl.treefile (the AA+3Di topology with
# amino-acid branch lengths re-estimated — the tree used for the B3 rate analysis). Rectangular phylogram,
# branches coloured by paralog: Pif1 clade (teal) vs Rrm3 clade (purple). This is the long-branch signal
# that misleads amino-acid trees: the Rrm3 paralogs sit on visibly longer branches (faster-evolving, more
# derived) than Pif1 — which is why an amino-acid-only tree gets pulled toward a spurious deep placement
# (long-branch attraction), and why adding structure (Fig 3A) is needed to recover the true node.
#
# Toolchain: R 4.4.3 arm64 + ggtree/treeio + ape + svglite. Vector-first; gs-outline for PDF.
# Run from repo root:  /usr/local/bin/Rscript docs/figures/scripts/fig3b_branchlength.R
# Date: 2026-07-08.

suppressMessages({ library(ape); library(ggtree); library(ggplot2); library(svglite) })

man <- read.csv("manifest.csv", stringsAsFactors = FALSE)
tm  <- read.delim("data/seqs/tip_map.tsv", stringsAsFactors = FALSE)
tip2acc <- setNames(tm$accession, tm$tip_label)
acc2sub <- setNames(man$subphylum, man$accession)
subof   <- function(tips) { s <- acc2sub[tip2acc[tips]]; s[is.na(s)] <- "?"; s }
teal <- "#009E73"; purple <- "#CC79A7"; coral <- "#D55E00"

prep <- function(path) {
  tr <- read.tree(path)
  ht <- tr$tip.label[subof(tr$tip.label) == "Craniata"]
  root(tr, outgroup = ht, resolve.root = TRUE)
}

# paralog identity from the AA+3Di topology (same as Fig 3A)
a3d  <- prep("results/seq_tree/pif1_aa3di.treefile")
dup0 <- getMRCA(a3d, which(grepl("P07271", a3d$tip.label) | grepl("P38766", a3d$tip.label)))
kids <- a3d$edge[a3d$edge[, 1] == dup0, 2]
tipset <- function(k) if (k <= Ntip(a3d)) a3d$tip.label[k] else extract.clade(a3d, k)$tip.label
k_by <- sapply(kids, function(k) any(grepl("P07271", tipset(k))))
PIF1_TIPS <- tipset(kids[which(k_by)]); RRM3_TIPS <- tipset(kids[which(!k_by)])

# extract the Saccharomycotina duplication clade WITH amino-acid branch lengths
bl  <- prep("results/seq_tree/aa_bl.treefile")
dup <- getMRCA(bl, which(grepl("P07271", bl$tip.label) | grepl("P38766", bl$tip.label)))
sub <- extract.clade(bl, dup)
cat(sprintf("clade: %d tips (AA branch lengths)\n", Ntip(sub)))

# rate: mean root(dup)-to-tip distance per paralog
rt   <- setNames(node.depth.edgelength(sub)[seq_len(Ntip(sub))], sub$tip.label)
r_p  <- rt[names(rt) %in% PIF1_TIPS]; r_r <- rt[names(rt) %in% RRM3_TIPS]
ratio <- mean(r_r) / mean(r_p)
cat(sprintf("root-to-tip: Pif1 %.3f | Rrm3 %.3f | ratio %.2fx\n", mean(r_p), mean(r_r), ratio))

grp  <- list("Pif1 clade" = intersect(sub$tip.label, PIF1_TIPS),
             "Rrm3 clade" = intersect(sub$tip.label, RRM3_TIPS))
subg <- groupOTU(sub, grp, group_name = "paralog")
anc  <- data.frame(label = sub$tip.label,
                   anchor = ifelse(grepl("P07271", sub$tip.label), "ScPif1",
                            ifelse(grepl("P38766", sub$tip.label), "ScRrm3", NA)),
                   stringsAsFactors = FALSE)

p <- ggtree(subg, aes(colour = paralog), linewidth = 0.28) %<+% anc +
  scale_colour_manual(name = NULL,
    breaks = c("Pif1 clade", "Rrm3 clade"),
    values = c("0" = "grey70", "Pif1 clade" = teal, "Rrm3 clade" = purple),
    guide  = guide_legend(override.aes = list(linewidth = 1.4))) +
  geom_tippoint(aes(subset = !is.na(anchor)), colour = coral, size = 1.6) +
  geom_tiplab(aes(subset = !is.na(anchor), label = anchor), colour = "grey15",
              size = 2.3, offset = 0.03, fontface = "bold") +
  geom_treescale(x = 0, y = -3, width = 0.2, fontsize = 2.2, linesize = 0.3, offset = 2) +
  labs(title = "Branch lengths within the Saccharomycotina duplication clade",
       subtitle = sprintf("Rrm3 evolves %.2fx faster than Pif1 (mean root-to-tip %.2f vs %.2f subs/site) - the long-branch signal",
                          ratio, mean(r_r), mean(r_p))) +
  theme_tree() +
  theme(plot.title    = element_text(size = 9, face = "bold"),
        plot.subtitle = element_text(size = 6.6, colour = "grey30"),
        legend.position = "right", legend.text = element_text(size = 7),
        legend.key.size = unit(3.4, "mm"),
        plot.margin   = margin(2, 6, 2, 2))

ggsave("docs/figures/fig3b_branchlength.svg", p, width = 120, height = 132, units = "mm",
       device = svglite::svglite)
cat("wrote docs/figures/fig3b_branchlength.svg\n")
