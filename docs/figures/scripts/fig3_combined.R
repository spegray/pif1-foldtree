#!/usr/bin/env Rscript
# fig3_combined.R  --  Figure 3 (A + B), the headline resolution + the long-branch cause.
#   A: two human-rooted fan cladograms (amino acids alone vs amino acids + 3Di) of all 957 PIF1-family
#      proteins; the Pif1/Rrm3 duplication node (coral) moves from the base of Fungi (950 tips) to
#      Saccharomycotina (197 tips). Tips coloured by paralog: Pif1 clade (teal) / Rrm3 clade (purple).
#   B: the 197-tip Saccharomycotina duplication clade as an amino-acid phylogram (aa_bl.treefile) - the
#      Rrm3 paralogs sit on longer branches (1.22x faster; root-to-tip 1.07 vs 0.88 subs/site), the
#      long-branch-attraction signal that misleads amino-acid trees and that structure corrects.
# Standalone panel scripts: fig3_aa_vs_aa3di_tree.R (A), fig3b_branchlength.R (B).
# Run from repo root:  /usr/local/bin/Rscript docs/figures/scripts/fig3_combined.R   ·  Date: 2026-07-08.

suppressMessages({ library(ape); library(ggtree); library(ggplot2); library(patchwork); library(svglite) })

man <- read.csv("manifest.csv", stringsAsFactors = FALSE)
tm  <- read.delim("data/seqs/tip_map.tsv", stringsAsFactors = FALSE)
tip2acc <- setNames(tm$accession, tm$tip_label); acc2sub <- setNames(man$subphylum, man$accession)
subof <- function(t) { s <- acc2sub[tip2acc[t]]; s[is.na(s)] <- "?"; s }
teal <- "#009E73"; purple <- "#CC79A7"; coral <- "#D55E00"
prep <- function(p) { tr <- read.tree(p); root(tr, outgroup = tr$tip.label[subof(tr$tip.label) == "Craniata"], resolve.root = TRUE) }

a3d  <- prep("results/seq_tree/pif1_aa3di.treefile")
dup0 <- getMRCA(a3d, which(grepl("P07271", a3d$tip.label) | grepl("P38766", a3d$tip.label)))
kids <- a3d$edge[a3d$edge[, 1] == dup0, 2]
tipset <- function(k) if (k <= Ntip(a3d)) a3d$tip.label[k] else extract.clade(a3d, k)$tip.label
kb <- sapply(kids, function(k) any(grepl("P07271", tipset(k))))
PIF1_TIPS <- tipset(kids[which(kb)]); RRM3_TIPS <- tipset(kids[which(!kb)])

fan <- function(tr, title, highlight) {
  dat <- data.frame(label = tr$tip.label,
                    paralog = ifelse(tr$tip.label %in% PIF1_TIPS, "Pif1 clade",
                              ifelse(tr$tip.label %in% RRM3_TIPS, "Rrm3 clade", NA)),
                    anchor  = ifelse(grepl("P07271", tr$tip.label), "ScPif1",
                              ifelse(grepl("P38766", tr$tip.label), "ScRrm3", NA)))
  dup <- getMRCA(tr, which(grepl("P07271", tr$tip.label) | grepl("P38766", tr$tip.label)))
  nd  <- length(extract.clade(tr, dup)$tip.label)
  frac <- round(mean(subof(extract.clade(tr, dup)$tip.label) == "Saccharomycotina") * 100)
  p <- ggtree(tr, layout = "fan", branch.length = "none", open.angle = 16, linewidth = 0.13, colour = "grey55") %<+% dat
  if (highlight) p <- p + geom_hilight(node = dup, fill = coral, alpha = 0.12, to.bottom = TRUE)
  p + geom_tippoint(aes(subset = !is.na(paralog), colour = paralog), size = 0.55, stroke = 0) +
    geom_point2(aes(subset = (node == dup), colour = "Pif1/Rrm3 duplication"), size = 2.4) +
    geom_tiplab2(aes(subset = !is.na(anchor), label = anchor), colour = "grey15", size = 2, offset = 1, fontface = "bold") +
    scale_colour_manual(name = NULL, breaks = c("Pif1 clade", "Rrm3 clade", "Pif1/Rrm3 duplication"),
      values = c("Pif1 clade" = teal, "Rrm3 clade" = purple, "Pif1/Rrm3 duplication" = coral),
      guide = guide_legend(override.aes = list(size = c(2.4, 2.4, 2.4)))) +
    labs(title = title, subtitle = sprintf("MRCA -> %d of 957 tips (%d%% Saccharomycotina)", nd, frac)) +
    theme(plot.title = element_text(hjust = 0.5, size = 8.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size = 6.2, colour = "grey30"), plot.margin = margin(1, 1, 1, 1))
}

bl  <- prep("results/seq_tree/aa_bl.treefile")
dup <- getMRCA(bl, which(grepl("P07271", bl$tip.label) | grepl("P38766", bl$tip.label)))
sub <- extract.clade(bl, dup)
rt  <- setNames(node.depth.edgelength(sub)[seq_len(Ntip(sub))], sub$tip.label)
r_p <- rt[names(rt) %in% PIF1_TIPS]; r_r <- rt[names(rt) %in% RRM3_TIPS]; ratio <- mean(r_r) / mean(r_p)
subg <- groupOTU(sub, list("Pif1 clade" = intersect(sub$tip.label, PIF1_TIPS),
                           "Rrm3 clade" = intersect(sub$tip.label, RRM3_TIPS)), group_name = "paralog")
anc <- data.frame(label = sub$tip.label, anchor = ifelse(grepl("P07271", sub$tip.label), "ScPif1",
                  ifelse(grepl("P38766", sub$tip.label), "ScRrm3", NA)))
pB <- ggtree(subg, aes(colour = paralog), linewidth = 0.28) %<+% anc +
  scale_colour_manual(values = c("0" = "grey70", "Pif1 clade" = teal, "Rrm3 clade" = purple), guide = "none") +
  geom_tippoint(aes(subset = !is.na(anchor)), colour = coral, size = 1.5) +
  geom_tiplab(aes(subset = !is.na(anchor), label = anchor), colour = "grey15", size = 2.2, offset = 0.03, fontface = "bold") +
  geom_treescale(x = 0, y = -3, width = 0.2, fontsize = 2.1, linesize = 0.3, offset = 2) +
  labs(title = "Branch lengths within the Saccharomycotina duplication clade",
       subtitle = sprintf("Rrm3 evolves %.2fx faster than Pif1 (root-to-tip %.2f vs %.2f subs/site) - the long-branch signal", ratio, mean(r_r), mean(r_p))) +
  theme_tree() + theme(plot.title = element_text(size = 8.5, face = "bold"),
                       plot.subtitle = element_text(size = 6.2, colour = "grey30"), plot.margin = margin(2, 4, 2, 2))

pAA  <- fan(prep("results/seq_tree/pif1.treefile"), "amino acids alone", FALSE) + labs(tag = "A") +
        theme(plot.tag = element_text(face = "bold", size = 12))
pA3D <- fan(a3d, "amino acids + 3Di structure", TRUE)
pB   <- pB + labs(tag = "B") + theme(plot.tag = element_text(face = "bold", size = 12))

fig <- ((pAA | pA3D) / pB) +
  plot_layout(heights = c(1, 1.35), guides = "collect") +
  plot_annotation(
    caption = paste0("Human-rooted trees of 957 PIF1-family proteins; tips coloured by paralog (Pif1 vs Rrm3 clade).\n",
                     "3Di structure moves the duplication from Fungi to Saccharomycotina (A); Rrm3's long branches (B) are the LBA signal."),
    theme = theme(plot.caption = element_text(size = 6, colour = "grey40", hjust = 0.5))) &
  theme(legend.position = "bottom", legend.text = element_text(size = 7), legend.key.size = unit(3, "mm"))

ggsave("docs/figures/fig3_combined.svg", fig, width = 180, height = 208, units = "mm", device = svglite::svglite)
cat(sprintf("wrote docs/figures/fig3_combined.svg  (A: fans; B: %d-tip clade, Rrm3/Pif1 = %.2fx)\n", Ntip(sub), ratio))
