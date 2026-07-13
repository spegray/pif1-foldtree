#!/usr/bin/env Rscript
# fig5_timetree.R  --  Figure 5: the PIF1/RRM3 duplication dated on the Shen 2018 chronogram.
#
# The budding-yeast time tree of Shen et al. 2018 (Cell; 332 genomes, fossil-calibrated, branch units
# x100 Myr) on a geological time axis. The duplication clade's species (MRCA of ScPif1/ScRrm3 in the
# AA+3Di tree) are matched onto the Shen tips by genus+species; their MRCA (subtends 317/332 tips) dates
# the event to [crown, stem] = [326, 383] Mya. Coral band = that bracket; the whole-genome duplication
# (~100 Mya) is marked for scale. So the Pif1/Rrm3 split is Devonian-Carboniferous, on the Saccharomycotina
# stem, ~3-4x older than the WGD. Reproduces workflow/17 (crown 326, stem 383, root/BYCA 404 Mya).
#
# Gotchas (heed): edge lengths x100 BEFORE revts(); revts() (a ggtree fn) BEFORE coord_geo(neg=TRUE);
# the Shen tree has no node CIs -> the event is a shaded band, not an HPD bar.
# Toolchain: R 4.4.3 arm64 + ggtree + deeptime + ape + svglite. Vector-first; gs-outline for PDF.
# Run from repo root:  /usr/local/bin/Rscript docs/figures/scripts/fig5_timetree.R   ·  Date: 2026-07-08.

suppressMessages({ library(ape); library(ggtree); library(deeptime); library(ggplot2); library(svglite) })

norm <- function(x) { x <- gsub("\\(.*?\\)", " ", x)
  p <- tolower(unlist(strsplit(gsub("[^A-Za-z0-9 ]", " ", x), " +"))); p <- p[p != ""]
  if (length(p) >= 2) paste(p[1], p[2], sep = "_") else if (length(p)) p[1] else "" }

man <- read.csv("manifest.csv", stringsAsFactors = FALSE)
tm  <- read.delim("data/seqs/tip_map.tsv", stringsAsFactors = FALSE)
tip2acc <- setNames(tm$accession, tm$tip_label); acc2sub <- setNames(man$subphylum, man$accession)
tip2org <- setNames(tm$organism, tm$tip_label)
subof   <- function(t) { s <- acc2sub[tip2acc[t]]; s[is.na(s)] <- "?"; s }
coral <- "#D55E00"; teal <- "#009E73"

# duplication-clade species (AA+3Di tree, human-rooted) -> normalized names
a3d <- read.tree("results/seq_tree/pif1_aa3di.treefile")
a3d <- root(a3d, outgroup = a3d$tip.label[subof(a3d$tip.label) == "Craniata"], resolve.root = TRUE)
dupA <- getMRCA(a3d, which(grepl("P07271", a3d$tip.label) | grepl("P38766", a3d$tip.label)))
dup_orgs <- unique(vapply(tip2org[extract.clade(a3d, dupA)$tip.label], norm, character(1)))

# Shen 2018 chronogram; MRCA of matched species; crown/stem ages
tt <- read.tree("data/species_tree/shen2018_timetree.newick")
tt_norm <- vapply(tt$tip.label, norm, character(1))
matched  <- tt$tip.label[tt_norm %in% dup_orgs]
DUP <- getMRCA(tt, matched)
dep <- node.depth.edgelength(tt); rootdepth <- max(dep[seq_len(Ntip(tt))])
crown <- (rootdepth - dep[DUP]) * 100
stem  <- (rootdepth - dep[tt$edge[tt$edge[, 2] == DUP, 1]]) * 100
cat(sprintf("Shen %d tips; matched %d; MRCA node %d subtends %d; crown %.0f stem %.0f root %.0f Mya\n",
            Ntip(tt), length(matched), DUP, length(extract.clade(tt, DUP)$tip.label), crown, stem, rootdepth * 100))

tt$edge.length <- tt$edge.length * 100                 # now in Myr, BEFORE revts()
NT <- Ntip(tt)

p <- revts(ggtree(tt, linewidth = 0.12, colour = "grey55")) +
  annotate("rect", xmin = -stem, xmax = -crown, ymin = 0, ymax = NT + 1, fill = coral, alpha = 0.16) +
  geom_vline(xintercept = -100, linetype = "dashed", colour = "grey45", linewidth = 0.3) +
  geom_point2(aes(subset = (node == DUP)), colour = coral, size = 2.6) +
  annotate("text", x = -(stem + crown) / 2, y = NT * 0.60,
           label = "PIF1/RRM3\nduplication\n~326-383 Mya", colour = coral, fontface = "bold",
           size = 2.5, lineheight = 0.95) +
  annotate("text", x = -100, y = NT * 0.5, label = "whole-genome duplication (~100 Mya)",
           colour = "grey35", size = 2.1, angle = 90, hjust = 0.5, vjust = -0.5) +
  scale_x_continuous(breaks = seq(-400, 0, 50), labels = abs(seq(-400, 0, 50)),
                     name = "millions of years ago") +
  labs(title = "The PIF1/RRM3 duplication dated on the budding-yeast time tree",
       subtitle = "Shen et al. 2018 chronogram (332 genomes): the split is Devonian-Carboniferous, on the Saccharomycotina stem, ~3-4x older than the WGD") +
  coord_geo(dat = "periods", neg = TRUE, abbrv = TRUE, size = 2.1, height = unit(1.1, "line"),
            xlim = c(-412, 6), ylim = c(-6, NT + 4), expand = FALSE) +
  theme_tree2() +
  theme(plot.title    = element_text(size = 9, face = "bold"),
        plot.subtitle = element_text(size = 6.3, colour = "grey30"),
        axis.title.x  = element_text(size = 7.5), axis.text.x = element_text(size = 6.5),
        plot.margin   = margin(3, 4, 2, 3))

ggsave("docs/figures/fig5_timetree.svg", p, width = 174, height = 120, units = "mm", device = svglite::svglite)
cat("wrote docs/figures/fig5_timetree.svg\n")
