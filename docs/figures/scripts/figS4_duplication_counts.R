#!/usr/bin/env Rscript
# figS4_duplication_counts.R  --  Figure S4 (companion to the thirdkind render): the reconciliation's
# per-node duplication counts, ranked. Cross-checks Figure 4 quantitatively: the two interpreted events -
# the Saccharomycotina ancestor (node_448, D=43) and the independent mushroom ancestor (node_578, D=23) -
# tower over every other node (all D<=16), and the long tail of small-D nodes is GeneRax absorbing
# residual deep gene-tree noise (global D=476; 140 of 166 D>0 nodes have D<5). This is why Fig 4 draws
# the two events solid and fades the rest.
#
# Source: results/reconciliation/aa3di_grafted/run/ (per_species_event_counts.txt + starting_species_tree.newick).
# Toolchain: R 4.4.3 arm64 + ape + ggplot2 + svglite. Vector-first; gs-outline for PDF.
# Run from repo root:  /usr/local/bin/Rscript docs/figures/scripts/figS4_duplication_counts.R   ·  Date: 2026-07-08.

suppressMessages({ library(ape); library(ggplot2); library(svglite) })

man     <- read.csv("manifest.csv", stringsAsFactors = FALSE)
tax2sub <- tapply(man$subphylum, man$taxid, function(x) x[1])
st <- read.tree("results/reconciliation/aa3di_grafted/run/species_trees/starting_species_tree.newick")
coral <- "#D55E00"

L  <- readLines("results/reconciliation/aa3di_grafted/run/per_species_event_counts.txt")
pf <- function(line, key) { m <- regmatches(line, regexpr(paste0(key, "=[0-9]+"), line))
                            if (length(m)) as.integer(sub(paste0(key, "="), "", m)) else 0L }
nm <- sub(" .*", "", L); D <- sapply(L, pf, "D"); names(D) <- nm
D  <- D[D > 0]
n_total <- length(D); n_lt5 <- sum(D < 5)

topn <- 15
top  <- names(sort(D, decreasing = TRUE))[seq_len(topn)]

clade_sub <- function(lbl) {                       # majority subphylum subtended by a node
  idx <- which(st$node.label == lbl)
  if (length(idx) == 1) {
    d <- extract.clade(st, Ntip(st) + idx)$tip.label
    s <- tax2sub[as.character(d)]; s <- s[!is.na(s) & s != ""]
    if (length(s)) names(sort(table(s), decreasing = TRUE))[1] else "unresolved"
  } else { v <- tax2sub[lbl]; if (is.na(v) || v == "") "single lineage" else as.character(v) }
}

df <- data.frame(node = top, D = as.integer(D[top]),
                 clade = vapply(top, clade_sub, character(1)), stringsAsFactors = FALSE)
df$event  <- df$node %in% c("node_448", "node_578")
df$clabel <- ifelse(df$node == "node_448", "Saccharomycotina ancestor",
              ifelse(df$node == "node_578", "mushroom ancestor", df$clade))
df$node   <- factor(df$node, levels = rev(top))    # highest D at top

cat(sprintf("top %d nodes; node_448 D=%d, node_578 D=%d; %d nodes total, %d with D<5\n",
            topn, D["node_448"], D["node_578"], n_total, n_lt5))

p <- ggplot(df, aes(x = D, y = node, fill = event)) +
  geom_col(width = 0.68) +
  geom_text(aes(label = clabel, colour = event, fontface = ifelse(event, "bold", "plain")),
            hjust = -0.06, size = 2.4) +
  scale_fill_manual(values = c(`FALSE` = "grey78", `TRUE` = coral), guide = "none") +
  scale_colour_manual(values = c(`FALSE` = "grey45", `TRUE` = coral), guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.42))) +
  labs(x = "duplications inferred at node (D)", y = NULL,
       title = "Duplication counts per species-tree node",
       subtitle = "Two events dominate the reconciliation; every other node is fine-scale noise",
       caption = sprintf("Top %d of %d nodes with D>0 (global D=476); %d have D<5. Grey bar labels = subphylum the node subtends.",
                         topn, n_total, n_lt5)) +
  theme_classic(base_size = 8.5) +
  theme(axis.text.y  = element_blank(), axis.ticks.y = element_blank(),
        plot.title    = element_text(face = "bold", size = 9.5),
        plot.subtitle = element_text(size = 6.8, colour = "grey30"),
        plot.caption  = element_text(size = 5.8, colour = "grey40", hjust = 0),
        panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.3),
        plot.margin = margin(3, 3, 3, 3))

ggsave("docs/figures/figS4_duplication_counts.svg", p, width = 128, height = 82, units = "mm",
       device = svglite::svglite)
cat("wrote docs/figures/figS4_duplication_counts.svg\n")
