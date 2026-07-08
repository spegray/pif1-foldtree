#!/usr/bin/env Rscript
# fig4_reconciliation.R  --  Figure 4: the duplication placement + recurrent losses on the species tree.
#
# The 719-taxon grafted phylogenomic species tree (Shen 2018 budding-yeast topology spliced into the NCBI
# backbone), with GeneRax 2.0.4 (UndatedDL) events mapped on. Coral nodes = duplications, area proportional
# to the per-node count D. The two interpretable events are the Saccharomycotina ancestor (node_448, D=43,
# the top duplication node in the tree) and an INDEPENDENT Agaricomycotina/mushroom ancestor (node_578,
# D=23); smaller coral dots elsewhere are reconciliation absorbing residual deep gene-tree noise as many
# small events (global D=476, SL=3618 - the interpretable signal is the two big nodes, not the totals).
# Within each duplication clade, tips are marked by present-day copy state: two-copy species retained both
# paralogs (filled), single-copy species dropped one (open) - the recurrent single-paralog loss.
#
# Source: results/reconciliation/aa3di_grafted/run/ (per_species_event_counts.txt + starting_species_tree.newick).
# Toolchain: R 4.4.3 arm64 + ggtree/treeio + ape + svglite. Vector-first; gs-outline for PDF.
# Run from repo root:  /usr/local/bin/Rscript docs/figures/scripts/fig4_reconciliation.R   ·  Date: 2026-07-08.

suppressMessages({ library(ape); library(ggtree); library(ggplot2); library(svglite) })

man     <- read.csv("manifest.csv", stringsAsFactors = FALSE)
tax2sub <- tapply(man$subphylum, man$taxid, function(x) x[1])
copies  <- table(man$taxid)                                   # proteins per species = copy number
teal <- "#009E73"; blue <- "#0072B2"; coral <- "#D55E00"

run <- "results/reconciliation/aa3di_grafted/run"
st  <- read.tree(file.path(run, "species_trees/starting_species_tree.newick"))

# --- parse per-node duplication counts ---
L  <- readLines(file.path(run, "per_species_event_counts.txt"))
pf <- function(line, key) { m <- regmatches(line, regexpr(paste0(key, "=[0-9]+"), line))
                            if (length(m)) as.integer(sub(paste0(key, "="), "", m)) else NA }
nm <- sub(" .*", "", L); Dv <- sapply(L, pf, "D"); names(Dv) <- nm

nodenum  <- function(lbl) Ntip(st) + which(st$node.label == lbl)
n448 <- nodenum("node_448"); n578 <- nodenum("node_578")     # Saccharomycotina / mushroom ancestors

all_lab <- c(st$tip.label, st$node.label)
df <- data.frame(node = seq_len(Ntip(st) + st$Nnode), label = all_lab, stringsAsFactors = FALSE)
df$D <- as.integer(Dv[df$label]); df$D[is.na(df$D)] <- 0L
df$lab <- NA_character_
df$lab[df$node == n448] <- "Saccharomycotina  D=43"
df$lab[df$node == n578] <- "Agaricomycotina  D=23"
# tip copy state, only within the two duplication clades
dup_tips <- union(extract.clade(st, n448)$tip.label, extract.clade(st, n578)$tip.label)
df$copy  <- NA_character_
ti <- df$node <= Ntip(st)
df$copy[ti] <- ifelse(df$label[ti] %in% dup_tips,
                      ifelse(as.integer(copies[df$label[ti]]) >= 2, "retained both", "lost one paralog"),
                      NA)

cat(sprintf("species tree %d tips; node_448 D=%d, node_578 D=%d; dup-clade tips %d\n",
            Ntip(st), Dv["node_448"], Dv["node_578"], length(dup_tips)))

p <- ggtree(st, layout = "fan", branch.length = "none", open.angle = 14, linewidth = 0.12,
            colour = "grey60") %<+% df +
  geom_hilight(node = n448, fill = teal, alpha = 0.13, to.bottom = TRUE) +
  geom_hilight(node = n578, fill = blue, alpha = 0.13, to.bottom = TRUE) +
  geom_tippoint(aes(subset = !is.na(copy), shape = copy, colour = copy), size = 0.7, stroke = 0.35) +
  geom_point2(aes(subset = (D > 0), size = D, alpha = D), colour = coral, stroke = 0) +
  geom_text2(aes(subset = !is.na(lab), label = lab), size = 2.4, fontface = "bold",
             colour = "grey15", hjust = 1.05) +
  scale_size_area(name = "duplications (D)", max_size = 7, breaks = c(5, 20, 43)) +
  scale_alpha_continuous(range = c(0.2, 0.85), guide = "none") +
  scale_shape_manual(name = "copy state", values = c("retained both" = 16, "lost one paralog" = 1)) +
  scale_colour_manual(name = "copy state", values = c("retained both" = "grey20", "lost one paralog" = "grey55")) +
  guides(shape = guide_legend(order = 2), colour = guide_legend(order = 2),
         size = guide_legend(order = 1)) +
  labs(title = "PIF1/RRM3 duplications and losses on the fungal species tree",
       subtitle = "Two independent duplications: the Saccharomycotina ancestor (D=43) and the mushroom ancestor (D=23)",
       caption = paste0("719-taxon grafted phylogenomic tree (GeneRax UndatedDL). Coral node area = per-node duplications; the two large nodes are the\n",
                        "interpreted events, the many faint ones are fine-scale reconciliation noise (global D=476). Open tips within a clade = single-copy reversion (loss).")) +
  theme(plot.title    = element_text(size = 9, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 6.6, colour = "grey30", hjust = 0.5),
        plot.caption  = element_text(size = 6, colour = "grey40", hjust = 0.5),
        legend.position = "right", legend.text = element_text(size = 6.5),
        legend.title = element_text(size = 7), legend.key.size = unit(3.2, "mm"),
        plot.margin = margin(2, 2, 2, 2))

ggsave("docs/figures/fig4_reconciliation.svg", p, width = 174, height = 132, units = "mm",
       device = svglite::svglite)
cat("wrote docs/figures/fig4_reconciliation.svg\n")
