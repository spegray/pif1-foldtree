#!/usr/bin/env Rscript
# figS1_family_filter.R  --  Figure S1: the family-filter decision (InterPro IPR048293 vs Pfam PF05970).
#
# The bare Pfam PF05970 ("PIF1-like helicase") domain is ALSO carried by Helitron rolling-circle
# transposon helicases, inflating fungal hit counts several-fold. Switching the family definition to
# InterPro IPR048293 ("PIF1_RRM3_pfh1" = the CELLULAR PIF1/RRM3/Pfh1 group) removes that contamination
# while keeping every anchor. This panel shows the per-phylum protein counts before (raw PF05970 gather)
# and after (IPR048293 cellular). Counts are DERIVED from the actual files, not hardcoded:
#   before = data/seqs/_superseded_PF05970/candidates*.tsv (minus header)
#   after  = final manifest.csv, proteins per phylum
# Result: Ascomycota 2255->681, Basidiomycota 2393->235, Mucoromycota 742->36, Chytridiomycota 33->4.
#
# Toolchain: R 4.4.3 arm64 + ggplot2 + svglite. Vector-first; gs-outline for PDF.
# Run from repo root:  /usr/local/bin/Rscript docs/figures/scripts/figS1_family_filter.R   ·  Date: 2026-07-08.

suppressMessages({ library(ggplot2); library(svglite) })

grey <- "grey62"; teal <- "#009E73"
sup  <- "data/seqs/_superseded_PF05970"
nlines_minus_header <- function(f) length(readLines(f)) - 1L    # each candidates file has 1 header

before <- c(
  Ascomycota      = nlines_minus_header(file.path(sup, "candidates.tsv")),            # default gather = Ascomycota
  Basidiomycota   = nlines_minus_header(file.path(sup, "candidates_basidiomycota.tsv")),
  Mucoromycota    = nlines_minus_header(file.path(sup, "candidates_mucoromycota.tsv")),
  Chytridiomycota = nlines_minus_header(file.path(sup, "candidates_chytridiomycota.tsv")))

m <- read.csv("manifest.csv", stringsAsFactors = FALSE)
asco <- c("Saccharomycotina", "Pezizomycotina", "Taphrinomycotina")
bas  <- c("Agaricomycotina", "Ustilaginomycotina", "Pucciniomycotina", "Wallemiomycotina")
muc  <- c("Mucoromycotina", "Glomeromycotina", "Mortierellomycotina")
phyl <- function(s) ifelse(s %in% asco, "Ascomycota", ifelse(s %in% bas, "Basidiomycota",
                    ifelse(s %in% muc, "Mucoromycota", ifelse(s == "Craniata", "Homo", "Chytridiomycota"))))
after_tab <- table(phyl(m$subphylum))
after <- setNames(as.integer(after_tab[names(before)]), names(before))

ph <- names(before)
df <- data.frame(
  phylum = factor(rep(ph, 2), levels = ph[order(before, decreasing = TRUE)]),
  filter = factor(rep(c("Pfam PF05970 (raw)", "InterPro IPR048293 (cellular)"), each = length(ph)),
                  levels = c("Pfam PF05970 (raw)", "InterPro IPR048293 (cellular)")),
  n = c(before, after))
fold <- data.frame(phylum = factor(ph, levels = levels(df$phylum)),
                   lab = sprintf("%.0f%% removed", 100 * (before - after) / before),
                   y = before)

cat("before:", paste(sprintf("%s=%d", ph, before), collapse = " "), "\n")
cat("after :", paste(sprintf("%s=%d", ph, after), collapse = " "), "\n")

p <- ggplot(df, aes(x = phylum, y = n, fill = filter)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.66) +
  geom_text(aes(label = format(n, big.mark = ",")), position = position_dodge(width = 0.72),
            vjust = -0.4, size = 2.3, colour = "grey25") +
  geom_text(data = fold, aes(x = phylum, y = y, label = lab), inherit.aes = FALSE,
            vjust = -1.9, size = 2.2, fontface = "italic", colour = "grey45") +
  scale_fill_manual(name = NULL, values = c("Pfam PF05970 (raw)" = grey,
                                            "InterPro IPR048293 (cellular)" = teal)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.16))) +
  labs(x = NULL, y = "PIF1-family proteins",
       title = "Family filter: InterPro IPR048293 removes Helitron-transposon contamination",
       caption = paste0("Bare Pfam PF05970 is also carried by Helitron rolling-circle transposon helicases.\n",
                        "InterPro IPR048293 isolates the cellular PIF1/RRM3/Pfh1 group and keeps every anchor.")) +
  theme_classic(base_size = 8.5) +
  theme(plot.title    = element_text(face = "bold", size = 9),
        plot.caption  = element_text(size = 6, colour = "grey40", hjust = 0),
        legend.position = "top", legend.text = element_text(size = 7.5),
        axis.text.x   = element_text(colour = "black"),
        plot.margin   = margin(3, 4, 3, 3))

ggsave("docs/figures/figS1_family_filter.svg", p, width = 128, height = 84, units = "mm",
       device = svglite::svglite)
cat("wrote docs/figures/figS1_family_filter.svg\n")
