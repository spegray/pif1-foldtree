#!/usr/bin/env Rscript
# fig2_copynumber.R  --  Figure 2: PIF1-family copy number by fungal subphylum.
#
# Provenance: reads manifest.csv (957 cellular PIF1-family proteins / 728 species); copy number per
#   species = number of proteins sharing a taxid. Bins species into 1 / 2 / 3+ copies and shows the
#   per-subphylum proportion. Reproduces manuscript R1 (Saccharomycotina 91/105 two-or-more-copy;
#   Pezizomycotina 430/447 single; Taphrinomycotina 11/11 single) and surfaces the convergent
#   Agaricomycetes (mushroom) multi-copy expansion that copy-counting alone cannot distinguish from
#   shared ancestry with the yeast pair.
# Toolchain: R 4.4.3 (native arm64) + ggplot2 4.0.2 + svglite; Okabe-Ito colours; vector-first.
#   Output SVG -> (rsvg-convert) PDF -> (gs -dNoOutputFonts) outlined submission PDF; QC with pdffonts.
# Run from repo root:  /usr/local/bin/Rscript docs/figures/scripts/fig2_copynumber.R
# Date: 2026-07-08.

suppressMessages({ library(ggplot2); library(svglite) })

m <- read.csv("manifest.csv", stringsAsFactors = FALSE)

# --- one row per species: subphylum + copy count ---
tx  <- unique(m$taxid)
sp  <- data.frame(taxid = tx, stringsAsFactors = FALSE)
sp$subphylum <- m$subphylum[match(sp$taxid, m$taxid)]
sp$copies    <- as.integer(table(m$taxid)[as.character(sp$taxid)])
sp$class     <- factor(ifelse(sp$copies == 1, "1 copy",
                       ifelse(sp$copies == 2, "2 copies", "3+ copies")),
                       levels = c("1 copy", "2 copies", "3+ copies"))

# --- keep well-sampled subphyla (n >= 10 species); omit human outgroup + sparse lineages ---
keep <- names(which(table(sp$subphylum) >= 10))
sp2  <- sp[sp$subphylum %in% keep & sp$subphylum != "", ]
omitted <- nrow(sp) - nrow(sp2)

# phylogenetic-ish top-to-bottom order (Ascomycota, Basidiomycota, Mucoromycota); star clade on top
ord_top <- c("Saccharomycotina", "Pezizomycotina", "Taphrinomycotina",
             "Agaricomycotina", "Ustilaginomycotina", "Pucciniomycotina",
             "Mucoromycotina", "Glomeromycotina")
ord_top <- ord_top[ord_top %in% keep]
sp2$subphylum <- factor(sp2$subphylum, levels = rev(ord_top))   # rev: geom draws first level at bottom

nlab <- as.data.frame(table(sp2$subphylum)); names(nlab) <- c("subphylum", "n")

okabe <- c("1 copy" = "#0072B2", "2 copies" = "#E69F00", "3+ copies" = "#D55E00")

p <- ggplot(sp2, aes(y = subphylum, fill = class)) +
  geom_bar(position = "fill", width = 0.72) +
  geom_text(data = nlab, aes(y = subphylum, x = 1.015, label = paste0("n=", n)),
            inherit.aes = FALSE, hjust = 0, size = 2.3, colour = "grey25") +
  scale_fill_manual(values = okabe) +
  scale_x_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.09))) +
  coord_cartesian(clip = "off") +   # let the n= labels draw into the right margin, not get clipped
  labs(x = "species (proportion)", y = NULL, fill = "PIF1-family copies",
       title = "PIF1-family copy number by fungal subphylum",
       caption = paste0(length(ord_top), " subphyla with n>=10 species shown; ",
                        omitted, " outgroup/sparsely-sampled species omitted.")) +
  theme_classic(base_size = 8.5) +
  theme(plot.title    = element_text(face = "bold", size = 9.5),
        plot.caption  = element_text(size = 6, colour = "grey40", hjust = 0),
        axis.text.y   = element_text(colour = "black"),
        axis.text.x   = element_text(colour = "black"),
        legend.position = "bottom",
        legend.key.size = unit(3.4, "mm"),
        legend.title  = element_text(size = 7.5),
        legend.margin = margin(t = -2),
        panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.3),
        plot.margin   = margin(5, 16, 4, 4))

ggsave("docs/figures/fig2_copynumber.svg", p, width = 120, height = 74, units = "mm",
       device = svglite::svglite)

cat("wrote docs/figures/fig2_copynumber.svg  (", nrow(sp2), "species,",
    length(ord_top), "subphyla;", omitted, "omitted )\n")
