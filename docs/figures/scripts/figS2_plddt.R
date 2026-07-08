#!/usr/bin/env Rscript
# figS2_plddt.R  --  Figure S2: structure-source consistency + helicase-core confidence.
#
# Two worries about the structural signal: (1) mixing predictors (AlphaFold DB, 828, vs ColabFold AF2,
# 129) could bias the trees; (2) low-confidence AlphaFold regions could feed noise into the 3Di.
# This panel addresses both. Full-length mean pLDDT is near-identical across the two predictors
# (predictor-consistent, ~63), and the region actually used - the helicase core - is uniformly
# HIGH-confidence (median ~88) in BOTH sources, well above the pLDDT-70 "confident" line. So the 3Di
# signal comes from well-modelled residues regardless of predictor. (Defuses review item M6.)
#
# Source: manifest.csv (full-length mean_plddt + structure_source), results/reviewer/core_plddt.tsv (core pLDDT).
# Toolchain: R 4.4.3 arm64 + ggplot2 + svglite. Vector-first; gs-outline for PDF.
# Run from repo root:  /usr/local/bin/Rscript docs/figures/scripts/figS2_plddt.R   ·  Date: 2026-07-08.

suppressMessages({ library(ggplot2); library(svglite) })

m  <- read.csv("manifest.csv", stringsAsFactors = FALSE)
cp <- read.delim("results/reviewer/core_plddt.tsv", stringsAsFactors = FALSE)
lev <- c("AFDB", "AF2_ColabFold")
xlab <- c(AFDB = "AlphaFold DB\n(n=828)", AF2_ColabFold = "ColabFold AF2\n(n=129)")
col  <- c(AFDB = "#0072B2", AF2_ColabFold = "#E69F00")

dd <- rbind(
  data.frame(source = m$structure_source, plddt = m$mean_plddt, region = "full-length model"),
  data.frame(source = cp$source, plddt = cp$core_plddt, region = "helicase core"))
dd <- dd[dd$source %in% lev & !is.na(dd$plddt), ]
dd$src    <- factor(dd$source, levels = lev)
dd$region <- factor(dd$region, levels = c("full-length model", "helicase core"))
meds <- aggregate(plddt ~ region + src, dd, median)
cat("medians:\n"); print(meds)

p <- ggplot(dd, aes(src, plddt, fill = src)) +
  geom_hline(yintercept = c(70, 90), linetype = "dotted", colour = "grey72", linewidth = 0.3) +
  geom_violin(width = 0.85, alpha = 0.45, colour = NA, scale = "width") +
  geom_boxplot(width = 0.14, outlier.size = 0.25, linewidth = 0.3, fill = "white") +
  geom_text(data = meds, aes(x = src, y = 98, label = round(plddt)), inherit.aes = FALSE,
            size = 2.5, fontface = "bold", colour = "grey20") +
  facet_wrap(~region) +
  scale_fill_manual(values = col, guide = "none") +
  scale_x_discrete(labels = xlab) +
  scale_y_continuous(breaks = seq(30, 100, 10)) +
  coord_cartesian(ylim = c(30, 100)) +
  labs(x = NULL, y = "mean pLDDT",
       title = "Structure-source consistency and helicase-core confidence",
       subtitle = "The two predictors agree, and the helicase core (the region actually used) is uniformly high-confidence",
       caption = "Bold numbers = medians. Dotted lines = pLDDT 70 (AlphaFold 'confident') and 90 ('very high'). Core median ~88 in both sources.") +
  theme_bw(base_size = 8.5) +
  theme(plot.title    = element_text(face = "bold", size = 9),
        plot.subtitle = element_text(size = 6.3, colour = "grey30"),
        plot.caption  = element_text(size = 5.9, colour = "grey40", hjust = 0),
        strip.text    = element_text(face = "bold", size = 8),
        strip.background = element_rect(fill = "grey93", colour = NA),
        panel.grid.minor = element_blank(),
        axis.text.x   = element_text(colour = "black", size = 7),
        plot.margin   = margin(3, 4, 3, 3))

ggsave("docs/figures/figS2_plddt.svg", p, width = 128, height = 86, units = "mm",
       device = svglite::svglite)
cat("wrote docs/figures/figS2_plddt.svg\n")
