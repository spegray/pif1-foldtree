#!/usr/bin/env Rscript
# analyze_lba.R -- interpret the outgroup-removal (LBA) tree from workflow/18.
# Question: with the long-branch outgroups (human + Mucoromycota + Chytridiomycota) removed and only the
# nearer Basidiomycota kept to root, does the AA-only tree collapse the ScPif1/ScRrm3 MRCA to a
# Saccharomycotina clade (as AA+3Di does, 197 tips), or leave it deep (as the full AA tree does, 950 tips)?
suppressMessages(library(ape))

tr  <- read.tree("results/lba/lba_noEDF.treefile")
tm  <- read.delim("data/seqs/tip_map.tsv", stringsAsFactors = FALSE)
man <- read.csv("manifest.csv", stringsAsFactors = FALSE)
grp     <- setNames(tm$group, tm$tip_label)
tip2acc <- setNames(tm$accession, tm$tip_label)
acc2sub <- setNames(man$subphylum, man$accession)
subof   <- function(t) { s <- acc2sub[tip2acc[t]]; s[is.na(s)] <- "?"; s }

cat("tips:", Ntip(tr), "\n")
basidio <- tr$tip.label[grp[tr$tip.label] == "basidiomycota"]
cat("basidiomycota tips:", length(basidio), "  monophyletic (unrooted):", is.monophyletic(tr, basidio), "\n")

tr2 <- tryCatch(root(tr, outgroup = basidio, resolve.root = TRUE),
                error = function(e) { cat("[root] clade-root failed; rooting on one basidio tip\n")
                                      root(tr, outgroup = basidio[1], resolve.root = TRUE) })

anchors <- which(grepl("P07271", tr2$tip.label) | grepl("P38766", tr2$tip.label))
cat("anchors:", paste(tr2$tip.label[anchors], collapse = ", "), "\n")
mrca  <- getMRCA(tr2, anchors)
clade <- extract.clade(tr2, mrca)$tip.label
cat("\n=== MRCA(ScPif1, ScRrm3) clade ===\n")
cat("clade size:", length(clade), "tips\n")
print(sort(table(subof(clade)), decreasing = TRUE))
lab <- tr2$node.label[mrca - Ntip(tr2)]
cat("MRCA node support (SH-aLRT/UFBoot):", lab, "\n")

sac <- tr2$tip.label[subof(tr2$tip.label) == "Saccharomycotina"]
cat("\nSaccharomycotina tips in tree:", length(sac), "  monophyletic:", is.monophyletic(tr2, sac), "\n")
frac_sac <- mean(subof(clade) == "Saccharomycotina")
cat(sprintf("fraction of MRCA clade that is Saccharomycotina: %.1f%% (%d/%d)\n",
            100 * frac_sac, sum(subof(clade) == "Saccharomycotina"), length(clade)))
cat("\nINTERPRETATION: clade ~ 197 Saccharomycotina  => deep AA placement was long-branch attraction.\n")
cat("                clade still deep/large + mixed  => deep signal is not an outgroup-length artifact.\n")
