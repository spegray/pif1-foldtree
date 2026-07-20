#!/usr/bin/env Rscript
# extract_ancestor.R -- pull the pre-duplication Saccharomycotina PIF1 core from the ASR (workflow/19).
# Finds the duplication node (MRCA of ScPif1/ScRrm3), extracts its marginal ancestral states + per-site
# posterior confidence, gap-masks columns absent across the clade, maps the Arg-wedge (ScPif1 R324) column,
# and writes a fold-ready sequence + a confidence table + a summary.
suppressMessages(library(ape))

readFasta <- function(p) {
  ls <- readLines(p); out <- list(); nm <- NULL; buf <- character(0)
  for (l in ls) {
    if (startsWith(l, ">")) { if (!is.null(nm)) out[[nm]] <- paste(buf, collapse=""); nm <- sub("^>(\\S+).*","\\1",l); buf <- character(0) }
    else buf <- c(buf, l)
  }
  if (!is.null(nm)) out[[nm]] <- paste(buf, collapse=""); out
}

tm  <- read.delim("data/seqs/tip_map.tsv", stringsAsFactors=FALSE)
man <- read.csv("manifest.csv", stringsAsFactors=FALSE)
tip2acc <- setNames(tm$accession, tm$tip_label)
acc2sub <- setNames(man$subphylum, man$accession)
subof <- function(t){ s<-acc2sub[tip2acc[t]]; s[is.na(s)]<-"?"; s }

tr <- read.tree("results/asr/asr.treefile")
human <- tr$tip.label[subof(tr$tip.label) == "Craniata"]
tr <- root(tr, outgroup = human, resolve.root = TRUE)
anchors <- tr$tip.label[grepl("P07271", tr$tip.label) | grepl("P38766", tr$tip.label)]
dup <- getMRCA(tr, anchors)
dup_label <- tr$node.label[dup - Ntip(tr)]
clade_tips <- extract.clade(tr, dup)$tip.label
cat(sprintf("duplication node = %s ; clade = %d tips (%d Saccharomycotina)\n",
            dup_label, length(clade_tips), sum(subof(clade_tips)=="Saccharomycotina")))

# ancestral states for the dup node (awk-prefilter the 34 MB .state; cols: Node Site State p1..p20)
d <- read.table(pipe(sprintf("awk -F'\t' '$1==\"%s\"' results/asr/asr.state", dup_label)),
                header=FALSE, stringsAsFactors=FALSE)
d <- d[order(d[[2]]), ]
state   <- d[[3]]
maxpost <- apply(as.matrix(d[, 4:23]), 1, max)
ncol_aln <- length(state)

# gap fraction per column across the duplication clade (from the trimmed alignment)
trim <- readFasta("results/seq_tree/aln.trim.fasta")
clade_mat <- do.call(rbind, strsplit(unlist(trim[clade_tips]), ""))
gapfrac <- colMeans(clade_mat == "-")

# --- wedge (ScPif1 R324) column via untrimmed<->trimmed subsequence (trimAl dropped 1 ScPif1 residue) ---
sck <- names(trim)[grepl("P07271", names(trim))][1]
U <- gsub("-", "", readFasta("results/seq_tree/aln.fasta")[[sck]])   # 207 res = 236..442
Tg <- trim[[sck]]; Tng <- gsub("-", "", Tg)
uc <- strsplit(U,"")[[1]]; tc <- strsplit(Tng,"")[[1]]
di <- length(uc); for (i in seq_along(tc)) if (uc[i] != tc[i]) { di <- i; break }   # single deletion
u_idx <- 324 - 236 + 1                                                              # = 89
stopifnot(uc[u_idx] == "R")
t_idx <- if (u_idx < di) u_idx else if (u_idx > di) u_idx - 1L else NA
Tc <- strsplit(Tg,"")[[1]]; wedge_col <- which(Tc != "-")[t_idx]
wedge_state <- state[wedge_col]; wedge_post <- maxpost[wedge_col]

# fold-ready sequence: keep columns present across the clade (<50% gaps), drop rare-insertion columns
keep <- gapfrac < 0.5
anc_full <- paste(state, collapse="")
anc_fold <- paste(state[keep], collapse="")
# wedge position within the folded sequence
wedge_fold_pos <- sum(keep[1:wedge_col])

dir.create("results/asr", showWarnings=FALSE)
writeLines(c(sprintf(">ancestor_dupnode_Saccharomycotina_PIF1_core | %s | %d aa | pre-duplication (MRCA ScPif1/ScRrm3)", dup_label, nchar(anc_fold)),
             anc_fold), "results/asr/ancestor_dupnode.fasta")
write.table(data.frame(aln_col=seq_len(ncol_aln), state=state, max_posterior=round(maxpost,3),
                       clade_gap_frac=round(gapfrac,3), kept=keep),
            "results/asr/ancestor_confidence.tsv", sep="\t", quote=FALSE, row.names=FALSE)

cat(sprintf("\nfold-ready length: %d aa (dropped %d high-gap columns)\n", nchar(anc_fold), sum(!keep)))
cat(sprintf("mean per-site posterior: %.3f ; sites >=0.8: %.0f%% ; >=0.5: %.0f%%\n",
            mean(maxpost), 100*mean(maxpost>=0.8), 100*mean(maxpost>=0.5)))
cat(sprintf("WEDGE (ScPif1 R324, aln col %d, fold pos %d): ancestor = %s  (posterior %.3f)\n",
            wedge_col, wedge_fold_pos, wedge_state, wedge_post))
cat(sprintf("  -> %s\n", if (wedge_state %in% c("R","K")) "BASIC wedge retained in the ancestor"
                          else paste("ancestor wedge is", wedge_state, "(not basic) -- inspect")))
