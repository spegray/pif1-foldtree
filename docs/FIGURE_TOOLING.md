# Figure tooling for the PIF1/RRM3 manuscript

*Produced 2026-07-08 by a fan-out research workflow (6 category researchers → synthesis → adversarial
web fact-check → finalize; 9 agents, ~552K tokens). This EXTENDS the G4-RNA-seq figure-as-code stack
(`~/Desktop/G4-RNA-seq/analysis/figures/FIGURE_TOOLCHAIN_PLAN.md`) with the phylogenetics-specific layer.
Load-bearing claims (arm64 availability, license, file-format ingestion, tool currency) were web-verified
where possible — see "Confidence & provenance" at the end for what still needs a local smoke-test.*

---

## Executive summary

The existing G4-RNA-seq figure-as-code stack (native-arm64 R 4.4.3 with **ggtree/treeio/ggplot2/patchwork/svglite**, the `figs` miniforge env with **PyMOL 3.1.0**, and Inkscape/Ghostscript/rsvg for finishing) already covers ~80% of the nine PIF1/RRM3 figures — every tree, bar chart, and structure render lives there. The only phylo-specific gaps are four small, in-discipline additions: **deeptime** (geological time-axis bands for Fig 5, paired with ggtree's own `revts()`), **thirdkind** (recPhyloXML reconciliation render for Fig 4), **graphviz** (reproducible pipeline DAG for Fig 1), and **ggtreeExtra + TreeDist + phytools** for tip-aligned tracks, concordance stats, and tanglegram pre-detangling. A Python-side **ETE toolkit** cross-check (already-installable in the miniforge env) is the one optional add for Fig 4 QC. Philosophy in one line: **reuse the R-native-arm64 ggtree stack for everything, and add only reconciliation + deep-time tooling that ingests files already on disk** — keeping the vector-first / colorblind-safe / font-embedded / headless discipline intact.

## Per-figure tool map

| Figure | What it shows | Recommended tool | Output format | Why |
|---|---|---|---|---|
| **Fig 1** | Pipeline schematic: homologs → structures → corecut → {AA tree, FoldTree, AA+3Di} → species tree → reconciliation | **Graphviz `dot`** → finish in Inkscape | SVG → PDF | A pipeline *is* a DAG; `rankdir=LR` + cluster subgraphs auto-lay-out the 3-tree fork from a git-diffable source, emits real `<text>`. Note: `dot -Tsvg` references fonts by name and does **not** embed them — Fig 1 only passes the pdffonts gate after the Inkscape/Ghostscript conversion embeds them (see gotchas) |
| **Fig 2** | Copy number by subphylum (Saccharomycotina 91/105 two-copy; Pezizomycotina 430/447 one-copy; Taphrinomycotina 11/11; two-copy mushrooms) | **ggplot2 + patchwork** (already installed) | SVG / PDF | Aggregate counts across ~5 clades — a stacked/proportion bar is clearer and more honest than aligning a summary to a 957-tip tree; reads `manifest.csv` directly |
| **Fig 3** | HEADLINE: AA-only vs AA+3Di gene tree, same 957 taxa, dup node moves base-of-Fungi → Saccharomycotina | **ggtree + patchwork** (two collapsed facing panels + `geom_segment` connectors), pre-detangle with **phytools::cophylo** | SVG (svglite) → PDF | Only ggtree collapses the hairball to ~20 legible units, spotlights the one moving clade, and stays in the Okabe-Ito/vector discipline |
| **Fig 4** | Reconciliation: duplication + recurrent single-paralog losses on the species tree (GeneRax; node_448 Saccharomycotina, node_578 mushroom) | **ggtree + treeio** species-tree event map (primary); **thirdkind** for a full-reconciliation supplement / QC cross-check; optional **ETE toolkit** as a second scriptable QC opinion | SVG / PDF | ggtree gives a clean, CVD-safe, single-node-focused figure by joining the event counts onto the species Newick; thirdkind confirms it against the actual recPhyloXML. **Use the labeled grafted event-count file** (see gotchas) — the D=43/D=23 headline numbers are from `aa3di_grafted/run/per_species_event_counts.txt`, verified this session |
| **Fig 5** | Absolute date: duplication on the Shen 2018 chronogram, Devonian–Carboniferous | **ggtree `revts()` + deeptime `coord_geo()`** | SVG / PDF | ggtree's `revts()` reverses the axis to Mya; deeptime then draws authoritative ICS period bands under it via `coord_geo(neg=TRUE)`, staying 100% inside the ggplot/svglite pipeline. `revts()` is a ggtree function, not a deeptime one — call it **before** `coord_geo()` |
| **S1** | Family-filter validation (InterPro IPR048293 vs Pfam PF05970 / Helitron contamination counts) | **ggplot2 + patchwork** | SVG / PDF | Simple grouped bar of counts; no new tooling |
| **S2** | Structure-source + core pLDDT distributions (AFDB vs ColabFold) | **ggplot2** (density/violin split) from `results/reviewer/core_plddt.tsv`; optional **PyMOL** pLDDT-colored core inset | SVG / PDF (+ 600-dpi PNG inset) | Table already has `source`/`core_plddt`; PyMOL adds the concrete structure visual headlessly |
| **S3** | Per-metric FoldTree trees: fident / alntmscore / lddt | **ggtree + ggtreeExtra** (3-panel via patchwork) | SVG / PDF | Same renderer, three trees; `geom_fruit` for any per-tip metric strips |
| **S4** | Sequence-vs-structure concordance | **TreeDist** `ClusteringInfoDistance` for the honest scalar; polished full-tree plot hand-built in **ggtree** | SVG / PDF | Information-theoretic distance behaves correctly at 957 tips where raw RF saturates. Note: TreeDist's `VisualizeMatching` is built to annotate a **handful** of splits and will not render a 957-tip concordance — use it only for the scalar distance and a small illustrative matching, and build the full-tip figure in ggtree |

## The recommended phylo stack

**Tree rendering — ggtree + treeio + ggtreeExtra** (Bioconductor). *Why:* the only code-first, native-arm64, true-vector renderer that parses every file you have (`.treefile` with SH-aLRT/UFBoot, GeneRax `.nhx`, Shen `.newick`). `ggtree`/`treeio` already installed; add `ggtreeExtra`. ggtree also supplies `revts()`, used for the Fig 5 reversed time axis.
`BiocManager::install("ggtreeExtra")` · Artistic-2.0 / GPL-3 · native arm64.

**Tree comparison / tanglegram — phytools (prototype/detangle) + TreeDist (concordance stat).** *Why:* `cophylo(..., rotate=TRUE)` detangles tip order before you collapse; `TreeDist::ClusteringInfoDistance` gives the honest S4 number. Its `VisualizeMatching` companion is for a handful of splits only, so the scalar is the usable output at 957 tips — render the full figure in ggtree.
`install.packages(c("phytools","TreeDist"))` · GPL-2/3 · native arm64.

**Reconciliation — thirdkind** (Rust CLI). *Why:* speaks GeneRax recPhyloXML natively (verified: crate 3.13.5, published 2025-10-24, reads newick/phyloXML/recPhyloXML at reconciliation levels 1/2/3); the `pif1_reconciliated.xml` files are already on disk. Use for the Fig 4 supplement/QC; ggtree-from-NHX is the in-stack primary.
`cargo install thirdkind` · CeCILL-2.1 · builds native aarch64 (needs rustup + a C linker; clang already present).

**Reconciliation cross-check (optional) — ETE toolkit** (`ete3`/`ete4`, Python). *Why:* reads NHX and can lay out gene-in-species reconciliations programmatically inside the existing miniforge Python discipline — a fully-local, vector-capable second opinion alongside thirdkind for Fig 4 QC. Not a primary deliverable; a correctness check.
`mamba install -n figs -c conda-forge ete3` · GPL-3 · pure-Python (arm64-fine).

**Deep-time — deeptime** (CRAN). *Why:* `coord_geo(neg=TRUE)` paints ICS Devonian/Carboniferous bands under the reversed-Mya ggtree axis, fully in-pipeline. deeptime supplies `coord_geo()`/`coord_geo_radial()` and the ICS period data; the axis reversal itself is ggtree's `revts()`.
`install.packages("deeptime")` · GPL(>=3) · native arm64 (pure-R).

**Schematic — Graphviz** (`dot`). *Why:* reproducible auto-laid-out DAG → editable SVG with real text. `dot` does not embed fonts, so the SVG only becomes pdffonts-clean after the Inkscape/Ghostscript step.
`brew install graphviz` · EPL-1.0 · native arm64 bottle.

**Structure panel — PyMOL open-source 3.1.0** (already installed) + optional **foldseek** for the 3Di string. *Why:* headless ray render of the PIF1 core by pLDDT / per-residue 3Di; the only mac headless renderer (ChimeraX `--offscreen` uses OSMesa and is Linux-only, confirmed). foldseek is the sole source of the 3Di alphabet, but `results/seq_tree/3di.trim.fasta` already holds the aligned 3Di, so you almost certainly need no re-run. **If** you do need foldseek on Apple Silicon, the bioconda package has **no osx-arm64 build** — use the project's direct `foldseek-osx-universal.tar.gz` binary (which carries the M1/NEON support) or an osx-64 env under Rosetta; the mamba/bioconda one-liner will fail here (see Install delta).

## Install delta

Everything below is *net-new* (ggtree/treeio/ape/ggplot2/patchwork/PyMOL/matplotlib/Inkscape/gs/rsvg are already present and are not relisted).

```bash
# --- R, CRAN (native-arm64, plain install.packages) ---
/usr/local/bin/Rscript -e 'install.packages(c("deeptime","phytools","TreeDist"), repos="https://cloud.r-project.org")'

# --- R, Bioconductor (use BiocManager, NOT install.packages) ---
/usr/local/bin/Rscript -e 'BiocManager::install("ggtreeExtra")'

# --- Homebrew (native arm64 bottle) ---
brew install graphviz

# --- Rust toolchain + thirdkind (native aarch64; clang linker already present) ---
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y
~/.cargo/bin/cargo install thirdkind          # crate 3.13.5, CeCILL-2.1

# --- OPTIONAL: ETE toolkit in the figs env, for Fig 4 reconciliation QC only ---
mamba install -n figs -c conda-forge ete3

# --- OPTIONAL: foldseek, ONLY if you must regenerate the 3Di ---
# results/seq_tree/3di.trim.fasta already holds the aligned 3Di -> normally SKIP.
# Apple Silicon: bioconda has NO osx-arm64 foldseek build, so `mamba install ... foldseek` FAILS.
# Use the direct universal binary instead:
#   curl -LO https://mmseqs.com/foldseek/foldseek-osx-universal.tar.gz
#   tar xzf foldseek-osx-universal.tar.gz   # ./foldseek/bin/foldseek (arm64/NEON)
# ...or build an osx-64 env and run under Rosetta if you specifically want the conda package.
```

Verify each after install: `dot -V`, `~/.cargo/bin/thirdkind -h`, and in R `packageVersion("deeptime")`. Pin versions in your env lock for reproducibility (3Di states and support-label parsing are version-sensitive).

## Phylo-specific gotchas

- **Large-tree SVG bloat.** A 957-tip rectangular tree in svglite balloons to tens of MB (one DOM node per branch/label/point) and chokes Inkscape. Use `cairo_pdf` for large *final* trees (compact, font-embedded) and reserve svglite for small hand-edited panels; render truly huge supplementary trees as high-dpi raster if needed. (Smoke-test both paths on your actual trees — the crossover point is machine-dependent.)
- **Tip-font floor at 957 tips.** No font size is both readable and non-overlapping. **Never `geom_tiplab` all tips.** Use fan/circular layout, `collapse()`/`scaleClade()` to fold clades to triangles, label only the Pif1/Rrm3 clades with `geom_cladelab`, and encode identity as a `geom_fruit` color strip.
- **Support values as glyphs, not text.** Show SH-aLRT/UFBoot as thresholded node dots — `geom_point2(aes(subset = !isTip & ufboot>=95 & shalrt>=80))` — not printed numbers.
- **IQ-TREE combined labels.** Internal labels are a single `"96.9/95"` string (SH-aLRT/UFBoot). Use `treeio::read.iqtree()` (or `read.newick` + one `tidyr::separate`); plain `ape::read.tree` mishandles the `/`.
- **NHX / recPhyloXML parsing.** `treeio::read.nhx()` parses the GeneRax `[&&NHX:S=…:D=Y/N:B=…]` tags so `D=='Y'` gives duplication nodes for free (verified on-disk: 464 `D=Y` nodes, matching `eventCounts D:464`) — but treeio has **no recPhyloXML reader**, so drive the Fig 4 ggtree map from the species Newick + the event-count table, not the `.xml`.
- **GeneRax has no standalone loss column.** In UndatedDL output a single-paralog loss is encoded as **SL** (speciation-loss); the NHX carries no per-node loss tag at all (verified: `SL:3226`, plain `L:0`). "Recurrent losses" in Fig 4 = the `SL` column; duplications = `D`. Losses must therefore come from the recPhyloXML or the SL column.
- **Fig 4 event-count files disagree — use the labeled grafted file.** `pif1_speciesEventCounts.txt` has **4 unlabeled** columns per node, whereas `per_species_event_counts.txt` is labeled (`#S #SL #D #T #TL`); counts also differ across runs (NCBI `aa3di` vs `aa3di_grafted`). **Resolved this session:** the D=43 (Saccharomycotina, `node_448`) / D=23 (mushroom, `node_578`) headline numbers are the labeled `D=` column of `results/reconciliation/aa3di_grafted/run/per_species_event_counts.txt` (the *grafted* run), verified directly (node_448 is the tree-wide maximum). Build Fig 4 from that file; ignore the unlabeled `pif1_speciesEventCounts.txt` and the NCBI-run counts (which give the lower 9-dup value).
- **cairo_pdf silently rasterizes alpha/gradients.** `geom_hilight` (your Saccharomycotina highlight) uses a semi-transparent fill. Render Figs 3/4 via **svglite → SVG**, convert to font-embedded PDF with Ghostscript/Inkscape, and confirm with `pdffonts`. Do not go straight to `cairo_pdf` on highlighted trees.
- **Graphviz SVG is not self-embedding.** `dot -Tsvg` writes font *references*, not glyphs, so the raw schematic fails `pdffonts`. Fonts get embedded only in the downstream Inkscape/Ghostscript conversion — treat that step as mandatory for Fig 1, not optional polish.
- **thirdkind: arm64 build + bespoke style.** No Homebrew formula; install via cargo (needs rustup, clang linker already present). Its SVG uses its own layout/colors/fonts, so budget an Inkscape restyle to Okabe-Ito, and collapse the species tree (`-Y`) since the literal 940-tip embedding is unreadable as a main figure.
- **TreeDist `VisualizeMatching` does not scale.** Its API and vignette are built to annotate a small number of splits (the docs demo six). At 957 tips it cannot produce the concordance figure — take only the scalar `ClusteringInfoDistance` (and, if useful, a small illustrative matching) from TreeDist and build the full-tip S4 plot in ggtree.
- **Tanglegram connectors.** At 957 tips a full tanglegram is a hairball in every tool — better rotation won't save it. **Collapse first**, then color/weight only the connectors whose destination clade changes (one saturated Okabe-Ito color; everything else thin grey85). Verify the connector color/weight survives a CVD sim.
- **Shen tree units + axis order.** `data/species_tree/shen2018_timetree.newick` branch lengths are in **×100 Myr** (max single edge 2.886, tip edges ~2.1, root ≈ 4.04 → ~404 Mya, consistent with Shen 2018). Multiply `edge.length` by 100 **before** you call ggtree's `revts()`, and call `revts()` **before** `coord_geo(neg=TRUE)`, or the Devonian–Carboniferous band lands in the wrong place. The Shen Newick has no node-age CIs, so mark the event as a shaded band, not a `geom_range` HPD bar.
- **PyMOL renders are raster.** Structure images are ray-traced PNG — embed as 600-dpi insets inside otherwise-vector SVG/PDF (standard for Nature structure panels) and note it in QC. Default PyMOL palette is not Okabe-Ito; set colors explicitly.

## Start here

**Build Fig 2 (copy number by subphylum) first.** It is the highest-value, lowest-risk prototype: it needs **zero new tooling** (ggplot2 + patchwork are installed), exercises your full svglite → Inkscape → pdffonts → CVD gate end to end on real data, and produces a publishable panel in one sitting — a clean win before you wrestle the 957-tip trees. (If you'd rather de-risk the headline early, Fig 3 is the alternative, but it front-loads the collapse/detangle problem.)

Concrete build against real files in `~/Desktop/pif1-foldtree`:

```r
library(tidyverse); library(svglite)
m <- read_csv("manifest.csv")          # has subphylum, structure_source, taxid, family
cn <- m |>
  count(subphylum, taxid, name = "paralogs") |>
  mutate(copies = if_else(paralogs >= 2, "two-copy", "single-copy")) |>
  count(subphylum, copies)
ggplot(cn, aes(subphylum, n, fill = copies)) +
  geom_col(position = "fill") +
  scale_fill_manual(values = c(`single-copy` = "#0072B2", `two-copy` = "#D55E00")) +  # Okabe-Ito
  coord_flip() + theme_classic()
ggsave("docs/figures/fig2_copynumber.svg", width = 89, units = "mm", device = svglite)
```

Then finish and QC exactly as the existing plan prescribes: `inkscape … --export-type=pdf`, `pdffonts fig2.pdf` (every row `emb=yes`, no Type 3), and a CVD sim. Confirm the `subphylum` column values and per-species paralog logic against the manifest before trusting the counts, and reconcile the totals against the caption figures (91/105, 430/447, 11/11).

## Honest tradeoffs

- **iTOL / FigTree / Dendroscope (GUI) vs ggtree (code).** iTOL is genuinely the most polished big-annotated-tree tool and has a batch API, but it is web-based (unpublished data leaves the machine), gates key export behind a subscription, and doesn't fit a fully-local vector/font-embed gate. FigTree v1.4.4 literally can't find Java on Apple Silicon without a manual JDK and is unscriptable; Dendroscope is dormant Java with no confirmed arm64 build. All fail the reproducible/auditable requirement — keep FigTree only as a throwaway eyeball QC viewer.
- **thirdkind vs ggtree-from-NHX vs ETE for Fig 4.** thirdkind draws the *literal* reconciliation (genes-in-species with D/L glyphs) but at 940 tips it's a dense hairball needing aggressive `-Y` collapsing and Inkscape restyling to match house style. The ggtree species-event-map gives a clean, single-node-focused, CVD-safe figure you fully control — so ggtree is the **primary deliverable**, thirdkind the **supplement + correctness cross-check**, and ETE (Python, reads the same NHX) an optional **second independent QC opinion** — not the other way round.
- **Tanglegram vs side-by-side for Fig 3.** A true `cophylo`/`cophyloplot` tanglegram of 957 tips is unreadable no matter the rotation algorithm (the 2025/2026 displacement-optimized methods beat phytools on crossings but don't fit the R/vector stack, and don't matter once you collapse). The winning approach is **two collapsed facing ggtree panels** with the duplication node explicitly marked in each (base-of-Fungi left, within-Saccharomycotina right) and only the moving connectors highlighted — use `phytools::cophylo` merely as a fast prototype to eyeball the topology difference and pre-detangle tip order.
- **dendextend for S4 — rejected.** Its `tanglegram`/`entanglement` aesthetics are tempting, but it operates on `hclust`/dendrograms; coercing an ML gene tree via `chronos()` throws away the branch-length/rooting structure the reconciliation argument rests on and invites reviewer objections. Use **TreeDist** (`ClusteringInfoDistance`) instead — native R, information-theoretic, well-behaved at 957 tips (just remember its `VisualizeMatching` renderer is for a handful of splits, so the full figure comes from ggtree).
- **BioRender for Fig 1 — rejected.** GUI (breaks figure-as-code), and exports are not freely licensed (paid plan to publish; BioRender retains copyright to its icons and to modifications; the Sept-2024 CC-BY claim still has contradictions). Graphviz → Inkscape gives a diff-able, openly-licensed, reproducible schematic.
- **Mermaid vs Graphviz for Fig 1.** Mermaid source is more readable, but its default SVG wraps labels in `<foreignObject>` HTML that won't import cleanly into Inkscape/librsvg (black boxes, missing text) and won't embed fonts — you must set `htmlLabels:false`. Graphviz emits real `<text>` out of the box and its rank-based layout maps directly onto the pipeline topology, so it's the safer default for the QC gate.

## Confidence & provenance

**Web-verified load-bearing claims (safe to rely on as written):**
- **thirdkind** is crate **3.13.5** (published 2025-10-24), license **CeCILL-2.1**, and reads newick/phyloXML/recPhyloXML at reconciliation levels 1/2/3 — verified against the crates.io API and the `simonpenel/thirdkind` GitHub README. cargo (not Homebrew) is the correct install path.
- **ChimeraX `--offscreen`** uses OSMesa and is "currently supported only on Linux" — verified against the UCSF ChimeraX options docs. This is why PyMOL is the headless render path on this Mac.
- **Mermaid** emits `<foreignObject>` HTML labels by default that break in Inkscape/librsvg, and `htmlLabels:false` forces standard `<text>` — verified against Inkscape GitLab #3268, mermaid-cli #112, and librsvg #1050. The Graphviz-vs-Mermaid tradeoff stands.
- **deeptime** supplies `coord_geo()`/`coord_geo_radial()` and the ICS timescale data; the axis-reversal `revts()` is a **ggtree** function — verified against deeptime's "Adding geological timescales to phylogenies" vignette, which instructs calling ggtree's `revts()` first, then `neg=TRUE` in `coord_geo()`.
- **bioconda `foldseek` has no osx-arm64 build** (only linux-64 and osx-64) — verified on anaconda.org/bioconda/foldseek. The arm64/NEON support lives in the separate `foldseek-osx-universal.tar.gz` direct binary.

**Verified against on-disk files (specific to this machine's data):**
- GeneRax NHX tags are `[&&NHX:S=…:D=N:H=N:B=…]` with **464 `D=Y`** nodes, matching `eventCounts D:464`; losses appear only as **`SL` (SL:3226, plain L:0)** — confirming the "no standalone loss column; losses = SL" guidance.
- The Fig 4 headline counts (`node_448` D=43 Saccharomycotina — tree-wide max; `node_578` D=23 mushroom) are the labeled `results/reconciliation/aa3di_grafted/run/per_species_event_counts.txt` (grafted run), verified directly this session.
- `results/seq_tree/3di.trim.fasta` is present, so foldseek regeneration is normally unnecessary.
- `data/species_tree/shen2018_timetree.newick` has 662 branch lengths, max single edge 2.886, tip edges ~2.0–2.1 → root-to-tip depth ~4, i.e. the **×100 Myr** scaling (~400 Mya crown), consistent with Shen 2018.
- `manifest.csv` carries `subphylum`, `taxid`, `family`, and `structure_source`, as the Fig 2 build assumes.

**Remaining uncertain — smoke-test locally before relying:**
- **TreeDist `VisualizeMatching` at scale.** Documented for ~six splits; only the scalar `ClusteringInfoDistance` is safe at 957 tips — confirm the full figure is built in ggtree.
- **svglite vs cairo_pdf crossover** for the large trees, and the whole Fig 3/4 alpha-fill → SVG → font-embedded PDF path, should be run once end-to-end and checked with `pdffonts` before you trust them for the headline figures.
- **Toolchain build success** (`dot -V`, `~/.cargo/bin/thirdkind -h`, `packageVersion("deeptime")`, ETE import) should be verified post-install — the license/version facts are confirmed, but the local arm64 builds are not yet exercised on this machine.
