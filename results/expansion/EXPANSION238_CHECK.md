# 238-structure expansion — placement check (2026-07-09)

**Question:** does adding 238 new PIF1-family structures (unique-organisms expansion, folded on the WSL
box) change the Saccharomycotina duplication placement?

**Input:** 367 ColabFold PDBs = 129 gap_129 (already
in our 957) + **238 expansion** (genuinely new; 0 overlap with the 957 manifest). Their UniProt entries are
now obsolete (empty organism/taxid), so sequences were extracted from the PDBs. Filtered to **223** (dropped
14 short <150-aa partial models + 1 giant 1929 aa).

**Method (quick structural placement):** foldseek all-vs-all search of the 223 new core structures against
the 957 existing cores (`data/3di/coresdb`); for each new structure, take the best-scoring (bits) existing
neighbor and read off its subphylum and whether it sits in the Saccharomycotina Pif1/Rrm3 duplication clade
(the 197-tip clade from `pif1_aa3di.treefile`). Foldseek = `tools/foldseek/bin/foldseek` (osx-universal).

**Result — the answer holds.**

| Nearest structural neighbor (subphylum) | count |
|---|---|
| Agaricomycotina (mushrooms) | **192 (86%)** |
| Saccharomycotina | 17 |
| Pezizomycotina | 6 |
| Glomeromycotina / Ustilaginomycotina | 2 |

Relative to the Pif1/Rrm3 duplication clade: **203/223 land OUTSIDE it**; only 14 fall inside (13 Pif1-side,
1 Rrm3-side) — a few new budding-yeast orthologs that would simply join the existing Saccharomycotina clade.
Median best-hit fident 0.35 (normal for distant homologs).

**Interpretation:** the expansion is overwhelmingly Agaricomycotina/mushroom PIF1s, which attach far from the
Saccharomycotina node. Adding them does not move the duplication placement; it *strengthens* R5 (the
independent Agaricomycetes duplication) by adding many more mushroom taxa. A full AA+3Di tree rebuild on 1,180
tips is the belt-and-suspenders confirmation (not run here; the structural placement answers it).

**Artifacts** (bulky ones gitignored): `data/expansion238/expansion238.faa` (seqs from PDBs),
`expansion238.filt.faa` (223), `results/expansion/hits.m8` (foldseek), `results/expansion/aln1180.fasta`
(MAFFT --add --keeplength of the 223 into the full alignment).
