# af3_expansion — new-organism PIF1/RRM3 structures (batch 2, 2026-07-08)

248 unique-organism PIF1-family sequences (from A_unique_organisms_Need_Structures.fasta)
folded locally with localcolabfold (ColabFold 1.6.1 = AlphaFold2) on an RTX 4090 via WSL2.

- Predictor / source: AF2_ColabFold (remote MMseqs2 MSA)
- Settings: --num-models 1 --num-recycle 3 --num-seeds 1 --model-type alphafold2_ptm
  (>1500 aa also --max-msa 512:1024). Memory-safe: TF_FORCE_UNIFIED_MEMORY=0.

## Folded: 238 / 248  (this folder, one <ACCESSION>.pdb each)

## Excluded: 10 / 248 — GPU-model segfault on their MSA (localcolabfold 1.6.1 + driver 610.62).
## 9 are the same recurring 192-aa fragment across genomes; 1 is a 111-aa fragment. Not modeled:
  - A0A9P7MFY7	111 aa
  - A0A3L6MPD3	192 aa
  - A0A3L6MQ31	192 aa
  - A0A3L6MRC3	192 aa
  - A0A3L6MTF2	192 aa
  - A0A8H6G863	192 aa
  - A0A8H6LLN3	192 aa
  - A0A8H6LNV8	192 aa
  - A0A8J5UJQ7	192 aa
  - W9HE65	192 aa
