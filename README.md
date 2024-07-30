# Novocraft SV benchmarking

## Introduction

This repository contains a Snakemake workflow used to compare the performance of SV callers when different short-read aligners are used. Currently it is being designed to compare SV calls by Manta on BAM files produced with NovoAlign and BWA-MEM2, although it could be extended to compare other SV callers or aligners.

If NovoAlign performs well, an SV-calling pipeline could potentially be incorporated into Novocraft's products.

## Overview

The process has been implemented as a Snakemake workflow (found in `./workflow`) that performs the following steps:

- [x] Align the FASTQ reads using the different aligners.
- [ ] Call the structural variants on the different BAM files using the SV caller(s).
- [ ] Compare the SV calls to the truth set to determine performance characteristics.
- [ ] Compare the results for each aligner.

For each step, the default/recommended settings were used according to the documentation.

## Datasets

- Sample: [HG002](https://github.com/human-pangenomics/HG002_Data_Freeze_v1.0) (Whole-genome data, downsampled to ~30x PCR-free Illumina 150bp)
- Reference: [GIAB GRCh37](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/references/GRCh37/)
- SV truth set: [HG002 SVs Tier 1](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NIST_SV_v0.6/)

## Benchmarking

### Questions

Note: Focus first on insertions and deletions.

- Does Novoalign perform better overall?
  - Calculate recall, precision, and F1 scores (3 d.p.).
    | Aligner/Caller  | F1    | Recall | Precision |
    | --------------- | ----- | ------ | --------- |
    | Novoalign/Manta | 0.572 | 0.420  | 0.895     |
    | Bwa-mem2/Manta  | 0.534 | 0.384  | 0.894     |
- Does Novoalign perform better for certain types or sizes of SVs?
  - Group by type/size, calculate performance characteristics for each group, then compare.
  - Generate Venn diagrams or upset plots to compare the calls made from each caller.

## Acknowledgements

This work was part of my internship at [Novocraft](novocraft.com) (July 2024 - September 2024).

## References
