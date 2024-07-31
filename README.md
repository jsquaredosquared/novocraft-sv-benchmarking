# Novocraft SV benchmarking

## Introduction

This repository contains a Snakemake workflow used to compare the performance of SV callers when different short-read aligners are used. Currently it is being designed to compare SV calls by Manta on BAM files produced with NovoAlign and BWA-MEM2, although it could be extended to compare other SV callers or aligners.

If NovoAlign performs well, an SV calling pipeline could potentially be incorporated into Novocraft's products.

## Overview

### Datasets

| Sample                                                                                             | Reference                                                                                          | SV truth set                                                                                                                        |
| -------------------------------------------------------------------------------------------------- | -------------------------------------------------------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------- |
| [HG002](https://github.com/human-pangenomics/HG002_Data_Freeze_v1.0) (Illiumina WGS 150 bp PE 30x) | [GIAB GRCh37](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/references/GRCh37/) | [HG002 SVs Tier 1](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NIST_SV_v0.6/) |

### Methods

The process has been implemented as a Snakemake workflow (found in `./workflow`) that performs the following steps:

#### `01_align`

This step aligns the FASTQ reads for each sample listed in the "samples" section of the config file and produces a CRAM file for each aligner listed in the "aligners" section of the config file.

#### `02_call`

This step takes each CRAM file and calls structural variants using each SV caller listed in the "callers" section of the config file. For each caller, the developer's default and/or recommended settings were used. This workflow currently works with the following callers:

- Delly
- Manta

#### `03_benchmark`

This step uses `truvari bench` to compare each SV VCF file to the truth set to determine performance characteristics (recall, precision, F1 score).

#### `04_compare`

This step compares the performance characteristics for each aligner/caller pair.

## Results

### Questions

- Does Novoalign perform better overall?
  - [ ] Calculate recall, precision, and F1 scores (3 s.f.).
    | Aligner   | Caller | F1    | Recall | Precision |
    | --------- | ------ | ----- | ------ | --------- |
    | Novoalign | Manta  | 0.572 | 0.420  | 0.895     |
    | Bwa-mem2  | Manta  | 0.534 | 0.384  | 0.894     |
    | Novoalign | Delly  | 0.431 | 0.295  | 0.795     |
    | Bwa-mem2  | Delly  | 0.371 | 0.246  | 0.753     |
- Does Novoalign perform better for certain types or sizes of SVs?
  - [ ] Group by type/size, calculate performance characteristics for each group, then compare.
  - [ ] Generate Venn diagrams or upset plots to compare the calls made from each caller.

## Acknowledgements

This work was part of my internship at [Novocraft](novocraft.com) (July 2024 - September 2024).

## References
