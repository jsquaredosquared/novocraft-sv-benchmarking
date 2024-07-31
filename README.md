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

The benchmarking process has been implemented as a Snakemake workflow (found in `workflow`). The workflow can be configured by editing the `config.yaml` file found in `config`. The workflow performs the following steps:

#### `01_align`

This step aligns the FASTQ reads for each sample listed in the "samples" section of the config file and produces a CRAM file for each aligner listed in the "aligners" section of the config file. This workflow currently works with the following aligners:

- BWA-MEM2
- Novoalign

Other aligners can be added, provided that a samtools-indexed CRAM file is produced.

#### `02_call`

This step takes each CRAM file and calls structural variants using each SV caller listed in the "callers" section of the config file. For each caller, the default settings recommended by the developer were used. This workflow currently works with the following callers:

- Delly
- Dysgu
- Manta

Other callers can be added, provided that a bgzipped and tabix-indexed SV VCF file is produced.

#### `03_benchmark`

This step uses `truvari bench` to compare each SV VCF file to the truth set to calculate the performance characteristics (recall, precision, F1 score). This workflow calculates the overall performance characteristics, as well as the performance characteristics by SVTYPE (DEL, DUP, INS, INV) and SVLEN (50-100, 100-500, 500-1000, 1000+).

The default `truvari bench` settings were used, except for the following:

- Sequence comparison was turned of (`--pctseq`) because the SV VCF is not guaranteed to have sequence-resolved calls.
- The size of SVs was limited using the options `--sizemin 50` and `--sizemax 1_000_000`.
- Only variants with `FILTER == PASS` are considered (`--passonly`).

As per the default settings, SVs are considered the same if:

- they have the same SVTYPE (`--typeignore False`)
- the distance between their breakends is less than 500 bp (`--refdist 500`)
- the size of the smaller SV is at least 70% the size of the larger SV (`--pctsize 0.7`)

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
