# Novocraft SV benchmarking

## Introduction

This repository contains a Snakemake workflow used to compare the performance of SV callers when different short-read aligners are used. Currently it is being designed to compare SV calls by Manta on BAM files aligned with NovoAlign and BWA-MEM2, although it could be extended to compare other SV callers or aligners.

If NovoAlign performs well, an SV-calling pipeline could potentially be incorporated into Novocraft's products.

## Roadmap

- [ ] Align the FASTQ reads using the different aligners.
- [ ] Call the structural variants on the different BAM files using the SV caller(s).
- [ ] Compare the SV calls to the truth set to determine performance characteristics (e.g., F score).
- [ ] Compare the results for each aligner.

## Datasets

Sample: [HG002](https://github.com/human-pangenomics/HG002_Data_Freeze_v1.0) (Whole-genome data, downsampled to ~30x PCR-free Illumina 150bp)
Reference: [GIAB GRCh37](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/references/GRCh37/)
SV truth set: [HG002 SVs Tier 1](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NIST_SV_v0.6/)

## Benchmarking

## Acknowledgements

This work was part of my internship at [Novocraft](novocraft.com) (July 2024 - September 2024).
