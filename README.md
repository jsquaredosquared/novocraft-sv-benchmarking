# Novocraft SV benchmarking

## Introduction

This repository contains a Snakemake workflow used to compare the performance of SV callers when different short-read aligners are used. Originally it was designed to compare SV calls by Manta on BAM files produced with NovoAlign and BWA-MEM2, although it can be extended to compare other SV callers or aligners.

If NovoAlign performs well, an SV calling pipeline could potentially be incorporated into Novocraft's products.

## Datasets

| Sample                                                                                             | Reference                                                                                          | SV truth set                                                                                                                        |
| -------------------------------------------------------------------------------------------------- | -------------------------------------------------------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------- |
| [HG002](https://github.com/human-pangenomics/HG002_Data_Freeze_v1.0) (Illumina WGS 150 bp PE 30x) | [GIAB GRCh37](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/references/GRCh37/) | [HG002 SVs Tier 1](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NIST_SV_v0.6/) |

## Results

### Questions

- Does Novoalign perform better overall?
  | Aligner   | Caller | F1    | Recall | Precision |
  | --------- | ------ | ----- | ------ | --------- |
  | NovoAlign | Manta  | 0.572 | 0.420  | 0.895     |
  | BWA-MEM2  | Manta  | 0.534 | 0.384  | 0.894     |
  | NovoAlign | Dysgu  | 0.512 | 0.419  | 0.660     |
  | BWA-MEM2  | Dysgu  | 0.508 | 0.392  | 0.723     |
  | NovoAlign | Delly  | 0.431 | 0.295  | 0.795     |
  | BWA-MEM2  | Delly  | 0.371 | 0.246  | 0.753     |
- Does Novoalign perform better for certain types or sizes of SVs?
  - [ ] Group by type/size, calculate performance characteristics for each group, then compare.
  - [ ] Generate Venn diagrams or upset plots to compare the calls made from each caller.

## Methods

The benchmarking process has been implemented as a Snakemake workflow (found in `workflow/rules`). The workflow can be configured by editing the `config/config.yaml` file. Tools can be added if they are available in the working environment or can be accessed using Snakemake's conda or docker integrations. The command-line options for each tool can be modified according to your needs.

The default workflow performs the following steps:

### `01_align`

This step aligns the FASTQ reads for each sample listed in the "samples" section of the config file and produces a CRAM file for each aligner listed in the "aligners" section of the config file. This workflow currently works with the following aligners:

- [x] [BWA-MEM2](https://github.com/bwa-mem2/bwa-mem2)
- [x] [NovoAlign](https://www.novocraft.com/products/novoalign/)

Other aligners can be added, provided that a samtools-indexed CRAM file is produced.

The location of the reference genome can be specified in the config file. The various indexes required by the aligners must be produced using the relevant tools (e.g., `samtools index`, `novoindex`, `bwa-mem2 index`) and made available in the same directory.

### `02_call`

This step takes each CRAM file and calls structural variants using each SV caller listed in the "callers" section of the config file. For each caller, the default settings recommended by the developer were used. This workflow currently works with the following callers:

- [x] [Delly](https://github.com/dellytools/delly)
- [x] [Dysgu](https://github.com/kcleal/dysgu)
- [ ] Lumpy (via Smoove)
- [x] [Manta](https://github.com/Illumina/manta)
- [ ] TIDDIT

Other callers can be added, provided that they accept a CRAM file and a bgzipped, tabix-indexed SV VCF file is produced.

The output is generally in the form `outputs/{caller}/{sample}.{aligner}.{caller}.vcf.gz`.

### `03_benchmark`

This step uses `truvari bench` to compare each SV VCF file to the truth set to calculate the performance characteristics (recall, precision, F1 score). This workflow calculates the overall performance characteristics, as well as the performance characteristics by SVTYPE.

The default `truvari bench` settings were used, with the following exceptions:

- Sequence comparison was turned of (`--pctseq`) because the SV VCF is not guaranteed to have sequence-resolved calls.
- The size of SVs was limited using the options `--sizemin 50` (common definition of SV) and `--sizemax 1_000_000` (because larger events are more likely to be erroneous).
- Only variants with `FILTER == PASS` are considered (`--passonly`).

As per the default settings, 2 SVs are considered the same if:

- they have the same SVTYPE (`--typeignore False`)
- the distance between their breakends is less than 500 bp (`--refdist 500`)
- the size of the smaller SV is at least 70% the size of the larger SV (`--pctsize 0.7`)

The output is generally in the form `outputs/truvari/{sample}.{aligner}.{caller}[.{svtype|svlen}].truvari-bench.json`.

### `04_compare`

This step compares the performance characteristics for each aligner/caller pair.

## Acknowledgements

This work was part of my internship at [Novocraft](novocraft.com) (July 2024 - September 2024).

## References
