novoindex = config["novoindex"]


rule download_reference:
    output:
        "resources/reference-genome/GRCh37/hs37d5.fa"
    log:
        "logs/download_giab_grch37_ref.log",
    shell:
        "(wget ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/references/GRCh37/hs37d5.fa.gz -O {output}.gz "
        "&& gunzip {output}.gz "
        ")2> {log}"


rule novoindex_reference:
    input:
        "resources/reference-genome/GRCh37/hs37d5.fa",
    output:
        "resources/reference-genome/GRCh37/hs37d5.fa.nix",
    log:
        "logs/novoindex_giab_grch37_ref.log",
    shell:
        "{config[novoindex]} {output} {input} 2> {log}"


rule samtools_index_reference:
    input:
        "resources/reference-genome/GRCh37/hs37d5.fa",
    output:
        "resources/reference-genome/GRCh37/hs37d5.fa.fai",
    log:
        "logs/samtools_index_giab_grch37_ref.log",
    conda:
        "../envs/alignment.yaml"
    shell:
        "samtools faidx {input} 2> {log}"


rule bwamem_index_reference:
    input:
        "resources/reference-genome/GRCh37/hs37d5.fa",
    output:
        multiext(
            "resources/reference-genome/GRCh37/hs37d5.fa",
            ".0123",
            ".amb",
            ".ann",
            ".bwt.2bit.64",
            ".pac",
        ),
    log:
        "logs/bwamem_index_giab_grch37_ref.log",
    conda:
        "../envs/alignment.yaml"
    shell:
        "bwa-mem2 index {input} 2> {log}"


rule download_hg002_fastqs:
    output:
        "resources/samples/HG002/HG002_HiSeq30x_subsampled_R1.fastq.gz",
        "resources/samples/HG002/HG002_HiSeq30x_subsampled_R2.fastq.gz",
    log:
        "logs/download_hg002_fastqs.log",
    shell:
        "(wget https://s3-us-west-2.amazonaws.com/human-pangenomics/NHGRI_UCSC_panel/HG002/hpp_HG002_NA24385_son_v1/ILMN/downsampled/HG002_HiSeq30x_subsampled_R1.fastq.gz -O {output[0]} "
        "&& wget https://s3-us-west-2.amazonaws.com/human-pangenomics/NHGRI_UCSC_panel/HG002/hpp_HG002_NA24385_son_v1/ILMN/downsampled/HG002_HiSeq30x_subsampled_R2.fastq.gz -O {output[1]} "
        ")2> {log}"


rule download_hg002_tier1_sv_truth_set:
    output:
        multiext(
            "resources/sv-benchmarks/HG002/HG002_SVs_Tier1_v0.6.ALL",
            ".vcf.gz",
            ".vcf.gz.tbi",
            ".bed",
        ),
    log:
        "logs/download_hg002_truth_set.log",
    shell:
        "(wget ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NIST_SV_v0.6/HG002_SVs_Tier1_v0.6.vcf.gz -O {output[0]} "
        "&& wget ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NIST_SV_v0.6/HG002_SVs_Tier1_v0.6.vcf.gz.tbi -O {output[1]} "
        "&& wget ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NIST_SV_v0.6/HG002_SVs_Tier1_v0.6.bed -O {output[2]} "
        ")2> {log}"


rule download_delly_exclude:
    output:
        "resources/delly/human.hg19.excl.tsv",
    log:
        "logs/download_delly_exclude.log",
    shell:
        "wget https://github.com/dellytools/delly/raw/main/excludeTemplates/human.hg19.excl.tsv -O {output} 2> {log}"


# rule download_all:
#     input:
#         rules.download_reference.output,
#         rules.novoindex_reference.output,
#         rules.samtools_index_reference.output,
#         rules.bwamem_index_reference.output,
#         rules.download_hg002_fastqs.output,
#         rules.download_hg002_tier1_sv_truth_set.output,
#         rules.download_delly_exclude.output,
#     default_target: True
