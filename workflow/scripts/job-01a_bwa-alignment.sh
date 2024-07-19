docker run -v /export:/export dceoy/bwa-mem2 mem -t 56 \
/export/home/jeffrey/Documents/sv-benchmarking/resources/reference-genome/GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta \
/export/Projects/2024-SVCalling/GIAB-AltScaffold/V300029232_L04_read_1.fq.gz \
/export/Projects/2024-SVCalling/GIAB-AltScaffold/V300029232_L04_read_2.fq.gz \
> /export/home/jeffrey/Documents/sv-benchmarking/resources/bam-files/HG001.bwa-mem.sam \
2> /export/home/jeffrey/Documents/sv-benchmarking/logs/align_HG001_with_bwa-mem.log
