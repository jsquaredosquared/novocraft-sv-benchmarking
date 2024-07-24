# Index with novoindex

~/.local/bin/novocraft/novoindex \
/export/home/jeffrey/Documents/sv-benchmarking/resources/reference-genome/GRCh37/hs37d5.fa.nix \
/export/home/jeffrey/Documents/sv-benchmarking/resources/reference-genome/GRCh37/hs37d5.fa

# Index with samtools

docker run -v /export:/export staphb/samtools samtools faidx /export/home/jeffrey/Documents/sv-benchmarking/resources/reference-genome/GRCh37/hs37d5.fa

# Index with bwa

docker run -v /export:/export dceoy/bwa-mem2 index /export/home/jeffrey/Documents/sv-benchmarking/resources/reference-genome/GRCh37/hs37d5.fa