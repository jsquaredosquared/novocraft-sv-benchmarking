TRUTH_SET = "/export/home/jeffrey/Documents/sv-benchmarking/resources/sv-benchmarks/HG002/HG002_SVs_Tier1_v0.6.vcf.gz"

size_bins = [(50,100), (100,500), (500,1000), (1000,1_000_000)]

for min_len, max_len in size_bins:
    $(vembrane filter -o @(TRUTH_SET.replace(".vcf", ".{0}-{1}.vcf".format(min_len, max_len))) @('{0}<=INFO["SVLEN"][0]<={1}'.format(min_len,max_len)) @(TRUTH_SET))