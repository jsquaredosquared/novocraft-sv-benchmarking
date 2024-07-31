TRUTH_SET = "resources/sv-benchmarks/HG002/HG002_SVs_Tier1_v0.6.vcf.gz"

for svtype in ["DEL", "DUP", "INS", "INV"]:
    $(vembrane filter -o @(TRUTH_SET.replace(".vcf", f"{svtype}.vcf")) @('INFO["SVTYPE"]=="{0}"'.format(svtype)) @(TRUTH_SET))