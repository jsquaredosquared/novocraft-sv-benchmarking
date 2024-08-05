TRUTH_SET = $ARG1

for svtype in $ARG2.split(','):
    out_file = TRUTH_SET.replace("ALL", svtype)
    $(vembrane filter @('INFO["SVTYPE"]=="{0}"'.format(svtype)) @(TRUTH_SET) | bgzip -o @(out_file) -)
    tabix -p vcf index @(out_file)