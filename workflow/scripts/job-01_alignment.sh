snakemake -c 40 -s 01_align.smk \
--use-conda \
--executor cluster-generic \
--cluster-generic-submit-cmd "qsub -N sv-benchmarking_alignment -pe 40 -q slave.q -o /export/home/jeffrey/Documents/sv-benchmarking/logs/job-01_alignment.output.log -e /export/home/jeffrey/Documents/sv-benchmarking/logs/job-01_alignment.error.log -cwd" \
--cluster-generic-status-cmd "qstat -f" \
-j 2