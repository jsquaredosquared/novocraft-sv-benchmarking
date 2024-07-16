snakemake -s 01_align.smk \
-c 50 \
--sdm apptainer \
--executor cluster-generic \
--cluster-generic-submit-cmd "qsub -sync y -cwd -pe slave_pe 50 -q slave.q -o /export/home/jeffrey/Documents/sv-benchmarking/logs/job-01_alignment.output.log -e /export/home/jeffrey/Documents/sv-benchmarking/logs/job-01_alignment.error.log -N sv-benchmarking_alignment" \
-j 2