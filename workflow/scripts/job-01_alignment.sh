snakemake -s 01_align.smk \
-c 60 \
--executor cluster-generic \
--cluster-generic-submit-cmd "qsub -sync y -cwd -pe slave_pe 30 -q slave.q \
-o /export/home/jeffrey/Documents/sv-benchmarking/logs/job-01_alignment.output.log \
-e /export/home/jeffrey/Documents/sv-benchmarking/logs/job-01_alignment.error.log -N sv-benchmarking_alignment \
-M jeffrey@novocraft.com \
-m beas" \
-j 2
