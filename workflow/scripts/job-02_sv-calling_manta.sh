snakemake -s 02_call.smk \
--cores all \
--use-conda \
--executor cluster-generic \
--cluster-generic-submit-cmd "qsub -sync y -cwd -pe slave_pe 30 -q slave.q \
-N run-manta" \
--jobs 3