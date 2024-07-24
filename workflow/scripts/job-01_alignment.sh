snakemake -s 01_align.smk \
-c 60 \
--sdm conda \
--executor cluster-generic \
--cluster-generic-submit-cmd "qsub -V -sync y -cwd -pe slave_pe 30 -q slave.q -N alignment" \
-j unlimited