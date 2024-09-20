rule configure_manta:
    input:
        multiext("resources/alignment-files/{sample}.{aligner}", ".cram", ".cram.crai"),
    output:
        "outputs/manta/{sample}/{aligner}/runWorkflow.py",
    conda:
        "../envs/manta.yaml"
    log:
        "logs/{sample}_{aligner}.configure_manta.log",
    params:
        out_dir=lambda wildcards: f"outputs/manta/{wildcards.sample}/{wildcards.aligner}",
    shell:
        "configManta.py "
        "--bam {input[0]} "
        "--referenceFasta {config[reference]} "
        "--runDir {params.out_dir} "
        "2> {log}"


# TODO Use script to reformat inversions into single fields
rule run_manta:
    input:
        "outputs/manta/{sample}/{aligner}/runWorkflow.py",
    output:
        "outputs/manta/{sample}.{aligner}.manta.vcf",
    conda:
        "../envs/manta.yaml"
    log:
        "logs/{sample}_{aligner}.execute_manta.log",
    threads: 16
    shell:
        "(python {input} "
        "-j {threads} "
        "-e jeffrey@novocraft.com "
        "&& gunzip --stdout outputs/manta/{wildcards.sample}/{wildcards.aligner}/results/variants/diploidSV.vcf.gz > {output} "
        ")2> {log}"


rule run_dysgu:
    input:
        multiext("resources/alignment-files/{sample}.{aligner}", ".cram", ".cram.crai"),
    output:
        "outputs/dysgu/{sample}.{aligner}.dysgu.vcf",
    conda:
        "../envs/dysgu.yaml"
    log:
        "logs/{sample}.{aligner}.dysgu.log",
    threads: 16
    params:
        temp_dir="outputs/dysgu/{sample}.{aligner}.temp",
    shell:
        "dysgu run --clean -p {threads} "
        "{config[reference]} "
        "{params.temp_dir} "
        "{input[0]} "
        "> {output} "
        "2> {log}"


# No multithreading (https://github.com/dellytools/delly/issues/268#issuecomment-975385454).
rule run_delly:
    input:
        multiext("resources/alignment-files/{sample}.{aligner}", ".cram", ".cram.crai"),
    output:
        "outputs/delly/{sample}.{aligner}.delly.vcf",
    conda:
        "../envs/delly.yaml"
    log:
        "logs/{sample}.{aligner}.delly.log",
    threads: 1
    params:
        exclude_regions="resources/delly/human.hg19.excl.tsv",
    shell:
        "delly call -g {config[reference]} "
        "-x {params.exclude_regions} "
        "{input[0]} > {output} "
        "2> {log}"


# TODO: This didn't work. What was the problem again?
rule run_smoove:
    input:
        multiext("resources/alignment-files/{sample}.{aligner}", ".cram", ".cram.crai"),
    output:
        "outputs/lumpy/{sample}.{aligner}.lumpy.vcf",
    conda:
        "../envs/smoove.yaml"
    log:
        "logs/{sample}.{aligner}.lumpy.log",
    params:
        exclude_regions="resources/lumpy/ceph18.b37.lumpy.exclude.2014-01-15.bed",
    threads: 16
    shell:
        "(smoove call "
        "--name {wildcards.sample}.{wildcards.aligner} "
        "--genotype "
        "--removepr "
        "--duphold "
        "--exclude {params.exclude_regions} "
        "--fasta {config[reference]} "
        "-p {threads} "
        "--outdir ../../outputs/lumpy "
        "{input[0]} "
        "&& gunzip --stdout ../../outputs/lumpy/{wildcards.sample}-smoove.genotyped.vcf.gz > {output} "
        ")2> {log}"


# WARNING: WHAM only works with BAMs?
# TODO: Should you use the filter script provided by the devs?
rule run_wham:
    input:
        multiext("resources/alignment-files/{sample}.{aligner}", ".cram", ".cram.crai"),
    output:
        "outputs/wham/{sample}.{aligner}.wham.vcf",
    conda:
        "../envs/wham.yaml"
    log:
        "logs/{sample}.{aligner}.wham.log",
    params:
        exclude="GL000207.1,GL000226.1,GL000229.1,GL000231.1,GL000210.1,GL000239.1,GL000235.1,GL000201.1,GL000247.1,GL000245.1,GL000197.1,GL000203.1,GL000246.1,GL000249.1,GL000196.1,GL000248.1,GL000244.1,GL000238.1,GL000202.1,GL000234.1,GL000232.1,GL000206.1,GL000240.1,GL000236.1,GL000241.1,GL000243.1,GL000242.1,GL000230.1,GL000237.1,GL000233.1,GL000204.1,GL000198.1,GL000208.1,GL000191.1,GL000227.1,GL000228.1,GL000214.1,GL000221.1,GL000209.1,GL000218.1,GL000220.1,GL000213.1,GL000211.1,GL000199.1,GL000217.1,GL000216.1,GL000215.1,GL000205.1,GL000219.1,GL000224.1,GL000223.1,GL000195.1,GL000212.1,GL000222.1,GL000200.1,GL000193.1,GL000194.1,GL000225.1,GL000192.1,NC_007605",
    threads: 10
    shell:
        "whamg -f {input[0]} "
        "-a {config[reference]} "
        "-e {params.exclude} "
        "-x {threads} "
        "> {output} "
        "2> {log} "


# TODO: Find solution to error (https://github.com/SciLifeLab/TIDDIT/issues/109)
rule run_tiddit:
    input:
        multiext("resources/alignment-files/{sample}.{aligner}", ".cram", ".cram.crai"),
    output:
        "outputs/tiddit/{sample}.{aligner}.tiddit.vcf",
    conda:
        "../envs/tiddit.yaml"
    log:
        "logs/{sample}.{aligner}.tiddit.log",
    threads: 8
    shell:
        "tiddit --sv --bam {input[0]} "
        "-o outputs/tiddit/{wildcards.sample}.{wildcards.aligner}.tiddit "
        "--ref {config[reference]} "
        "--threads {threads} "
        "2> {log}"


rule run_insurveyor:
    input:
        multiext("resources/alignment-files/{sample}.{aligner}", ".cram", ".cram.crai"),
    output:
        fixed_cram="resources/alignment-files/{sample}.{aligner}.picard-fmi.cram",
        cram_index="resources/alignment-files/{sample}.{aligner}.picard-fmi.cram.crai",
        insurveyor_vcf="outputs/insurveyor/{sample}.{aligner}.insurveyor.vcf",
    params:
        workdir=lambda wildcards: f"outputs/insurveyor/{wildcards.sample}.{wildcards.aligner}",
    conda:
        "../envs/insurveyor.yaml"
    log:
        "logs/{sample}.{aligner}.insurveyor.log",
    threads: 16
    shell:
        """
        picard FixMateInformation -I {input[0]} -O {output.fixed_cram} -R {config[reference]} --CREATE_INDEX true
        
        insurveyor.py --threads {threads} {output.fixed_cram} {params.workdir} {config[reference]}
        
        gunzip {params.workdir}/out.pass.vcf.gz > {output.insurveyor_vcf}
        """


rule bgzip_and_index_sv_vcf:
    input:
        "outputs/{caller}/{sample}.{aligner}.{caller}.vcf",
    output:
        multiext("outputs/{caller}/{sample}.{aligner}.{caller}.vcf", ".gz", ".gz.tbi"),
    conda:
        "../envs/tabix.yaml"
    log:
        "logs/{sample}.{aligner}.{caller}.bgzip_tabix.log",
    shell:
        "(bgzip {input} && tabix -p vcf {output[0]})2> {log}"
