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


# NOTE: No multithreading (https://github.com/dellytools/delly/issues/268#issuecomment-975385454).
rule run_delly:
    input:
        multiext("resources/alignment-files/{sample}.{aligner}", ".cram", ".cram.crai"),
        "resources/delly/human.hg19.excl.tsv",
    output:
        "outputs/delly/{sample}.{aligner}.delly.vcf",
    conda:
        "../envs/delly.yaml"
    log:
        "logs/{sample}.{aligner}.delly.log",
    threads: 1
    shell:
        "delly call -g {config[reference]} "
        "-x {input[2]} "
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
        "--outdir outputs/lumpy "
        "{input[0]} "
        "&& gunzip --stdout outputs/lumpy/{wildcards.sample}-smoove.genotyped.vcf.gz > {output} "
        ")2> {log}"


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

# TODO: Since this only calls insertions, plot generation will need to be modified to accommodate this.
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

        mkdir {params.workdir}
        
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
