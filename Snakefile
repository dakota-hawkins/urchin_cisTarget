configfile: "config.yaml"


shell.prefix("module load bedops; module load bedtools; ")

from pathlib import Path


rule all:
    input:
        rvm_rankings=Path("output")
            .joinpath("databases", config["params"]["db_prefix"])
            .with_suffix(".regions_vs_motifs.rankings.feather"),
        rvm_scores=Path("output")
            .joinpath("databases", config["params"]["db_prefix"])
            .with_suffix(".regions_vs_motifs.scores.feather"),
        mvr_scores=Path("output")
            .joinpath("databases", config["params"]["db_prefix"])
            .with_suffix(".motifs_vs_regions.scores.feather"),


rule format_tf_motifs:
    input:
        rsat=config["input"]["transfac_motifs"],
    output:
        motif_dir=directory("output/motifs"),
        motif_ids="output/motif_ids.txt",
    script:
        "scripts/format_transfac_motifs.py"


rule get_chromosome_lengths:
    input:
        fasta=config["input"]["genome_fasta"],
    output:
        "output/sizes.genome",
    shell:
        "faidx {input.fasta} -i chromsizes > {output}"


rule extract_gene_sequences:
    input:
        gff=config["input"]["gff"],
        sequence_ids=config["input"]["sequences"],
        chr_size="output/sizes.genome",
        fasta=config["input"]["genome_fasta"],
    params:
        bp=config["params"]["upstream_bp"],
    output:
        "output/selected_fasta.fa",
    shell:
        """
        grep -f {input.sequence_ids} {input.gff} | 
        awk '$3=="gene" {{ print }}' | 
        gff2bed | 
        bedtools slop -i stdin -g {input.chr_size} -l {params.bp} -r 0 -s |
        bedtools getfasta -fi {input.fasta} -bed stdin -nameOnly > {output}
        """


rule create_database:
    input:
        fasta="output/selected_fasta.fa",
        motif_dir="output/motifs",
        motif_ids="output/motif_ids.txt",
    params:
        prefix=Path("output").joinpath("databases", config["params"]["db_prefix"]),
        cbust_loc=config["params"]["clusterbuster"],
        cisTarget_loc=config["params"]["cisTarget"],
    conda:
        "envs/create_cistarget_databases.yaml"
    threads: 16
    output:
        rvm_rankings=Path("output")
            .joinpath("databases", config["params"]["db_prefix"])
            .with_suffix(".regions_vs_motifs.rankings.feather"),
        rvm_scores=Path("output")
            .joinpath("databases", config["params"]["db_prefix"])
            .with_suffix(".regions_vs_motifs.scores.feather"),
        mvr_scores=Path("output")
            .joinpath("databases", config["params"]["db_prefix"])
            .with_suffix(".motifs_vs_regions.scores.feather"),
    shell:
        """
        {params.cisTarget_loc} -f {input.fasta} -M {input.motif_dir} -m {input.motif_ids} -o {params.prefix} -c {params.cbust_loc} -t {threads}
        """
