configfile: "config.yaml"

from pathlib import Path


rule all:
    input:
        gvm_rankings=Path("output")
            .joinpath("databases", config["params"]["db_prefix"])
            .with_suffix(".genes_vs_motifs.rankings.feather"),
        gvm_scores=Path("output")
            .joinpath("databases", config["params"]["db_prefix"])
            .with_suffix(".genes_vs_motifs.scores.feather"),
        mvg_scores=Path("output")
            .joinpath("databases", config["params"]["db_prefix"])
            .with_suffix(".motifs_vs_genes.scores.feather"),
        motif_table='output/cisTarget_motif2tf.tbl'


rule format_tf_motifs:
    input:
        rsat=config["input"]["motifs"],
    output:
        motif_dir=directory("output/motifs"),
        motif_ids="output/motif_ids.txt",
    params:
        is_jaspar_transfac=False
    script:
        "scripts/format_transfac_motifs.py"


rule get_chromosome_lengths:
    input:
        fasta=config["input"]["genome_fasta"],
    output:
        "output/sizes.genome",
    conda:
        "envs/create_cistarget_databases.yaml"
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
        temp("output/selected_sequences_temp.fa"),
    conda:
        "envs/create_cistarget_databases.yaml"
    shell:
        """
        grep -f {input.sequence_ids} {input.gff} | 
        awk '$3=="gene" {{ print }}' | 
        gff2bed | 
        bedtools slop -i stdin -g {input.chr_size} -l {params.bp} -r 0 -s |
        bedtools getfasta -fi {input.fasta} -bed stdin -name > {output}
        """

rule rename_sequences:
    input:
        fasta="output/selected_sequences_temp.fa",
        csv=config['input']['id_to_name']
    output:
        fasta="output/selected_sequences.fa"
    script:
        "scripts/rename_sequences.py"
    

rule create_database:
    input:
        fasta="output/selected_sequences.fa",
        motif_dir="output/motifs",
        motif_ids="output/motif_ids.txt",
    params:
        prefix=Path("output").joinpath("databases", config["params"]["db_prefix"]),
        cbust_loc=config["params"]["clusterbuster"],
        cisTarget_loc=config["params"]["cisTarget"],
        regex=config["params"]['gene_regex']
    conda:
        "envs/create_cistarget_databases.yaml"
    threads: 16
    output:
        gvm_rankings=Path("output")
            .joinpath("databases", config["params"]["db_prefix"])
            .with_suffix(".genes_vs_motifs.rankings.feather"),
        gvm_scores=Path("output")
            .joinpath("databases", config["params"]["db_prefix"])
            .with_suffix(".genes_vs_motifs.scores.feather"),
        mvg_scores=Path("output")
            .joinpath("databases", config["params"]["db_prefix"])
            .with_suffix(".motifs_vs_genes.scores.feather"),
    shell:
        "{params.cisTarget_loc} -f {input.fasta} -M {input.motif_dir} "
        "-m {input.motif_ids} -o {params.prefix} -c {params.cbust_loc} "
        "-t {threads} -g \"{params.regex}\""

rule create_motif2tf_table:
    input:
        motif_dir='output/motifs/',
        motif2gene=config['input']['motif_to_tf']
    output:
        table='output/cisTarget_motif2tf.tbl'
    script:
        "scripts/motif2gene.py"
