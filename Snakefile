configfile: "config.yaml"

from pathlib import Path

# TODO now multiple gene models for genes / TFs
# do the gene#1 during renaming, add extra motifs + append extra motif -> gene table

def get_kbp(wc):
    lookup = {f"{bp // 1000}kb": bp for bp in config['params']['upstream_bp']}
    return lookup[wc.kbp]

rule all:
    input:
        gvm_rankings=expand(
            Path("output")
            .joinpath("databases", config["params"]["db_prefix"] + '{kbp}')
            .with_suffix(".genes_vs_motifs.rankings.feather"),
            kbp=[f"{bp // 1000}kb" for bp in config['params']['upstream_bp']]
        ),
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
        bp=lambda wc: get_kbp(wc),
    output:
        temp("output/{kbp}/selected_sequences_temp.fa"),
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
        fasta="output/{kbp}/selected_sequences_temp.fa",
        csv=config['input']['id_to_name']
    output:
        fasta="output/{kbp}/selected_sequences.fa"
    script:
        "scripts/rename_sequences.py"


rule create_database:
    input:
        fasta="output/{kbp}/selected_sequences.fa",
        motif_dir="output/motifs",
        motif_ids="output/motif_ids.txt",
    params:
        prefix=Path("output").joinpath(
            "databases", 
            config["params"]["db_prefix"] + "{kbp}"
        ),
        cbust_loc=config["params"]["clusterbuster"],
        cisTarget_loc=config["params"]["cisTarget"],
        regex="#[0-9]+$"
    conda:
        "envs/create_cistarget_databases.yaml"
    threads: 16
    output:
        gvm_rankings=Path("output")
            .joinpath("databases", config["params"]["db_prefix"] + "{kbp}")
            .with_suffix(".genes_vs_motifs.rankings.feather"),
        gvm_scores=Path("output")
            .joinpath("databases", config["params"]["db_prefix"] + "{kbp}")
            .with_suffix(".genes_vs_motifs.scores.feather"),
        mvg_scores=Path("output")
            .joinpath("databases", config["params"]["db_prefix"] + "{kbp}")
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
