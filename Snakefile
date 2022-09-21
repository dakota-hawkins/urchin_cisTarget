configfile: "config.yaml"


from pathlib import Path

# TODO now multiple gene models for genes / TFs
# do the gene#1 during renaming, add extra motifs + append extra motif -> gene table


def get_kbp(wc):
    lookup = {f"{bp // 1000}kb": bp for bp in config["params"]["upstream_bp"]}
    return lookup[wc.kbp]


OUTDIR = Path(config["output"]["dir"])


rule all:
    input:
        gvm_rankings=expand(
            Path("output")
            .joinpath("databases", config["params"]["db_prefix"] + "{kbp}")
            .with_suffix(".genes_vs_motifs.rankings.feather"),
            kbp=[f"{bp // 1000}kb" for bp in config["params"]["upstream_bp"]],
        ),
        motif_table=OUTDIR.joinpath("cisTarget_motif2tf.tbl"),


rule get_chromosome_lengths:
    input:
        fasta=config["input"]["genome_fasta"],
    output:
        OUTDIR.joinpath("sizes.genome"),
    conda:
        "envs/create_cistarget_databases.yaml"
    shell:
        "faidx {input.fasta} -i chromsizes > {output}"


# It was suggested to get N kb upstream, 250bp ds of TSS, and possibly first intron,
# right now we are getting N kb upstream, and the entire gene model including exons
rule extract_gene_sequences:
    input:
        gff=config["input"]["gff"],
        gene_table=config["input"]["id_to_name"],
        chr_size=OUTDIR.joinpath("sizes.genome"),
        fasta=config["input"]["genome_fasta"],
    params:
        up_kb=lambda wc: get_kbp(wc),
        ds_tss_bp=250,
    output:
        upstream=temp(OUTDIR.joinpath("{kbp}", "upstream.fa")),
        selected_bed=temp(OUTDIR.joinpath("{kbp}", "selected.bed")),
        gene_bed=temp(OUTDIR.joinpath("{kbp}", "gene.bed")),
        first_introns=temp(OUTDIR.joinpath("{kbp}", "first_introns.fa")),
    conda:
        "envs/create_cistarget_databases.yaml"
    shell:
        """
        tail -n +2 {input.gene_table} | cut -d',' -f1 |
        grep -f - {input.gff} | 
        gff2bed > {output.selected_bed}; 
        awk '$8 == "gene"' {output.selected_bed} > {output.gene_bed}; 
        bedtools flank -i {output.gene_bed} -g {input.chr_size} -l {params.up_kb} -r 0 -s |
        bedtools slop -i stdin -g {input.chr_size} -l 0 -r {params.ds_tss_bp} -s |
        bedtools getfasta -fi {input.fasta} -bed stdin -name > {output.upstream}
        awk '$8 == "exon"' {output.selected_bed} | 
        bedtools subtract -a {output.gene_bed} -b stdin | 
        awk 'n[$4]++<1' | 
        bedtools getfasta -fi {input.fasta} -bed stdin -name > {output.first_introns}
        """


rule rename_and_combine_sequences:
    input:
        upstream=OUTDIR.joinpath("{kbp}", "upstream.fa"),
        introns=OUTDIR.joinpath("{kbp}", "first_introns.fa"),
        csv=config["input"]["id_to_name"],
    output:
        fasta=OUTDIR.joinpath("{kbp}", "selected_sequences.fa"),
    script:
        "scripts/rename_sequences.py"


def get_db_prefix(wc):
    out = Path("output").joinpath("databases", config["params"]["db_prefix"] + wc.kbp)
    return out


rule extract_motif_ids:
    input:
        config["input"]["motif_to_tf"],
    output:
        OUTDIR.joinpath("motif_ids.txt"),
    shell:
        "tail -n +2 {input} | cut -d',' -f1 > {output}"


rule create_database:
    input:
        fasta=OUTDIR.joinpath("{kbp}", "selected_sequences.fa"),
        motif_dir=config["input"]["motif_dir"],
        motif_ids=OUTDIR.joinpath("motif_ids.txt"),
    params:
        prefix=lambda wc: get_db_prefix(wc),
        cbust_loc=config["params"]["clusterbuster"],
        cisTarget_loc=config["params"]["cisTarget"],
        regex="#[0-9].*$",
    conda:
        "envs/create_cistarget_databases.yaml"
    threads: 16
    output:
        gvm_rankings=OUTDIR.joinpath(
            "databases", config["params"]["db_prefix"] + "{kbp}"
        ).with_suffix(".genes_vs_motifs.rankings.feather"),
        gvm_scores=OUTDIR.joinpath(
            "databases", config["params"]["db_prefix"] + "{kbp}"
        ).with_suffix(".genes_vs_motifs.scores.feather"),
        mvg_scores=OUTDIR.joinpath(
            "databases", config["params"]["db_prefix"] + "{kbp}"
        ).with_suffix(".motifs_vs_genes.scores.feather"),
    shell:
        "{params.cisTarget_loc} -f {input.fasta} -M {input.motif_dir} "
        "-m {input.motif_ids} -o {params.prefix} -c {params.cbust_loc} "
        '-t {threads} -g "{params.regex}"'


rule create_motif2tf_table:
    input:
        motif_dir=config["input"]["motif_dir"],
        motif2gene=config["input"]["motif_to_tf"],
    output:
        table=OUTDIR.joinpath("cisTarget_motif2tf.tbl"),
    script:
        "scripts/motif2gene.py"
