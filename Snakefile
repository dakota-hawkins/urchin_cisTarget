configfile: "config.yaml"
shell.prefix("module load bedops; module load bedtools; ")

from pathlib import Path

rule all:
    input:
        Path("output").joinpath(config['params']['db_prefix']).with_suffix(".rankings.feather")
#TODO:
# extract genomic reads -5kb + full gene model
# change format_rsat_motifs to write one PWM per file
# write motif names to separate filesrast
#  --- is in a sort of TRANSFAC file format, but BIOPython doesn't like
rule format_tf_motifs:
    input:
        rsat=config['input']['rsat_motifs']
    output:
        motif_dir=directory('output/motifs'),
        motif_ids='output/motif_ids.txt'
    script:
        "scripts/format_transfac_motifs.py"


rule filter_gff_and_bed:
    input:
        gff=config['input']['gff'],
        sequence_ids=config['input']['sequences'],
        chr_size=config['input']['chr_sizes']
    params:
        bp=config['params']['upstream_bp']
    output:
        "output/sequences.bed"
    shell:
        """
        grep -f {input.sequence_ids} {input.gff} | 
        awk '$3=="gene" {{ print }}' | 
        gff2bed | 
        bedtools slop -i stdin -g {input.chr_size} -l {params.bp} -r 0 -s > {output}
        """

rule extract_sequences:
    input:
        bed="output/sequences.bed",
        fasta=config['input']['genome_fasta']
    output:
        "output/selected_fasta.fa"
    shell:
        "bedtools getfasta -fi {input.fasta} -bed {input.bed} -nameOnly > {output}"


rule create_database:
    input:
        fasta="output/selected_fasta.fa",
        motif_dir='output/motifs',
        motif_ids='output/motif_ids.txt'
    params:
        prefix=Path("output").joinpath(config['params']['db_prefix']),
        cbust_loc=config['params']['clusterbuster'],
        cisTarget_loc=config['params']['cisTarget']
    conda:
        "envs/create_cistarget_databases.yaml"
    threads: 16 
    output:
        rankings=Path("output").joinpath(config['params']['db_prefix']).with_suffix(".rankings.feather"),
        scores=Path("output").joinpath(config['params']['db_prefix']).with_suffix(".scores.feather")
    shell:
        """
        {params.cisTarget_loc} -f {input.fasta} -M {input.motif_dir} -m {input.motif_ids} -o {params.prefix} -c {params.cbust_loc} -t {threads}
        """


