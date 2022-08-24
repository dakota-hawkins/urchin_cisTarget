configfile: "config.yaml"

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