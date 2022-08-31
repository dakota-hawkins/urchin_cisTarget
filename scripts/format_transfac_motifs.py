import re
import pathlib
from Bio import motifs


def read_next_line(regex, line):
    return regex.search(line) is None and line != ""


def parse_jaspar_transfac(tf_file):
    dbd_motifs = []
    motif_ids = []
    desc_regex = re.compile("^DE[\s]+")
    id_regex = re.compile("accession\:[\s]+")
    consensus_regex = re.compile("consensus.strict\:[\s]+")
    matrix_regex = re.compile("^[0-9]+[\s]+")
    with open(tf_file, "r") as f:
        line = f.readline()
        while line != "":
            while read_next_line(desc_regex, line):
                line = f.readline()
            desc = desc_regex.split(line)[-1]
            for i in range(2):
                line = f.readline()
            matrix = ""
            while line[:2] != "XX" and line != "":
                matrix += matrix_regex.split(line)[1]
                line = f.readline()
            # skip to CC lines
            line = f.readline()
            # read description lines, check for consensus description.
            while line[:2] == "CC":
                if consensus_regex.search(line) is not None:
                    desc = consensus_regex.split(line)[-1]
                if id_regex.search(line) is not None:
                    motif_ids.append(
                        id_regex.split(line)[-1].replace("\n", "_motif")
                    )
                line = f.readline()
            # edge case for last line of the file
            if line != '':
                dbd_motifs.append(">" + desc + matrix)
    
    assert len(dbd_motifs) == len(motif_ids)
    for motif, name in zip(dbd_motifs, motif_ids):
        with open(motif_dir.joinpath(name).with_suffix('.cb'), "w") as f:
            f.write(motif)
        with open(id_file, "a") as f:
            f.write(name + '\n')

if __name__ == "__main__":
    try:
        snakemake
    except NameError:
        snakemake = None
    if snakemake is not None:
        tf_motifs = pathlib.Path(snakemake.input["rsat"])
        motif_dir = pathlib.Path(snakemake.output["motif_dir"])
        if not motif_dir.exists():
            motif_dir.mkdir()
        id_file = pathlib.Path(snakemake.output["motif_ids"])
        if id_file.exists():
            id_file.unlink()

        extension_to_format = {'.cb': 'clusterbuster', '.jaspar': 'jaspar', '.transfac': "transfac"}  
        if snakemake.params['is_jaspar_transfac'] and tf_motifs.is_file():
            dbd_motifs, motif_ids = parse_jaspar_transfac(tf_motifs)
        elif tf_motifs.is_dir():
            for each in tf_motifs.glob('*'):
                with open(each, 'r') as handle:
                    motif = motifs.read(handle, extension_to_format[each.suffix])
                    # force name to consensus sequence per clusterbuster
                    motif.name = motif.consensus
                    motif_id = each.stem.replace('.', '_')
                    with open(motif_dir.joinpath(motif_id).with_suffix('.cb'), 'w') as handle:
                        handle.write(motif.format('clusterbuster'))
                    with open(id_file, 'a') as handle:
                        handle.write(motif_id + '\n')
        else:
            raise ValueError("Unsupported motif input")
        
