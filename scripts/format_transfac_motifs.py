import re
import pathlib


def read_next_line(regex, line):
    return regex.search(line) is None and line != ""


if __name__ == "__main__":
    try:
        snakemake
    except NameError:
        snakemake = None
    if snakemake is not None:
        tf_file = snakemake.input["rsat"]
        motif_dir = pathlib.Path(snakemake.output["motif_dir"])
        if not motif_dir.exists():
            motif_dir.mkdir()
        id_file = pathlib.Path(snakemake.output["motif_ids"])
        if id_file.exists():
            id_file.unlink()
        motifs = []
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
                            id_regex.split(line)[-1]
                            .replace("/", "_")
                            .replace(",", ";")
                            .replace("'", "")
                            .replace("\n", "_motif\n")
                        )
                    line = f.readline()
                motifs.append(">" + desc + matrix)

                # line = f.readline()

        for motif, name in zip(motifs, motif_ids):
            with open(motif_dir.joinpath(name + ".txt"), "w") as f:
                f.write(motif)
            with open(id_file, "a") as f:
                f.write(name)
