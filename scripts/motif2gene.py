import pandas as pd
from Bio import motifs
import pathlib

HEADER = [
    "#motif_id",
    "motif_name",
    "motif_description",
    "source_name",
    "source_version",
    "gene_name",
    "motif_similarity_qvalue",
    "similar_motif_id",
    "similar_motif_description",
    "orthologous_identity",
    "orthologous_gene_name",
    "orthologous_species",
    "description",
]

ENTRY = [
    "{}",
    "{}",
    "{}",
    "jaspar",
    "1.1",
    "{}",
    "0.000000",
    "None",
    "None",
    "1.000000",
    "None",
    "None",
    "orthology",
]


def get_line(row, motif_dir):
    with open(motif_dir.joinpath(row.matrix_id).with_suffix(".cb")) as mhandle:
        motif = motifs.read(mhandle, "clusterbuster")
    return (
        "\t".join(ENTRY).format(row.matrix_id, row.base_id, motif.consensus, row.gene)
        + "\n"
    )


if __name__ == "__main__":
    try:
        snakemake
    except NameError:
        snakemake = None
    gene_motif_df = pd.read_csv(snakemake.input["motif2gene"])
    motif_dir = pathlib.Path(snakemake.input["motif_dir"])
    outfn = pathlib.Path(snakemake.output["table"])
    lines = [
        get_line(gene_motif_df.loc[idx, :], motif_dir) for idx in gene_motif_df.index
    ]
    print(len(lines))
    with open(outfn, "w") as handle:
        handle.write("\t".join(HEADER) + "\n")
        for each in lines:
            handle.write(each)
