from lib2to3.pgen2.pgen import generate_grammar
from tokenize import Name
from Bio import SeqIO
import pandas as pd


def set_gene_name(row):
    if row["name"] == "none":
        return row.gene_id
    return row["name"]


if __name__ == "__main__":
    try:
        snakemake
    except NameError:
        snakemake = None
    if snakemake is not None:
        records = list(SeqIO.parse(snakemake.input["fasta"], "fasta"))
        seq_df = pd.read_csv(snakemake.input["csv"])
        name_counter = {x: 0 for x in seq_df.name.unique()}
        for each in records:
            gid = each.id.split("::")[0]
            gene_name = seq_df.loc[gid, "name"]
            name_counter[gene_name] += 1
            each.id = f"{gene_name}#{name_counter[gene_name]}"
            each.name = f"{gene_name}#{name_counter[gene_name]}"
        SeqIO.write(records, snakemake.output["fasta"], "fasta")
