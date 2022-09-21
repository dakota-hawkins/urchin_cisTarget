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
        upstream_seq = list(SeqIO.parse(snakemake.input["upstream"], "fasta"))
        introns = {
            x.id.split("::")[0]: x
            for x in SeqIO.parse(snakemake.input["introns"], "fasta")
        }
        seq_df = pd.read_csv(snakemake.input["csv"])
        seq_df["name"] = seq_df.apply(lambda x: set_gene_name(x), axis=1)
        seq_df.set_index("gene_id", inplace=True)
        name_counter = {x: 0 for x in seq_df["name"].unique()}
        for each in upstream_seq:
            gid = each.id.split("::")[0]
            gene_name = seq_df.at[gid, "name"]
            name_counter[gene_name] += 1
            if gid in introns:
                each.seq += introns[gid].seq
            each.id = f"{gene_name}#{name_counter[gene_name]}"
            each.name = f"{gene_name}#{name_counter[gene_name]}"
        SeqIO.write(upstream_seq, snakemake.output["fasta"], "fasta")
