# urchin_cisTarget
Workflow to generate cisTarget database for urchins

This is a snakemake workflow to generate a [cisTarget database](https://github.com/aertslab/create_cisTarget_databases) for the SCENIC workflows.
# Configuration

The workflow assumes access to a compiled `ClusterBuster` binary.

Modify the configuration (`config.yaml`) as appropriate for your data, with the following example entries:

```yaml
input:
  transfac_motifs: data/jaspar_urochordates_core_and_unvalidated.tf  # motifs in the transfact formats
  gff: /projectnb/bradham/data/ReferenceSequences/wray-genome/L_var_clean.gff  # gff file with genome annotations
  sequences: data/grn_genes.txt  # sequences of interest to include in the database (TFs, targets, etc.)
  genome_fasta: data/ReferenceSequences/wray-genome/Lvar_genome.fasta  # genome fasta
params:
  upstream_bp: 5000  # number of bases upstream from TSS to grab
  db_prefix: "lvar_grn_genes"  # output prefix
  clusterbuster: 'src/cbust'  # location of ClusterBluster
  cisTarget: create_cisTarget_databases/create_cistarget_motif_databases.py  # location of database creation script
```

**Notes:**
If you already have motifs in the `ClusterBuster` format, you can simply put them into the directory `output/motifs`. You will
then need to generate a motif id file for these motifs and save it as "output/motif_ids.txt"

# Running
Simply issue the command:

`snakemake --cores <n> --use-conda`
