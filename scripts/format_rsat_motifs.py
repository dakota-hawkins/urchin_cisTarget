import re

if __name__ == '__main__':
    try:
        snakemake
    except NameError:
        snakemake = None
    if snakemake is not None:
        tf_file = snakemake.input['rsat']
        out_file = snakemake.output['tfs']
        motifs = []
        desc_regex = re.compile('^DE[\s]+')
        consensus_regex = re.compile("consensus.strict\:[\s]+")
        matrix_regex = re.compile("^[0-9]+[\s]+")
        with open(tf_file, 'r') as f:
            line = f.readline()
            while line != '':
                while desc_regex.search(line) is None and line != '':
                    line = f.readline()
                desc = desc_regex.split(line)[-1]
                for i in range(2):
                    line = f.readline()
                matrix = ''
                while line[:2] != 'XX' and line != '':
                    matrix += matrix_regex.split(line)[1]
                    line = f.readline()
                line = f.readline()
                # read description lines, check for consensus description. 
                while consensus_regex.search(line) is None and line[:2] != 'XX' and line != '':
                    line = f.readline()
                if consensus_regex.search(line) is not None:
                    desc = consensus_regex.split(line)[-1]
                motifs.append(">" + desc + matrix)
                line = f.readline()
        with open(out_file, 'w') as f:
            f.write("".join(motifs))