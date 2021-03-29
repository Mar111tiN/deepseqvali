import os
from os import system as run
from script_utils import show_output


def main(s):
    input = s.input
    output = s.output
    threads = s.threads

    extension = os.path.splitext(input[0])[1]
    if extension == ".fastq":
        # compress fastq as fastq.gz into workdir
        run(f"pigz -5 -p {threads} {input} > {output}")
    elif extension == ".gz":
        show_output(f"Creating symlink for file {input}")
        # create shortcut to fastq.gz in workdir/fastq
        if snakemake.config['setup']['demulti']:
            run(f"ln -s $(pwd)/{input} {output}")
        else:
            run(f"ln -s {input} {output}")
    elif extension == ".bz2":
        show_output(f"file extension {extension} --> unzipping with bzcat")
        # uncompress fastq.b2 and recompress to fastq.gz in workdir/fastq
        run(f"bzcat {input} | pigz -5 -p {threads} > {output}")


if __name__ == "__main__":
    main(snakemake)
