import yaml

# ############ SETUP ##############################
configfile: "configs/config_Lauf3.yaml"
# configfile: "configs/config.json"
workdir: config['workdir']
snakedir = os.path.dirname(workflow.snakefile)

# include helper functions
include: "includes/io.snk"
include: "includes/utils.snk"
include: "includes/demulti_utils.snk"


# retrieve the file_df with all the file paths from the samplesheet
# if demulti is true, short_df is the index_df (for picards barcodes.txt)
mut_df, sample_df = get_files(config)
chrom_list = get_chrom_list(config)

# print("sample", sample_df)
# print("mutations", mut_df)

reads=["R1", "R2"] if config['setup']['PE'] else ["R1"]
picard_reads = [1, 2] if config['setup']['PE'] else [1]
# ############## INCLUDES ################################################
include: "includes/demulti.snk"
include: "includes/fastq.snk"
include: "includes/map.snk"
include: "includes/processBAM.snk"
include: "includes/dedup.snk"
include: "includes/VAF.snk"
include: "includes/coverage.snk"

# specified wildcards have to match the regex
wildcard_constraints:
    # eg sample cannot contain _ or / to prevent ambiguous wildcards
    sample = "[^_/.]+",
    type = "[^_/.]+",
    read = "[^_/.]+",
    read_or_index = "[^_/.]+",
    chrom = "(chr)?[0-9XY]+",
    samples_index = "[^_/.0-9]+"

# ############## MASTER RULE ##############################################

rule all:
    input:
        # expand("fastq/{sample}_{read}.fastq.gz", sample=sample_df.index, read=reads),
        "QC/index-fastQC.html",
        "QC/samples-fastQC.html",
        # expand("IGVnav/{sample}.txt", sample=sample_df.index),
        # expand("IGVnav/{sample}.offTarget.txt", sample=sample_df.index),
        "results.txt"

###########################################################################
# print out of installed tools
onstart:
    print("    EXOM SEQUENCING PIPELINE STARTING.......")
    if config['setup']['PE']:
        print('fastq:', sample_df.loc[:, ['sR1', 'sR2']])
    else:
        print('fastq:', sample_df.loc[sample_df['samples_index'] == "samples", ['sR1']])
    ##########################
    path_to_config = os.path.join(config['workdir'], "config.yaml")
    with open(path_to_config, 'w+') as stream:
        yaml.dump(config, stream, default_flow_style=False)
    # create logs folder


onsuccess:
    # shell("export PATH=$ORG_PATH; unset ORG_PATH")
    print("Workflow finished - everything ran smoothly")

    if config['setup']['cleanup']:
        shell("rm -f picard/fastq/*barcode*fastq.gz")
        shell("rm -f fastqc/*_fastqc.html fastqc/*.sub")