import yaml

# ############ SETUP ##############################
configfile: "configs/config_coop2.yaml"
# configfile: "configs/config.json"
workdir: config['workdir']
snakedir = os.path.dirname(workflow.snakefile)

# include helper functions
include: "includes/io.snk"
include: "includes/utils.snk"



# retrieve the file_df with all the file paths from the samplesheet
sample_df, short_df, mut_df = get_files(config['inputdirs'], config['samples']['samplesheet'])
chrom_list = get_chrom_list(config)
# print(sample_df)
# print(mut_df)
# ############## INCLUDES ################################################
include: "includes/fastq.snk"
include: "includes/map.snk"
include: "includes/processBAM.snk"
include: "includes/dedup.snk"
include: "includes/umi_filter.snk"
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


# ############## MASTER RULE ##############################################

rule all:
    input:
        "QC/fastQC.html",
        expand("IGVnav/{sample}.txt", sample=sample_df.index),
        expand("IGVnav/{sample}.offTarget.txt", sample=sample_df.index),
        "results.txt"

###########################################################################
# print out of installed tools
onstart:
    print("    EXOM SEQUENCING PIPELINE STARTING.......")

    print('fastq:', short_df.loc[:, ['R1', 'R2']])
    ##########################
    path_to_config = os.path.join(config['workdir'], "config.yaml")
    with open(path_to_config, 'w+') as stream:
        yaml.dump(config, stream, default_flow_style=False)
    # create logs folder


onsuccess:
    # shell("export PATH=$ORG_PATH; unset ORG_PATH")
    print("Workflow finished - everything ran smoothly")
