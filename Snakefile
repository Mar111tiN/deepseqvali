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
global_dfs = get_files(config)
chrom_list = get_chrom_list(config)

print("sample", global_dfs['fastq'])
print("short", global_dfs['short_fastq'])

if config['setup']['demulti']:
    print("sample_index", global_dfs['samples'])
    if config['setup']['check_all_indices']:
        print("all_indices", global_dfs['index'])

# print(index_df)

reads=[1,2] if config['setup']['PE'] else [1]
# ############## INCLUDES ################################################
include: "includes/demulti.snk"
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
    samples_index = "[^_/.0-9]+"

# ############## MASTER RULE ##############################################

rule all:
    input:
        expand("fastq/{sample}.{read}.fastq.gz", sample=global_dfs['samples']['sample'], read=reads),
        expand("fastq/index/{sample}.{read}.fastq.gz", sample=global_dfs['index']['sample'], read=reads)
        # "QC/fastQC.html",
        # expand("IGVnav/{sample}.txt", sample=sample_df.index),
        # expand("IGVnav/{sample}.offTarget.txt", sample=sample_df.index),
        # "results.txt"

###########################################################################
# print out of installed tools
onstart:
    print("    EXOM SEQUENCING PIPELINE STARTING.......")

    print('fastq:', global_dfs['short_fastq'].loc[:, ['R1', 'R2']])
    ##########################
    path_to_config = os.path.join(config['workdir'], "config.yaml")
    with open(path_to_config, 'w+') as stream:
        yaml.dump(config, stream, default_flow_style=False)
    # create logs folder


onsuccess:
    # shell("export PATH=$ORG_PATH; unset ORG_PATH")
    print("Workflow finished - everything ran smoothly")
