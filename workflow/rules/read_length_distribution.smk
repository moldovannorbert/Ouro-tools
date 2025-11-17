import pandas as pd

configfile: "../../config/config.yaml"  # Set config file.

ProjDirPath = config["OutPath"] + "/" + config["ProjName"] # Set project directory.
TmpDirPath = config["TmpDir"] + "/" + config["ProjName"]  # Set TEMPDIR.

# Read sample sheet.
Samplesheet = pd.read_csv(config["Samplesheet"],
                          delim_whitespace=True)

rule read_length_distribution_all:
    input:
        expand(ProjDirPath + "/read_length_distribution/{sample}.tsv",
               sample=Samplesheet["sample_name"])

rule read_length_distribution:
    input:
        ProjDirPath + "/mapping_QC/{sample}/NanoPlot-data.tsv.gz"
    output:
        ProjDirPath + "/read_length_distribution/{sample}.tsv"
    conda: "../envs/mapQC_env.yaml"
    benchmark:
        ProjDirPath + "/benchmarks/read_length_distribution/{sample}.tsv"
    log:
        ProjDirPath + "/logs/read_length_distribution/{sample}.log"
    shell:
        """
        python3 ../scripts/read_length_distribution.py \
        -i {input} \
        -o {output}
        """