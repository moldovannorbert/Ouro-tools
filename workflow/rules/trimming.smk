import pandas as pd

configfile: "../../config/config.yaml"  # Set config file.

ProjDirPath = config["OutPath"] + "/" + config["ProjName"] # Set project directory.
TmpDirPath = config["TmpDir"] + "/" + config["ProjName"]  # Set TEMPDIR.

# Read sample sheet.
Samplesheet = pd.read_csv(config["Samplesheet"],
                          delim_whitespace=True)

rule all_trim:
        input:
            expand(TmpDirPath + "/trimming/{porechop_sample}.fastq.gz",
                   zip,
                   porechop_sample=Samplesheet["sample_name"])

rule trim_porechop:
    input:
        config["SamplePath"] + "/{porechop_sample}.fastq.gz"
    output:
        TmpDirPath + "/trimming/{porechop_sample}.fastq.gz"
    params:
        end_bases_trimmed = config["ExtraEndTrim"]
    threads: config["ThreadNr"]
    conda: "../envs/trimming_env.yaml"
    benchmark:
        (ProjDirPath + "/benchmarks/trimming/{porechop_sample}.tsv")
    log:
        (ProjDirPath + "/logs/trimming/{porechop_sample}.log")
    shell:
        """
        porechop -i {input} \
                 -o {output} \
                 -t {threads} \
                 --extra_end_trim {params.end_bases_trimmed} \
                 &> {log}
        """