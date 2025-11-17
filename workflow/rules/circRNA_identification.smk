import pandas as pd

configfile: "../../config/config.yaml"  # Set config file.

ProjDirPath = config["OutPath"] + "/" + config["ProjName"] # Set project directory.
TmpDirPath = config["TmpDir"] + "/" + config["ProjName"]  # Set TEMPDIR.

# Read sample sheet.
Samplesheet = pd.read_csv(config["Samplesheet"],
                          delim_whitespace=True)

rule all_circRNA_identification:
        input:
            expand(ProjDirPath + "/circRNA_identification/{sample}/{sample}.cand_circ.fa",
                   zip,
                   sample=Samplesheet["sample_name"]),
            ProjDirPath + "/circRNA_identification/" + config["ProjName"] + ".resultssheet.csv",
            ProjDirPath + "/Collapsed_results/Collapse_" + config["ProjName"] + ".info"

rule call_circRNA:
    input:
        TmpDirPath + "/trimming/{sample}.fastq.gz"
    output:
        ProjDirPath + "/circRNA_identification/{sample}/{sample}.cand_circ.fa"
    params:
        outDir = ProjDirPath + "/circRNA_identification/{sample}/",
        refFa = config["RefPath"],
        refGTF = config["AnnotPath"],
        prefix = "{sample}"
    threads: config["ThreadNr"]
    conda: "../envs/circRNA_identification_env.yaml"
    benchmark:
        (ProjDirPath + "/benchmarks/call_circRNA/{sample}.tsv")
    log:
        (ProjDirPath + "/logs/call_circRNA/{sample}.log")
    shell:
        """
        CIRI-long call -i {input} \
                       -o {params.outDir} \
                       -r {params.refFa} \
                       -a {params.refGTF} \
                       -p {params.prefix} \
                       -t {threads}\
                       &> {log}
        """

rule create_results_sheet:
    input:
        expand(ProjDirPath + "/circRNA_identification/{sample}/{sample}.cand_circ.fa",
               zip,
               sample=Samplesheet["sample_name"])
    output:
        ProjDirPath + "/circRNA_identification/" + config["ProjName"] + ".resultssheet.csv"
    run:
        sampleL = Samplesheet["sample_name"]
        pathL = [ProjDirPath + "/circRNA_identification/" + sample + "/" + sample + ".cand_circ.fa" for sample in sampleL]
        resultsSheetDf = pd.DataFrame({"sample_name": sampleL,
                                       "results_path": pathL})
        resultsSheetDf.to_csv("".join(ProjDirPath + "/circRNA_identification/" + config["ProjName"] + ".resultssheet.csv"),
                              sep=" ",
                              header=False,
                              index=False)

rule colaps_isoforms:
    input:
        ProjDirPath + "/circRNA_identification/" + config["ProjName"] + ".resultssheet.csv"
    output:
        ProjDirPath + "/Collapsed_results/Collapse_" + config["ProjName"] + ".info"
    params:
        outDir = ProjDirPath + "/Collapsed_results",
        refFa = config["RefPath"],
        refGTF= config["AnnotPath"],
        prefix = config["ProjName"]
    threads: config["ThreadNr"]
    conda: "../envs/circRNA_identification_env.yaml"
    benchmark:
        (ProjDirPath + "/benchmarks/collapse_circRNA/" + config["ProjName"] + ".tsv")
    log:
        (ProjDirPath + "/logs/collapse_circRNA/" + config["ProjName"] + ".log")
    shell:
        """
        CIRI-long collapse -i {input} \
                           -o {params.outDir} \
                           -r {params.refFa} \
                           -a {params.refGTF} \
                           -p {params.prefix} \
                           -t {threads}\
                            &> {log}
        """