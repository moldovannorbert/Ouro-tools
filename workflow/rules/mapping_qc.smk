import pandas as pd

configfile: "../../config/config.yaml"  # Set config file.

ProjDirPath = config["OutPath"] + "/" + config["ProjName"] # Set project directory.
TmpDirPath = config["TmpDir"] + "/" + config["ProjName"]  # Set TEMPDIR.

# Read sample sheet.
Samplesheet = pd.read_csv(config["Samplesheet"],
                          delim_whitespace=True)

rule mappingQC_all:
    input:
        directory(expand(ProjDirPath + "/mapping_QC/{sample}",
                  sample=Samplesheet["sample_name"]))

rule map:
    input:
        TmpDirPath + "/trimming/{sample}.fastq.gz"
    output:
        TmpDirPath + "/mapping/{sample}.bam"
    params:
        ref = config["RefPath"]
    threads: config["ThreadNr"]
    conda: "../envs/mapQC_env.yaml"
    benchmark:
        ProjDirPath + "/benchmarks/mapping/{sample}.tsv"
    log:
        ProjDirPath + "/logs/mapping/{sample}.log"
    shell:
        """
        minimap2 -ax splice {params.ref} \
                 {input} \
                 -t {threads} \
                 2> {log} |
                 samtools sort \
                 -@ {threads} \
                 -o {output} 2>> {log}
        """

rule filter_and_index_mappingQC:
    input:
        TmpDirPath + "/mapping/{sample}.bam"
    output:
        bam = TmpDirPath + "/filter_index/{sample}.bam",
        bai = TmpDirPath + "/filter_index/{sample}.bam.bai"
    threads: config["ThreadNr"]
    params:
        qual = 5
    conda: "../envs/mapQC_env.yaml"
    benchmark:
        ProjDirPath + "/benchmarks/filter_index/{sample}.tsv"
    shell:
        """
        samtools view {input} \
        -h `# include header` \
        -q {params.qual} `#> only reads with high quality` \
        -F 4 `# not unmapped` \
        -F 256 `# no secondary alignment; thus primary alignment` \
        -F 1024 `# no PCR duplicate` \
        -F 2048 `# no supplementary (=chimeric) alignment` \
        -@ {threads} \
        -b | tee {output.bam} | \
        \
        samtools index - `# stdin` \
        -@ {threads} \
        {output.bai}
        """

rule NanoPlot:
        input:
            TmpDirPath + "/filter_index/{sample}.bam"
        output:
            directory(ProjDirPath + "/mapping_QC/{sample}")
        threads: config["ThreadNr"]
        conda: "../envs/mapQC_env.yaml"
        benchmark:
            ProjDirPath + "/benchmarks/mapping_QC/{sample}.tsv"
        log:
            ProjDirPath + "/logs/mapping_QC/{sample}.log"
        shell:
            """
            NanoPlot --bam {input} \
                     -o {output} \
                     --raw \
                     --alength \
                     -t {threads} \
                     --huge \
                     2>> {log}
            """