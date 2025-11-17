import pandas as pd

configfile: "../../config/config.yaml"  # Set config file.

ProjDirPath = config["OutPath"] + "/" + config["ProjName"] # Set project directory.
TmpDirPath = config["TmpDir"] + "/" + config["ProjName"]  # Set TEMPDIR.

# Read sample sheet.
Samplesheet = pd.read_csv(config["Samplesheet"],
                          delim_whitespace=True)

rule all_dedup:
    input:
        expand(TmpDirPath + "/deduplication/{sample}.fastq.gz",
               zip,
               sample=Samplesheet["sample_name"]),
        expand(ProjDirPath + "/deduplication_stats/{sample}_dedup_stats.tsv",
               zip,
               sample=Samplesheet["sample_name"])

rule deduplicate_reads:
    input:
        TmpDirPath + "/trimming/{sample}.fastq.gz"
    output:
        reads=TmpDirPath + "/deduplication/{sample}.fastq.gz",
        stats=ProjDirPath + "/deduplication_stats/{sample}_dedup_stats.tsv"
    params:
        similarity=config["DedupSimilarity"],
        end_length=config["DedupEndLength"]
    threads: config["ThreadNr"]
    conda: "../envs/deduplication_env.yaml"
    log:
        ProjDirPath + "/logs/deduplication/{sample}.log"
    benchmark:
        ProjDirPath + "/benchmarks/deduplication/{sample}.tsv"
    shell:
        """
        python ../scripts/deduplicate_reads.py \
            -i {input} \
            -o {output.reads} \
            --stats {output.stats} \
            --similarity {params.similarity} \
            --end-length {params.end_length} \
            -t {threads} \
            &> {log}
        """
