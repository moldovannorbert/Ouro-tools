import pandas as pd

configfile: "../../config/config.yaml"  # Set config file.

ProjDirPath = config["OutPath"] + "/" + config["ProjName"] # Set project directory.
TmpDirPath = config["TmpDir"] + "/" + config["ProjName"]  # Set TEMPDIR.


rule all_ires_prediction:
    input: 
        ProjDirPath + "/ires_prediction/ires_prediction_results.csv",
        ProjDirPath + "/ires_prediction/ireslike_prediction_results.csv",
        ProjDirPath + "/ires_prediction/ires_prediction_results_merged.csv"

rule extract_sequence:
    input:
        circRNA_info = ProjDirPath + "/Collapsed_results/Collapse_" + config["ProjName"] + ".info",
        ref_genome = config["RefPath"]
    output:
        ProjDirPath + "/ires_prediction/CircRNA_sequences.fasta"
    conda: 
        "../envs/ires_prediction_env.yaml"
    benchmark:
        ProjDirPath + "/benchmarks/extract_sequence.tsv"
    log:
        ProjDirPath + "/logs/extract_sequence.log"
    shell:
        """
        python ../scripts/extract_circRNA_sequences.py \
            --circRNA_info {input.circRNA_info} \
            --ref_genome {input.ref_genome} \
            --output {output} 2> {log}
        """

rule predict_ires:
    input:
        ProjDirPath + "/ires_prediction/CircRNA_sequences.fasta"
    output:
        ProjDirPath + "/ires_prediction/ires_prediction_results.csv"
    conda: 
        "../envs/ires_prediction_env.yaml"
    benchmark: 
        ProjDirPath + "/benchmarks/predict_ires.tsv"
    log:
        ProjDirPath + "/logs/predict_ires.log"
    shell:
        """
        python ../scripts/DeepIRES/DeepIRES.py \
            -i {input} \
            -o {output} 2> {log}
        """

rule predict_ireslike:
    input:
        fasta = ProjDirPath + "/ires_prediction/CircRNA_sequences.fasta"
    output:
        results = ProjDirPath + "/ires_prediction/ireslike_prediction_results.csv",
        summary = ProjDirPath + "/ires_prediction/ireslike_prediction_results_summary.csv"
    params:
        hexamers = config["IRES-like_motifs"]
    conda:
        "../envs/ires_prediction_env.yaml"  
    threads: 8
    benchmark:
        ProjDirPath + "/benchmarks/predict_ireslike.tsv"
    log:
        ProjDirPath + "/logs/predict_ireslike.log"
    shell:
        """
        python ../scripts/find_ires_like.py \
            -i {input.fasta} \
            -o {output.results} \
            -x "{params.hexamers}" \
            -t {threads} \
            -v 2> {log}
        """

rule merge_ires_prediction:
    input:
        ires_prediction = ProjDirPath + "/ires_prediction/ires_prediction_results.csv",
        ireslike_prediction = ProjDirPath + "/ires_prediction/ireslike_prediction_results.csv",
        circRNA_sequences = ProjDirPath + "/coding_evaluation/CircRNA_sequences_duplicated.fasta",
        cpc2_results = ProjDirPath + "/coding_evaluation/cpc2_results.tsv"
    output:
        merged_results = ProjDirPath + "/ires_prediction/ires_prediction_results_merged.csv",
        orf_summary = ProjDirPath + "/ires_prediction/ires_orf_summary.csv"
    conda:
        "../envs/ires_prediction_env.yaml"
    threads: 8
    benchmark:
        ProjDirPath + "/benchmarks/merge_ires_prediction.tsv"
    log:
        ProjDirPath + "/logs/merge_ires_prediction.log"
    shell:
        """
        python ../scripts/merge_ires_predictions.py \
            --ires {input.ires_prediction} \
            --ireslike {input.ireslike_prediction} \
            --fasta {input.circRNA_sequences} \
            --cpc2 {input.cpc2_results} \
            --output {output.merged_results} \
            --summary {output.orf_summary} \
            --threads {threads} \
            -v 2> {log}
        """
    