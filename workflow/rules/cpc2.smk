import pandas as pd

configfile: "../../config/config.yaml"  # Set config file.

ProjDirPath = config["OutPath"] + "/" + config["ProjName"] # Set project directory.
TmpDirPath = config["TmpDir"] + "/" + config["ProjName"]  # Set TEMPDIR.


rule all_cpc2:
    input:
        ProjDirPath + "/coding_evaluation/cpc2_results.tsv"

rule extract_sequence_cpc2:
    input:
        circRNA_info = ProjDirPath + "/Collapsed_results/Collapse_" + config["ProjName"] + ".info",
        ref_genome = config["RefPath"]
    output:
        ProjDirPath + "/coding_evaluation/CircRNA_sequences.fasta"
    conda: 
        "../envs/cpc2_env.yaml"
    benchmark:
        ProjDirPath + "/benchmarks/extract_sequence_cpc2.tsv"
    log:
        ProjDirPath + "/logs/extract_sequence_cpc2.log"
    shell:
        """
        python ../scripts/extract_circRNA_sequences.py \
            --circRNA_info {input.circRNA_info} \
            --ref_genome {input.ref_genome} \
            --output {output} 2> {log}
        """

rule duplicate_sequences_cpc2:
    input:
        fasta = ProjDirPath + "/coding_evaluation/CircRNA_sequences.fasta",
        setup_complete = ProjDirPath + "/.cpc2_setup_complete"
    output:
        ProjDirPath + "/coding_evaluation/CircRNA_sequences_duplicated.fasta"
    conda:
        "../envs/cpc2_env.yaml"
    threads: config["ThreadNr"]
    benchmark:
        ProjDirPath + "/benchmarks/duplicate_sequences_cpc2.tsv"
    log:
        ProjDirPath + "/logs/duplicate_sequences_cpc2.log"
    shell:
        """
        python ../scripts/duplicate_sequences.py \
            -i {input.fasta} \
            -o {output} \
            -t {threads} \
            -v 2> {log}
        """

rule predict_coding_potential:
    input:
        ProjDirPath + "/coding_evaluation/CircRNA_sequences_duplicated.fasta"
    output:
        ProjDirPath + "/coding_evaluation/cpc2_results.tsv"
    conda:
        "../envs/cpc2_env.yaml"
    threads: config["ThreadNr"]
    benchmark:
        ProjDirPath + "/benchmarks/predict_coding_potential.tsv"
    log:
        ProjDirPath + "/logs/predict_coding_potential.log"
    shell:
        """
        # CPC2_output_peptide.py adds .txt extension automatically, so we pass path without extension
        output_base=$(echo {output} | sed 's/\\.tsv$//')
        python ../scripts/CPC2_output_peptide.py -i {input} -o $output_base --ORF 2> {log}
        
        # Rename .txt to .tsv
        txt_file="$output_base.txt"
        if [ -f "$txt_file" ]; then
            mv "$txt_file" {output}
        else
            echo "Error: Expected output file $txt_file not found" >&2
            exit 1
        fi
        
        # Ensure output file exists and is not empty
        if [ ! -s {output} ]; then
            echo "Error: Output file {output} is missing or empty" >&2
            exit 1
        fi
        """