import os
import pandas as pd
import glob

configfile: "../../config/config.yaml"

# Define project output directories and file paths.
ProjDirPath = config["OutPath"] + "/" + config["ProjName"]
TmpDirPath = config["TmpDir"] + "/" + config["ProjName"]

# Get organism and genome_version from config (or use defaults if not provided).
organism = config.get("organism", "human")
genome_version = config.get("genome_version", "hg38")

def get_circ_ids(input_file):
    """Extract circular RNA IDs from info file"""
    circ_ids = []
    with open(input_file) as f:
        for line in f:
            fields = line.strip().split("\t")
            if len(fields) < 9:
                continue
            attrs = fields[8].strip()
            strand = fields[6].strip()
            # Skip entries with None strand
            if strand in [".", "?", "", "None"]:
                continue
            import re
            circ_id_match = re.search(r'circ_id\s+"([^"]+)"', attrs)
            if circ_id_match:
                circ_id = circ_id_match.group(1).strip()
                circ_ids.append(circ_id)
    return list(set(circ_ids))  # Return unique IDs

rule miRNA_predict_all:
    input:
        ProjDirPath + "/miRNA_binding_site_predictions/" + config["ProjName"] + ".circr_output.csv"

rule extract_circ_exons:
    input:
        info_file = ProjDirPath + "/Collapsed_results/Collapse_" + config["ProjName"] + ".info"
    output:
        expand(TmpDirPath + "/extract_circ_exons/{circ_id}.bed", 
               circ_id=get_circ_ids(ProjDirPath + "/Collapsed_results/Collapse_" + config["ProjName"] + ".info"))
    benchmark:
        ProjDirPath + "/benchmarks/extract_circ_exons/" + config["ProjName"] + ".tsv"
    log:
        ProjDirPath + "/logs/extract_circ_exons/" + config["ProjName"] + ".log"
    run:
        import re
        import os
        import time
        
        # Create the output directory
        out_dir = os.path.dirname(str(output[0]))
        os.makedirs(out_dir, exist_ok=True)
        
        # Log start
        with open(log[0], "w") as log_file:
            log_file.write(f"Starting extract_circ_exons at {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
            
        # Parse the input file once and organize by circ_id
        circ_data = {}
        skipped_count = 0
        with open(input.info_file) as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                fields = line.split("\t")
                if len(fields) < 9:
                    continue
                    
                chrom = fields[0].strip()
                strand = fields[6].strip()
                
                if strand in [".", "?", "", "None"]:
                    skipped_count += 1
                    continue
                    
                attrs = fields[8].strip()
                circ_id_match = re.search(r'circ_id\s+"([^"]+)"', attrs)
                if not circ_id_match:
                    continue
                    
                circ_id = circ_id_match.group(1).strip()
                
                if circ_id not in circ_data:
                    circ_data[circ_id] = {
                        "chrom": chrom,
                        "strand": strand,
                        "exons": []
                    }
                
                # Extract isoform exons
                isoform_matches = re.findall(r'isoform\s+"([^"]+)"', attrs)
                if isoform_matches:
                    isoform_str = ";".join(isoform_matches)
                    isoform_entries = re.split(r'[;,]', isoform_str)
                    for exon in isoform_entries:
                        exon = exon.strip()
                        if not exon:
                            continue
                        coords = exon.split('-')
                        if len(coords) == 2:
                            circ_data[circ_id]["exons"].append((coords[0].strip(), coords[1].strip()))
        
        # Write individual bed files directly to final location
        for circ_id, data in circ_data.items():
            output_file = os.path.join(out_dir, f"{circ_id}.bed")
            with open(output_file, "w") as f_out:
                for start, end in data["exons"]:
                    bed_line = f"{data['chrom']}\t{start}\t{end}\t{circ_id}\t.\t{data['strand']}"
                    f_out.write(bed_line + "\n")
                f_out.flush()
                os.fsync(f_out.fileno())
            
            # Verify file exists and has content
            if not os.path.exists(output_file) or os.path.getsize(output_file) == 0:
                raise Exception(f"Failed to create or write to {output_file}")
            
            # Add small delay to allow filesystem to sync
            time.sleep(0.1)
        
        # Log completion
        with open(log[0], "a") as log_file:
            log_file.write(f"Completed extract_circ_exons at {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
            log_file.write(f"Created {len(circ_data)} bed files\n")
            log_file.write(f"Skipped {skipped_count} entries with invalid strand information\n")

rule run_circr:
    input:
        bed_file = TmpDirPath + "/extract_circ_exons/{circ_id}.bed"
    output:
        csv_file = TmpDirPath + "/circr_results/{circ_id}.csv"
    params:
        organism = config["organism"],
        genome_version = config["genome_version"],
        annot_path = config["AnnotPath"],
        ref_path = config["RefPath"],
        rRNA_annot = config["rRNA_annot"],
        miRNA_annot = config["miRNA_annot"],
        miR_family_info = config["miR_family_info"],
        AGO_annot = config["AGO_annot"],
        validated_interactions = config["validated_interactions"],
        circbase_annot = config["circbase_annot"],
        tmp_dir = TmpDirPath + "/circr_tmp/{circ_id}"
    threads: config["ThreadNr"]
    conda:
        "../envs/circr_env.yaml"
    benchmark:
        ProjDirPath + "/benchmarks/circr/{circ_id}.tsv"
    log:
        ProjDirPath + "/logs/circr/{circ_id}.log"
    shell:
        """
        mkdir -p $(dirname {output.csv_file})
        mkdir -p {params.tmp_dir}
        python ../scripts/Circr.py \
               -i {input.bed_file} \
               -c \
               -s {params.organism} \
               -v {params.genome_version} \
               --gtf {params.annot_path} \
               --genome {params.ref_path} \
               --rRNA {params.rRNA_annot} \
               --miRNA {params.miRNA_annot} \
               --miR_family_info {params.miR_family_info} \
               --AGO {params.AGO_annot} \
               --validated_interactions {params.validated_interactions} \
               --circbase_annot {params.circbase_annot} \
               --threads {threads} \
               --tmp_dir {params.tmp_dir} \
               -o {output.csv_file} \
               --log {log} \
               --keep_tmp
        """

rule merge_circr_results:
    input:
        csv_files = expand(TmpDirPath + "/circr_results/{circ_id}.csv", 
                          circ_id=get_circ_ids(ProjDirPath + "/Collapsed_results/Collapse_" + config["ProjName"] + ".info"))
    output:
        merged_file = ProjDirPath + "/miRNA_binding_site_predictions/" + config["ProjName"] + ".circr_output.csv"
    benchmark:
        ProjDirPath + "/benchmarks/merge_circr_results/" + config["ProjName"] + ".tsv"
    log:
        ProjDirPath + "/logs/merge_circr_results/" + config["ProjName"] + ".log"
    run:
        import os
        import pandas as pd
        
        # Create output directory if it doesn't exist
        os.makedirs(os.path.dirname(output.merged_file), exist_ok=True)
        
        # Get list of valid input files
        valid_files = [f for f in input.csv_files if os.path.exists(f) and os.path.getsize(f) > 0]
        
        if not valid_files:
            # Create an empty file with header if no valid input files exist
            pd.DataFrame(columns=["Chrom", "Start", "End", "miRNA Name", "Circ Name", "Strand", 
                                "Seed Category", "ID", "Software Matched", "Validated", "AGO", "circBase ID"]).to_csv(
                output.merged_file, index=False)
            with open(log[0], "w") as f:
                f.write("No valid input files found. Created empty output file with header.\n")
        else:
            # Read and concatenate all valid files
            dfs = []
            for f in valid_files:
                try:
                    df = pd.read_csv(f)
                    if not df.empty:
                        dfs.append(df)
                except Exception as e:
                    with open(log[0], "a") as f:
                        f.write(f"Error reading file {f}: {str(e)}\n")
            
            if dfs:
                # Concatenate all dataframes
                merged_df = pd.concat(dfs, ignore_index=True)
                # Remove duplicates if any
                merged_df = merged_df.drop_duplicates()
                # Sort by chromosome and start position
                merged_df = merged_df.sort_values(["Chrom", "Start"])
                # Save to output file
                merged_df.to_csv(output.merged_file, index=False)
                
                with open(log[0], "w") as f:
                    f.write(f"Successfully merged {len(dfs)} files\n")
                    f.write(f"Total rows in output: {len(merged_df)}\n")
            else:
                # Create an empty file with header if no valid data found
                pd.DataFrame(columns=["Chrom", "Start", "End", "miRNA Name", "Circ Name", "Strand", 
                                    "Seed Category", "ID", "Software Matched", "Validated", "AGO", "circBase ID"]).to_csv(
                    output.merged_file, index=False)
                with open(log[0], "w") as f:
                    f.write("No valid data found in input files. Created empty output file with header.\n")