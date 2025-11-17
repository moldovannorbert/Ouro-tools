# Ouro-tools: Pipeline for the Analysis of circRNA from Long-Read Nanopore Sequencing

## Overview
Ouro-tools is an end-to-end, Snakemake-driven workflow for the discovery, quantification and functional annotation of circular RNAs (circRNAs) from Oxford Nanopore long-read sequencing data.  Starting from raw FASTQ files, the pipeline performs read trimming, alignment, quality control, circRNA calling, isoform collapsing, and downstream analyses such as miRNA binding-site discovery and IRES prediction.  All steps are automatically executed inside reproducible Conda environments to guarantee portability and full provenance tracking.

Key features:
* Automatic orchestration with Snakemake ≥8
* Reproducible Conda environments for every rule
* Support for multiplexed projects via a simple sample sheet
* Modular rules that can be re-used or skipped as required
* Rich HTML and text reports for QC, mapping statistics and functional annotation

## Repo Contents
```
.
├── config/                # Default project configuration (YAML)
├── workflow/
│   ├── envs/             # Conda recipes (one per rule)
│   ├── rules/            # Snakemake rules – main workflow logic
│   └── scripts/          # Custom Python/R helper scripts
├── R/                     # R markdown notebooks for statistics & visualisation
├── test/                  # Minimal dataset and utility scripts for CI/testing
├── LICENSE                # MIT license
└── README.md              # You are here
```

## System Requirements
1. Operating system: Linux or macOS (tested on Ubuntu ≥20.04)
2. Conda ≥4.10 (Miniconda or Mamba strongly recommended)
3. Snakemake ≥8 installed in the base environment
4. Python ≥3.8 and R ≥4.0 (handled automatically through Conda)
5. ~20 GB disk space for reference indexes & intermediate files
6. ≥16 GB RAM and 8 CPU threads recommended for speed; all thread counts are user-configurable.

External tools pulled into the workflow include, but are not limited to: minimap2, samtools, CIRI-long, NanoPlot, CPC2, DeepIRES, and Circr (which integrates RNAhybrid, miRanda, and TargetScan for miRNA binding site prediction).  They are installed on-the-fly via the `workflow/envs/*.yaml` environments.

**Note:** Mamba is a prerequisite for this workflow and can be installed following the [Mamba installation guide](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html).

## Installation Guide
```bash
# 1. Clone the repository
$ git clone https://github.com/moldovannorbert/Ouro-tools.git
$ cd Ouro-tools

# 2. Create and activate a clean conda environment that only holds snakemake
$ mamba create -n ouro-tools snakemake>=8 python=3.10 -y
$ mamba activate ouro-tools
```
All rule-specific dependencies will be created the first time the rule is executed and cached under `$HOME/.snakemake/conda`.

## User Guide
1. **Prepare the configuration**
   * Edit `config/config.yaml` to point to your project output directory, temporary directory, reference genome FASTA, GTF annotation and other run-time parameters.
   * Set the `SamplePath` key to the directory containing your FASTQ files. Files should be named as `{sample_name}.fastq.gz` where `{sample_name}` matches the sample names in your sample sheet.
2. **Create a sample sheet** with a `sample_name` column:
   ```
   sample_name
   SAMPLE1
   SAMPLE2
   ```
   The workflow will locate FASTQ files using `SamplePath` + `{sample_name}.fastq.gz`. Set the path to your sample sheet file via the `Samplesheet` key in `config.yaml`.
3. **Acquire external reference files**  
   The workflow relies on several reference resources that are **not** shipped with the repository. Download or prepare the following and supply their absolute paths in `config/config.yaml`:

   * **Reference genome (FASTA)** – The workflow uses the FASTA file directly with minimap2; no pre-built index is required. Set the `RefPath` field to the absolute path of your reference genome FASTA file.

   * **Reference gene annotation (GTF)** – Set the `AnnotPath` field.

   * **miRNA binding-site prediction bundle** – pre-compiled resources containing mature miRNA sequences, rRNA/AGO tracks, validated interactions, etc. The miRNA prediction module uses Circr.py which integrates RNAhybrid, miRanda, and TargetScan for comprehensive miRNA binding site prediction. Download the [archive](https://drive.google.com/drive/folders/1zJVyzEFAMtvZTTueWRocxXs63jUxsl-U) and unpack it. Then fill in the corresponding keys in `config.yaml` ( `rRNA_annot`, `miRNA_annot`, `miR_family_info`, `AGO_annot`, `validated_interactions`, `circbase_annot` ).

   Keeping the files on a shared/high-speed filesystem is recommended if you run the pipeline on a cluster.
4. **Run the workflow modules sequentially**  
   Ouro-tools separates major processing stages into individual Snakemake files so that you can checkpoint and inspect results after every step.  Execute them **in the following order**, each time pointing `--snakefile` to the corresponding file under `workflow/rules/`:

   | Order | Module (.smk)                   | Short description |
   |-------|--------------------------------|-------------------|
   | 1     | `trimming.smk`                 | Adapter removal and optional end-trimming of ONT reads. Outputs cleaned FASTQ files in `tmp/trimming/`. |
   | 2     | `mapping_qc.smk`               | Align trimmed reads to the reference genome with minimap2, filter low-quality alignments and generate NanoPlot QC reports. |
   | 3     | `read_length_distribution.smk` | Compute and plot read-length distributions to assess library quality. |
   | 4     | `circRNA_identification.smk`   | Detect circRNAs with CIRI-long and collapse isoforms across samples. |
   | 5     | `cpc2.smk`                     | Evaluate coding potential of circRNAs using CPC2. Outputs coding potential scores and ORF predictions. |
   | 6     | `miRNA_predict.smk`            | Predict miRNA binding sites on called circRNAs using Circr (RNAhybrid, miRanda, and TargetScan) with validated interaction databases. |
   | 7     | `IRES_predict.smk`             | Scan circRNA sequences for internal ribosome entry sites (IRES) with DeepIRES and IRES-like motif heuristics. Merges results with CPC2 ORF predictions. |

   Example command pattern (bash loop):
   ```bash
   for module in trimming mapping_qc read_length_distribution circRNA_identification cpc2 miRNA_predict IRES_predict; do
       snakemake --snakefile workflow/rules/${module}.smk \
                 -j 16 --use-conda
   done
   ```
   Replace `-j 16` with the number of cores you wish to allocate; add `--profile <cluster>` if submitting to an HPC scheduler.

5. **Inspect results** – by default all results are written under:
   ```
   <OutPath>/<ProjName>/
   ├── circRNA_identification/     # per-sample CIRI-long results
   ├── Collapsed_results/          # merged circRNA isoforms
   ├── coding_evaluation/          # CPC2 coding potential scores and ORF predictions
   ├── miRNA_binding_site_predictions/  # Circr output with miRNA binding sites
   ├── ires_prediction/            # IRES predictions (DeepIRES, IRES-like motifs, merged with ORF data)
   ├── mapping_QC/                 # NanoPlot QC reports (HTML)
   ├── read_length_distribution/   # read length distribution statistics
   └── benchmarks/                 # run-time & resource usage stats
   ```

6. **Visualisation & downstream analyses**
   The R notebooks under the `R/` directory can be rendered after the workflow has completed to produce summary statistics and publication-ready figures.

```

### Common flags
* `--printshellcmds` – show the wrapped shell commands
* `--keep-going` – continue independent jobs after an error
* `--rerun-incomplete` – re-run tasks with incomplete output

## License
Ouro-tools is released under the terms of the MIT license (see `LICENSE`).

## Citation
If you use Ouro-tools in your research, please cite:
> Wever B., *et al.* **Ouro-seq: Improved Recovery of Full-Length circRNAs from Samples with Limited RNA Content Using Short-Amplicon Suppression**. 2025. *In preparation*.

Please also cite the third-party tools wrapped by this workflow, notably Snakemake, CIRI-long, minimap2, samtools, NanoPlot, CPC2, Circr (RNAhybrid, miRanda, TargetScan), and DeepIRES. 