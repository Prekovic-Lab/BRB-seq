## BRB-seq Processing Pipeline (STARsolo)
![Image](https://img.shields.io/badge/snakemake-≥7.0-brightgreen.svg)
![Image](https://img.shields.io/badge/apptainer-containerized-blue.svg)

A highly efficient, scalable, and fully containerized Snakemake pipeline for processing BRB-seq (Bulk RNA Barcoding and sequencing) data.
Instead of traditional, slow FASTQ demultiplexing, this pipeline leverages STARsolo to simultaneously map and demultiplex reads in a single step. It outputs a clean, ready-to-use Genes × Samples count matrix alongside comprehensive Quality Control reports.

## Key Features
- Lightning Fast: Uses STARsolo for combined mapping and demultiplexing. FastQC is run on a subset of reads (1M) to provide instant library quality metrics without bottlenecking the pipeline.
- Ready for Analysis: Automatically translates STARsolo Matrix Market (.mtx) outputs into a standard Pandas/R-friendly TSV count matrix with your actual sample names.
- Comprehensive QC: Generates custom R-based PDF reports (PCA, Depth vs. Genes Detected, counts per sample) and aggregates STAR/FastQC logs using MultiQC.
- Flexible BAM Handling: Optionally generate BAM files and split them into individual, sample-specific BAM files (generate_bam: True).
- Fully Containerized: Almost zero software installation required. All tools (STAR, MultiQC, Samtools, R, Python) are automatically pulled via Docker/Apptainer.
- Failsafe: Built-in validation instantly halts the pipeline if duplicate sample names are detected in your metadata.

## Prerequisites
To run this pipeline, you only need:
- A conda environment that contains:
  - Snakemake  
  - [snakemake-executor-plugin-slurm](https://anaconda.org/channels/bioconda/packages/snakemake-executor-plugin-slurm/overview)
  - Apptainer
- A pre-built STAR index for your target organism.

## Repository Structure
``` text
├── Snakefile                                     # Main Snakemake workflow definition
├── config/
│   ├── config.yaml                               # Main configuration file (paths, toggles)
│   ├── template_demultiplex_file.tsv             # Example metadata format
│   └── alitheia_brb_barcodes_v5D_384_brb.txt     # BRB barcode library
├── rules/                                        # Snakemake rule modules
│   ├── prepare_inputs.smk                        # Metadata validation and whitelist generation
│   ├── starsolo_new.smk                          # Mapping, demultiplexing, BAM splitting
│   ├── qc_new.smk                                # FastQC, custom R plots, MultiQC
│   └── reformat_counts.smk                       # Converts STARsolo output to final TSV matrix
├── scripts/
│   ├── plot_qc.R                                 # Generates PCA and read depth distributions
│   ├── prepare_starsolo_inputs.py
│   └── reformat_starsolo.py
└── run_pipeline.sh                               # Wrapper script for Slurm/HPC submission
```

## Setup & Configuration
1. Clone the repository:
2. Prepare your metadata (Demultiplexing file):  
Create a tab-separated file based on `config/template_demultiplex_file.tsv`.  
**Required columns:** `Sample_name`, `Barcode`  
**Optional column:** `Experiment` (in cases you are running multiple experiments at once and you would like different qc plots/counts).
3. Edit the config/config.yaml file:  
Adjust the file paths to point to your data. Key parameters include:
``` yaml
# ── Input Files ──
lane_r1: "/path/to/R1.fq.gz"
lane_r2: "/path/to/R2.fq.gz"
brb_barcodes: "path/to/the/barcodes/sequence/file"
demultiplex_tsv: "/path/to/your_barcodes_samples_identifiers.tsv"    # see step 2 right above
star_genome_dir: "/path/to/STAR_index"
result_dir: "path/to/results/directory"    # will be created if it doesn't exist

# ── Pipeline Toggles ──
UMIs_dedup: True    # True = Deduplicate UMIs (Standard for BRB-seq)
generate_bam: False # Set to True if you need BAM files (increases runtime/storage)
```

## Running the Pipeline
A wrapper script `run_pipeline.sh` is provided to handle Apptainer cache directories and submit the Snakemake job to a cluster.  
If you are using Slurm (and have a Snakemake slurm profile configured) run the following:
``` bash
bash run_pipeline.sh
```
To run it locally or interactively (e.g., allocating 16 cores):  
``` bash
snakemake --use-singularity --cores 16
```

## Output Files
The pipeline generates a highly organized results/ directory based on your config:  

| Directory/File | Description |
| :--- | :--- |
| **counts/all_samples_counts.tsv** | The final output. A TSV file where rows are Genes and columns are your specific Sample Names. |
| **qc/BRB_seq_QC_Plots.pdf** | Custom PDF report containing UMI distributions, Genes detected, and PCA plots (grouped by Experiment if provided). |
| **qc/multiqc_report.html** | Aggregated report of STAR mapping metrics and FastQC. |
| **starsolo/** | Raw STAR logs and the raw Matrix Market (.mtx) outputs. |
| **bams/split_by_sample/** | (Optional) Individual .bam and .bai files for every demultiplexed sample. |
| **logs/ & benchmarks/** | Resource usage and error logs for every executed rule. |

## How it Works
Input Validation (`prepare_inputs.smk`): Python instantly checks your TSV for duplicate sample names. It then cross-references your barcodes against the BRB dictionary to generate a STARsolo whitelist.  
Alignment & Demux (`starsolo_new.smk`): STARsolo processes the paired-end reads. By default, it uses a 1-mismatch directional UMI deduplication algorithm (`--soloUMIdedup 1MM_Directional`). If generate_bam is True, it will output a massive aligned BAM and seamlessly split it into sample-specific BAMs using samtools split.  
Reformatting (`reformat_counts.smk`): A custom Pandas script loads the sparse matrix, maps the barcode sequences to your Sample_names, drops empty droplets, and exports a dense TSV file.  
Quality Control (`qc_new.smk`): To avoid FastQC running for hours on massive merged lanes, the pipeline subsets the first 1 million reads. Meanwhile, an R script generates summary plots directly from your output count matrix.

## License
MIT License - Feel free to use, modify, and distribute this pipeline.

## Contributing
Issues and Pull Requests are welcome! If you encounter any issues regarding the pipeline please open an issue or contact me at t.chalkiadakis@umcutrecht.nl.
