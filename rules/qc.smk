import os

RD = config["result_dir"]

# Dynamically extract the base filename (removes .fq.gz or .fastq.gz)
# e.g., turns "/path/to/V350405307_L02_read_1.fq.gz" into "V350405307_L02_read_1"
R1_BASE = os.path.basename(config["lane_r1"]).replace(".fastq.gz", "").replace(".fq.gz", "")
R2_BASE = os.path.basename(config["lane_r2"]).replace(".fastq.gz", "").replace(".fq.gz", "")

rule fastqc_lane:
    input:
        r1 = config["lane_r1"],
        r2 = config["lane_r2"]
    output:
        html_r1 = RD + f"/qc/fastqc/{R1_BASE}_fastqc.html",
        zip_r1  = RD + f"/qc/fastqc/{R1_BASE}_fastqc.zip",
        html_r2 = RD + f"/qc/fastqc/{R2_BASE}_fastqc.html",
        zip_r2  = RD + f"/qc/fastqc/{R2_BASE}_fastqc.zip"
    params:
        outdir = RD + "/qc/fastqc",
        # Pass the base names into params so the shell block can use them
        r1_base = R1_BASE,
        r2_base = R2_BASE
    resources:
        runtime = 240
    log:
        RD + "/logs/fastqc/fastqc_lane.log"
    benchmark:
        RD + "/benchmarks/fastqc/fastqc_lane.tsv"
    threads:
        config["threads"]["fastqc"]
    container:
        config["containers"]["fastqc"]
    shell:
        "fastqc {input.r1} {input.r2} -o {params.outdir} --threads {threads} > {log} 2>&1"
        
rule multiqc:
    input:
        # Instead of hardcoding, tell Snakemake to wait for the exact files fastqc_lane produces
        fastqc_zip_r1 = rules.fastqc_lane.output.zip_r1,
        fastqc_zip_r2 = rules.fastqc_lane.output.zip_r2,
        star_log      = RD + "/starsolo/Log.final.out"
    output:
        RD + "/qc/multiqc_report.html"
    params:
        outdir = RD + "/qc",
        dirs   = RD + "/qc/fastqc " + RD + "/starsolo"
    resources:
        runtime = 30
    log:
        RD + "/logs/multiqc/multiqc.log"
    benchmark:
        RD + "/benchmarks/multiqc/multiqc.tsv"
    container:
        config["containers"]["multiqc"]
    shell:
        "multiqc {params.dirs} -o {params.outdir} --force > {log} 2>&1"