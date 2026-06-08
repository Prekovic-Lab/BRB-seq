RD = config["result_dir"]

# Check the config file for the BAM setting (default to True if missing)
GENERATE_BAM = config.get("generate_bam", True)

rule starsolo:
    input:
        r1        = config["lane_r1"],
        r2        = config["lane_r2"],
        whitelist = RD + "/starsolo/whitelist.txt",
        genome    = config["star_genome_dir"],
        gtf       = config["gtf"]
    output:
        matrix   = RD + "/starsolo/Solo.out/Gene/raw/matrix.mtx",
        barcodes = RD + "/starsolo/Solo.out/Gene/raw/barcodes.tsv",
        features = RD + "/starsolo/Solo.out/Gene/raw/features.tsv",
        star_log = RD + "/starsolo/Log.final.out",
        star_summary = RD + "/starsolo/Solo.out/Gene/Summary.csv",
        # Conditionally require the BAM file as an output
        bam      = [RD + "/bams/Aligned.sortedByCoord.out.bam"] if GENERATE_BAM else []
    params:
        prefix   = RD + "/starsolo/",
        bam_dir  = RD + "/bams/",
        cb_len   = config["barcode_len"],
        umi_len  = config["umi_len"],
        
        # Convert Python boolean to a string we can use in bash
        make_bam = "True" if GENERATE_BAM else "False",
        
        # Dynamically build the BAM arguments for the STAR command
        bam_args = lambda wildcards, threads: (
            f"--outSAMtype BAM SortedByCoordinate "
            f"--outSAMattributes NH HI AS NM CR UR CB UB GX GN sS sQ sM "
            f"--limitBAMsortRAM 67112659435 "
            f"--outBAMsortingThreadN {threads}"
        ) if GENERATE_BAM else "--outSAMtype None"
        
    log:
        RD + "/logs/starsolo/starsolo.log"
    benchmark:
        RD + "/benchmarks/starsolo/starsolo.tsv"
    resources:
        mem_mb = 250000,
        runtime = 2880
    threads:
        config["threads"]["star"]
    container:
        config["containers"]["star"]
    shell:
        """
        umi_start=$(( {params.cb_len} + 1 ))

        STAR \
            --runMode alignReads \
            --runThreadN {threads} \
            --genomeDir {input.genome} \
            --sjdbGTFfile {input.gtf} \
            --readFilesIn {input.r2} {input.r1} \
            --readFilesCommand zcat \
            --soloType CB_UMI_Simple \
            --soloCBwhitelist {input.whitelist} \
            --soloCBstart 1 --soloCBlen {params.cb_len} \
            --soloUMIstart $umi_start --soloUMIlen {params.umi_len} \
            --clipAdapterType CellRanger4 \
            --soloFeatures Gene \
            --soloUMIdedup NoDedup \
            --soloBarcodeReadLength 0 \
            --soloCellFilter EmptyDrops_CR \
            --outFileNamePrefix {params.prefix} \
            {params.bam_args} \
        > {log} 2>&1

        # Only move the BAM file if we actually told STAR to create it
        if [ "{params.make_bam}" == "True" ]; then
            mv {params.prefix}Aligned.sortedByCoord.out.bam {params.bam_dir}
        fi
        """

rule samtools_index:
    input:
        RD + "/bams/Aligned.sortedByCoord.out.bam"
    output:
        RD + "/bams/Aligned.sortedByCoord.out.bam.bai"
    log:
        RD + "/logs/samtools_index/samtools_index.log"
    benchmark:
        RD + "/benchmarks/samtools_index/samtools_index.tsv"
    container:
        config["containers"]["samtools"]
    shell:
        "samtools index {input} > {log} 2>&1"