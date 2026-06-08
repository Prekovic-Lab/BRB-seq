RD = config["result_dir"]

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
        bam      = RD + "/bams/Aligned.sortedByCoord.out.bam",
        star_log = RD + "/starsolo/Log.final.out"
    params:
        prefix   = RD + "/starsolo/",
        bam_dir  = RD + "/bams/",
        cb_len   = config["barcode_len"],
        umi_len  = config["umi_len"]
    log:
        RD + "/logs/starsolo/starsolo.log"
    benchmark:
        RD + "/benchmarks/starsolo/starsolo.tsv"
    resources:
        mem_mb = 50000,
        runtime = 700
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
            --soloUMIdedup 1MM_Directional \
            --soloBarcodeReadLength 0 \
            --soloCellFilter EmptyDrops_CR \
            --outSAMtype BAM SortedByCoordinate \
            --outSAMattributes NH HI AS NM CR UR CB UB GX GN sS sQ sM \
            --outFileNamePrefix {params.prefix} \
            --limitBAMsortRAM 67112659435 \
            --outBAMsortingThreadN {threads} \
        > {log} 2>&1

        mv {params.prefix}Aligned.sortedByCoord.out.bam {params.bam_dir}
        """

# For extremely large fq.gz files, STARsolo can either be terminated if not a lot of resources were asked and/or during the sorting of bams if not enough memory was requested.
# For big files add this --limitBAMsortRAM 67112659435 \ to the STAR command above and increase the mem_mb resource to 70000 or more (I did up to 250GB and requested 2 days for 450GB R1 & 370GB R2). 
# --soloUMIdedup 1MM_Directional \ is another option for UMI deduplication, but it is very slow and can be turned off with --soloUMIdedup NoDedup \


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