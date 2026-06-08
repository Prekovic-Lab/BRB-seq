# scripts/plot_qc.R
# Suppress package startup messages to keep logs clean
suppressPackageStartupMessages({
    library(ggplot2)
})

# ── Logging ──────────────────────────────────────────────────────────────────
log_file <- snakemake@log[[1]]
con <- file(log_file, open = "wt")
sink(con, type = "output")
sink(con, type = "message")

# ── I/O ──────────────────────────────────────────────────────────────────────
counts_file <- snakemake@input[["counts"]]
output_pdf  <- snakemake@output[["report"]]

# Ensure output directory exists
dir.create(dirname(output_pdf), recursive = TRUE, showWarnings = FALSE)

# 1. Capture inputs and outputs from Snakemake
counts_file <- snakemake@input[["counts"]]
output_pdf  <- snakemake@output[["report"]]

# 2. Read the counts matrix
# Assumes standard format: rows = genes, cols = samples (with a header)
# check.names=FALSE prevents R from converting hyphens in sample names to dots
counts_mat <- read.table(counts_file, header=TRUE, row.names=1, sep="\t", check.names=FALSE)

# 3. Calculate QC metrics per sample
# Total counts (UMIs) per sample
total_umis <- colSums(counts_mat)
# Total unique genes detected (count > 0) per sample
genes_detected <- colSums(counts_mat > 0)

# Build a dataframe for ggplot
qc_df <- data.frame(
    Sample = names(total_umis),
    Total_UMIs = total_umis,
    Genes_Detected = genes_detected
)

# Order the dataframe by Total UMIs (Highest to Lowest) for a "Knee Plot" style drop-off
qc_df$Sample <- factor(qc_df$Sample, levels = qc_df$Sample[order(qc_df$Total_UMIs, decreasing = TRUE)])

# 4. Generate the PDF
pdf(output_pdf, width = 10, height = 7)

# --- PLOT 1: Histogram of Total Counts ---
p1 <- ggplot(qc_df, aes(x = Total_UMIs)) +
    geom_histogram(bins = 30, fill = "steelblue", color = "black", alpha = 0.8) +
    theme_minimal() +
    labs(title = "1. Distribution of Total Counts per Sample",
         subtitle = "Quickly see if most of your wells yielded a decent number of reads",
         x = "Total UMIs per Sample",
         y = "Frequency (Number of Samples)")

# --- PLOT 2: Bar Plot of Total Counts per Sample ---
p2 <- ggplot(qc_df, aes(x = Sample, y = Total_UMIs)) +
    geom_bar(stat = "identity", fill = "coral") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6)) +
    labs(title = "2. Total UMIs per Sample",
         subtitle = "Ordered highest to lowest. Look for massive drops or hogs.",
         x = "Sample ID",
         y = "Total UMIs")

# --- PLOT 3: Bar Plot of Genes Detected per Sample ---
p3 <- ggplot(qc_df, aes(x = Sample, y = Genes_Detected)) +
    geom_point(color = "seagreen", size = 2) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6)) +
    labs(title = "3. Unique Genes Detected per Sample",
         subtitle = "How many distinct genes had > 0 counts?",
         x = "Sample ID",
         y = "Genes Detected")

# --- PLOT 4: Saturation Check (Depth vs Genes) ---
p4 <- ggplot(qc_df, aes(x = Total_UMIs, y = Genes_Detected)) +
    geom_point(color = "purple", alpha = 0.7, size = 3) +
    theme_minimal() +
    labs(title = "4. Sequencing Depth vs Genes Detected",
         subtitle = "If this curve is completely flat, sequencing is highly saturated.",
         x = "Total UMIs",
         y = "Genes Detected")

# Print plots to the active PDF device
print(p1)
print(p2)
print(p3)
print(p4)

# Close the device to save the file
dev.off()

# ── Close log ────────────────────────────────────────────────────────────────
sink(type = "message")
sink()
close(con)