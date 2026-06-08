configfile: "config/config.yaml"

import pandas as pd

_demux = pd.read_csv(config["demultiplex_tsv"], sep="\t")
_demux["sample_id"] = (
    _demux["Sample_name"]
    .str.strip()
    .str.replace(r"[\s/\\]", "_", regex=True)
)
SAMPLES = _demux["sample_id"].tolist()

RD = config["result_dir"]

# Order of include is important! In some cases the output of one rule is an input to another, so we need to define the rules in the correct order.
include: "rules/prepare_inputs.smk"
include: "rules/starsolo_new.smk"
include: "rules/qc_new.smk"
include: "rules/reformat_counts.smk"

rule all:
    input:
        RD + "/counts/all_samples_counts.tsv",
        RD + "/qc/multiqc_report.html",
        RD + "/qc/BRB_seq_QC_Plots.pdf"