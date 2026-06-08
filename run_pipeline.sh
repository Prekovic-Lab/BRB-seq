#!/bin/bash

# 1. Define your custom temp directory
MY_TMP_DIR="/hpc/local/Rocky8/prekovic/theo/temporary/apptainer_cache"

# 2. Export variables so Apptainer/Singularity uses this space for building
export APPTAINER_TMPDIR="$MY_TMP_DIR"
export SINGULARITY_TMPDIR="$MY_TMP_DIR"
export TMPDIR="$MY_TMP_DIR"

# 3. (Optional) Set cache dir here if you don't use 'singularity-prefix' in the profile
export APPTAINER_CACHEDIR="$MY_TMP_DIR"
export SINGULARITY_CACHEDIR="$MY_TMP_DIR"

# 4. Run Snakemake
snakemake --profile my-slurm "$@"