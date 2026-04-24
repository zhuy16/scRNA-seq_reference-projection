#!/usr/bin/env bash
# =============================================================================
# run_loo_benchmark.sh — Leave-One-Out Patient Benchmark
# =============================================================================
# For each patient in the Yost 2019 BCC dataset:
#   1. Run CCA pipeline with that patient held out as query
#   2. Run scVI pipeline on the same split
#   3. Extract true/predicted labels into data/loo/fold_<patient>_cca.csv
#      and data/loo/fold_<patient>_scvi.csv
#
# After all folds complete, run benchmarking/benchmark_celltype.ipynb
# (§4 LOO section) to aggregate results across patients.
#
# Usage:
#   conda run -n scrnaseq bash run_loo_benchmark.sh
#   bash run_loo_benchmark.sh [patient1 patient2 ...]   # subset of patients
#
# Estimated runtime: ~11 min × 11 patients × 2 methods ≈ 4 h
#   (CCA and scVI run sequentially per fold; parallelise manually if needed)
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# ── Colours ───────────────────────────────────────────────────────────────────
GREEN='\033[0;32m'; YELLOW='\033[1;33m'; RED='\033[0;31m'; CYAN='\033[0;36m'; NC='\033[0m'
log()  { echo -e "${GREEN}[$(date '+%H:%M:%S')] $*${NC}"; }
warn() { echo -e "${YELLOW}[$(date '+%H:%M:%S')] WARN: $*${NC}"; }
fail() { echo -e "${RED}[$(date '+%H:%M:%S')] FAILED: $*${NC}"; exit 1; }

# ── macOS ARM OpenMP fix ───────────────────────────────────────────────────────
export DYLD_INSERT_LIBRARIES=/opt/anaconda3/envs/scrnaseq/lib/libomp.dylib
export KMP_DUPLICATE_LIB_OK=TRUE
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

# ── Patient list ──────────────────────────────────────────────────────────────
# All 11 Yost 2019 BCC patients (su001–su012, no su011)
ALL_PATIENTS=(su001 su002 su003 su004 su005 su006 su007 su008 su009 su010 su012)

# Allow overriding patients via command-line arguments
if [[ $# -gt 0 ]]; then
    PATIENTS=("$@")
else
    PATIENTS=("${ALL_PATIENTS[@]}")
fi

LOO_DIR="${SCRIPT_DIR}/data/loo"
mkdir -p "$LOO_DIR"
LOG_FILE="${LOO_DIR}/loo_benchmark.log"

log "LOO benchmark — patients: ${PATIENTS[*]}"
log "Outputs: ${LOO_DIR}/"
log "Full log: ${LOG_FILE}"

# ── Helper: read one value from params.yaml ───────────────────────────────────
yaml_val() {
    python3 -c "import yaml; d=yaml.safe_load(open('config/params.yaml')); print(d.get('$1',''))"
}

# ── Read shared params ────────────────────────────────────────────────────────
RAW_COUNTS_PATH=$(yaml_val raw_counts_path)
RAW_METADATA_PATH=$(yaml_val raw_metadata_path)
RAW_TCR_PATH=$(yaml_val raw_tcr_path)
CAUSHI_CD8_GZ=$(yaml_val caushi_cd8_gz)
CAUSHI_CLONOTYPES_PATH=$(yaml_val caushi_clonotypes_path)
N_DIMS=$(yaml_val n_dims)
MIN_GENES=$(yaml_val min_genes)
MAX_GENES=$(yaml_val max_genes)
MAX_PCT_MITO=$(yaml_val max_pct_mito)
MIN_CELLS_PER_CLONOTYPE=$(yaml_val min_cells_per_clonotype)
N_EXHAUSTED_CLONOTYPES=$(yaml_val n_exhausted_clonotypes)

NB_DIR="${SCRIPT_DIR}/notebooks"
OUT_DIR="${NB_DIR}/executed"
mkdir -p "$OUT_DIR"

# ── run_notebook helper (same pattern as run_cca_pipeline.sh) ─────────────────
run_notebook() {
    local nb_name="$1" kernel="$2"; shift 2
    local extra_params=("$@")
    local input="${NB_DIR}/${nb_name}"
    local output="${OUT_DIR}/${FOLD_STAMP}_${nb_name}"
    mkdir -p "$(dirname "$output")"
    [ -f "$input" ] || fail "Notebook not found: $input"
    local param_str=""
    for p in "${extra_params[@]}"; do param_str="${param_str} -p ${p}"; done
    # shellcheck disable=SC2086
    papermill "$input" "$output" --kernel "${kernel}" --report-mode $param_str 2>&1
}

# ── Extract labels from CCA RDS ───────────────────────────────────────────────
extract_cca_labels() {
    local rds_path="$1" out_csv="$2"
    Rscript - <<REOF
local({
  renv_lib <- file.path(getwd(), "renv", "library", "macos", "R-4.4",
                        "aarch64-apple-darwin20.0.0")
  if (dir.exists(renv_lib)) .libPaths(c(renv_lib, .libPaths()))
})
suppressPackageStartupMessages(library(Seurat))
obj  <- readRDS("${rds_path}")
meta <- obj@meta.data
cols <- intersect(c("patient","cluster","predicted.id"), colnames(meta))
out  <- meta[, cols, drop=FALSE]
colnames(out)[colnames(out) == "cluster"]     <- "true_label"
colnames(out)[colnames(out) == "predicted.id"] <- "pred_label"
write.csv(out, "${out_csv}", row.names=TRUE)
cat(sprintf("CCA labels saved: %d cells -> %s\n", nrow(out), "${out_csv}"))
REOF
}

# ── Extract labels from scVI h5ad ────────────────────────────────────────────
extract_scvi_labels() {
    local h5ad_path="$1" out_csv="$2"
    python3 -c "
import anndata as ad, pandas as pd, sys
adata = ad.read_h5ad('${h5ad_path}')
cols = [c for c in ['patient','cluster','predicted_label'] if c in adata.obs.columns]
out = adata.obs[cols].copy()
out.columns = [c if c != 'predicted_label' else 'pred_label' for c in out.columns]
out.columns = [c if c != 'cluster' else 'true_label' for c in out.columns]
out.to_csv('${out_csv}')
print(f'scVI labels saved: {len(out)} cells -> ${out_csv}')
"
}

# ── Main LOO loop ─────────────────────────────────────────────────────────────
TOTAL=${#PATIENTS[@]}
FOLD=0
FAILED_FOLDS=()

for PATIENT in "${PATIENTS[@]}"; do
    FOLD=$(( FOLD + 1 ))
    FOLD_STAMP="loo_${PATIENT}"
    CCA_CSV="${LOO_DIR}/fold_${PATIENT}_cca.csv"
    SCVI_CSV="${LOO_DIR}/fold_${PATIENT}_scvi.csv"

    log "═══════════════════════════════════════════════════════"
    log "Fold ${FOLD}/${TOTAL} — held-out patient: ${PATIENT}"
    log "═══════════════════════════════════════════════════════"

    FOLD_START=$(date +%s)

    # ── CCA pipeline ──────────────────────────────────────────────────────────
    log "[CCA] Step 0 — data acquisition (patient split: query=${PATIENT})"
    run_notebook "cca/00_data_acquisition.ipynb" ir \
        "raw_counts_path ${RAW_COUNTS_PATH}" \
        "raw_metadata_path ${RAW_METADATA_PATH}" \
        "raw_tcr_path ${RAW_TCR_PATH}" \
        "caushi_cd8_gz ${CAUSHI_CD8_GZ}" \
        "output_query data/cca/query_${PATIENT}.rds" \
        "output_reference data/cca/reference_no_${PATIENT}.rds" \
        "output_tcr data/cca/yost2019_bcc_tcr.rds" \
        "output_caushi_clonotypes ${CAUSHI_CLONOTYPES_PATH}" \
        "query_patient ${PATIENT}" \
        "n_dims ${N_DIMS}" \
        2>&1 | tee -a "$LOG_FILE"

    log "[CCA] Step 1 — preprocessing"
    run_notebook "cca/01_preprocessing.ipynb" ir \
        "query_path data/cca/query_${PATIENT}.rds" \
        "min_genes ${MIN_GENES}" \
        "max_genes ${MAX_GENES}" \
        "max_pct_mito ${MAX_PCT_MITO}" \
        "output_path data/cca/query_${PATIENT}.rds" \
        2>&1 | tee -a "$LOG_FILE"

    log "[CCA] Step 2 — reference projection"
    run_notebook "cca/02_reference_projection.ipynb" ir \
        "reference_path data/cca/reference_no_${PATIENT}.rds" \
        "query_path data/cca/query_${PATIENT}.rds" \
        "output_path data/cca/query_${PATIENT}_projected.rds" \
        2>&1 | tee -a "$LOG_FILE"

    log "[CCA] Extracting labels -> ${CCA_CSV}"
    extract_cca_labels "data/cca/query_${PATIENT}_projected.rds" "$CCA_CSV" \
        2>&1 | tee -a "$LOG_FILE"

    # ── scVI pipeline ─────────────────────────────────────────────────────────
    log "[scVI] Step 0 — convert RDS -> h5ad"
    run_notebook "scvi/00_convert.ipynb" python3 \
        "reference_rds data/cca/reference_no_${PATIENT}.rds" \
        "query_rds     data/cca/query_${PATIENT}.rds" \
        "output_ref_h5ad   data/scvi/reference_${PATIENT}.h5ad" \
        "output_query_h5ad data/scvi/query_${PATIENT}.h5ad" \
        2>&1 | tee -a "$LOG_FILE"

    log "[scVI] Step 1 — train scVI + scANVI"
    run_notebook "scvi/01_train_reference.ipynb" python3 \
        "reference_h5ad data/scvi/reference_${PATIENT}.h5ad" \
        "model_dir      models/scanvi_reference_${PATIENT}" \
        "labels_key     cluster" \
        "n_epochs_scvi  400" \
        "n_epochs_scanvi 20" \
        2>&1 | tee -a "$LOG_FILE"

    log "[scVI] Step 2 — project query"
    run_notebook "scvi/02_project_query.ipynb" python3 \
        "query_h5ad             data/scvi/query_${PATIENT}.h5ad" \
        "reference_h5ad_trained data/scvi/reference_${PATIENT}_trained.h5ad" \
        "model_dir              models/scanvi_reference_${PATIENT}" \
        "output_projected_h5ad  data/scvi/query_${PATIENT}_projected.h5ad" \
        "USE_SURGERY            False" \
        2>&1 | tee -a "$LOG_FILE"

    log "[scVI] Extracting labels -> ${SCVI_CSV}"
    extract_scvi_labels "data/scvi/query_${PATIENT}_projected.h5ad" "$SCVI_CSV" \
        2>&1 | tee -a "$LOG_FILE"

    FOLD_ELAPSED=$(( $(date +%s) - FOLD_START ))
    log "Fold ${PATIENT} complete in $((FOLD_ELAPSED/60))m $((FOLD_ELAPSED%60))s"
    log "  CCA:  ${CCA_CSV}"
    log "  scVI: ${SCVI_CSV}"
done

# ── Summary ───────────────────────────────────────────────────────────────────
log "═══════════════════════════════════════════════════════"
log "LOO benchmark complete."
log "Results in: ${LOO_DIR}/"
ls -lh "${LOO_DIR}"/fold_*.csv 2>/dev/null || warn "No fold CSVs found."

if [[ ${#FAILED_FOLDS[@]} -gt 0 ]]; then
    warn "Failed folds: ${FAILED_FOLDS[*]}"
else
    log "All ${TOTAL} folds succeeded."
    log ""
    log "Next step: open benchmarking/benchmark_celltype.ipynb and run §4 (LOO section)"
    log "or run:"
    log "  conda run -n scrnaseq papermill benchmarking/benchmark_celltype.ipynb \\"
    log "    benchmarking/benchmark_celltype_executed.ipynb \\"
    log "    -k python3 -p loo_dir data/loo"
fi
