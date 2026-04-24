#!/usr/bin/env bash
# =============================================================================
# run_scvi_pipeline.sh — scVI / scANVI Reference Projection Pipeline
# =============================================================================
# Method: scVI variational autoencoder + scANVI label transfer (Python)
# Notebooks 00–03: Python / scvi-tools (python3 kernel) — scVI-specific
# Notebooks 04–05: Python / scikit-learn (python3 kernel) — shared with CCA
# Prerequisite: data/cca/reference.rds and data/cca/query.rds must exist
#               (produced by run_cca_pipeline.sh NB00, or run NB00 standalone).
# Outputs: data/scvi/
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

TIMESTAMP=$(date +%Y%m%d_%H%M%S)
OUT_DIR="notebooks/executed"
mkdir -p "$OUT_DIR"

log() { echo "[$(date +%H:%M:%S)] $*"; }

# ── Read query patient from config ────────────────────────────────────────────
CONFIG="${SCRIPT_DIR}/config/params.yaml"
yaml_val() { python3 -c "import yaml; d=yaml.safe_load(open('$CONFIG')); print(d.get('$1',''))"; }
QUERY_PATIENT=$(yaml_val query_patient)
[[ -z "${QUERY_PATIENT}" ]] && QUERY_PATIENT="su009"
log "Query patient: ${QUERY_PATIENT}"

# Patient-stamped paths (must match CCA pipeline for the same fold)
CCA_REF="data/cca/reference_no_${QUERY_PATIENT}.rds"
CCA_QRY="data/cca/query_${QUERY_PATIENT}.rds"
SCVI_REF="data/scvi/reference_${QUERY_PATIENT}.h5ad"
SCVI_QRY="data/scvi/query_${QUERY_PATIENT}.h5ad"
SCVI_QRY_PROJ="data/scvi/query_${QUERY_PATIENT}_projected.h5ad"
SCVI_QRY_CLONO="data/scvi/query_${QUERY_PATIENT}_with_clonotypes.h5ad"
SCVI_RANKED="data/scvi/ranked_clonotypes_${QUERY_PATIENT}.csv"
SCVI_SELECTED="data/scvi/selected_tumor_reactive_tcrs_${QUERY_PATIENT}.csv"
SCVI_PPV="data/scvi/ppv_summary_${QUERY_PATIENT}.csv"

run_nb() {
    local nb="$1" kernel="$2"; shift 2
    # Accept absolute-ish paths (starting with notebooks/) or scvi-relative names
    if [[ "$nb" == notebooks/* ]]; then
        local input="$nb"
    else
        local input="notebooks/scvi/${nb}"
    fi
    local basename="${nb##*/}"
    local output="${OUT_DIR}/${TIMESTAMP}_query-${QUERY_PATIENT}_scvi/${basename}"
    mkdir -p "${OUT_DIR}/${TIMESTAMP}_query-${QUERY_PATIENT}_scvi"
    local params=""
    for kv in "$@"; do params="$params -p ${kv%% *} ${kv#* }"; done
    log "Running ${basename} (kernel: ${kernel}) ..."
    # shellcheck disable=SC2086
    conda run -n scrnaseq papermill "$input" "$output" -k "$kernel" $params 2>&1
    log "${basename} complete -> ${output}"
}

# ── Prevent macOS ARM OpenMP duplicate-library crash ─────────────────────────
# Pre-load the canonical libomp so dyld deduplicates all subsequent loads
export DYLD_INSERT_LIBRARIES=/opt/anaconda3/envs/scrnaseq/lib/libomp.dylib
export KMP_DUPLICATE_LIB_OK=TRUE
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

# ── Install scvi-tools if not present ─────────────────────────────────────────
conda run -n scrnaseq python -c "import scvi" 2>/dev/null || {
    log "Installing scvi-tools >= 1.1.0 ..."
    conda run -n scrnaseq pip install "scvi-tools>=1.1.0" "scanpy>=1.9" --quiet
    log "Installation complete."
}

# ── Step 0: Convert RDS -> h5ad ───────────────────────────────────────────────
run_nb "00_convert.ipynb" python3 \
    "reference_rds ${CCA_REF}" \
    "query_rds     ${CCA_QRY}" \
    "output_ref_h5ad   ${SCVI_REF}" \
    "output_query_h5ad ${SCVI_QRY}"

# ── Step 1: Train SCVI + SCANVI on reference ──────────────────────────────────
run_nb "01_train_reference.ipynb" python3 \
    "reference_h5ad ${SCVI_REF}" \
    "model_dir      models/scanvi_reference_${QUERY_PATIENT}" \
    "labels_key     cluster" \
    "n_epochs_scvi  400" \
    "n_epochs_scanvi 20"

# ── Step 2: Project query (pure inference, no retraining) ─────────────────────
run_nb "02_project_query.ipynb" python3 \
    "query_h5ad            ${SCVI_QRY}" \
    "model_dir             models/scanvi_reference_${QUERY_PATIENT}" \
    "output_projected_h5ad ${SCVI_QRY_PROJ}" \
    "USE_SURGERY           False"

# ── Step 3: Clonotype integration & exhaustion scoring ────────────────────────
run_nb "03_clonotype_exhaustion.ipynb" python3 \
    "query_h5ad             ${SCVI_QRY_PROJ}" \
    "exhaustion_panel_path  data/exhaustion_gene_panel.txt" \
    "output_clonotype_table ${SCVI_RANKED}" \
    "output_h5ad            ${SCVI_QRY_CLONO}" \
    "min_cells_per_clonotype 2" \
    "n_exhausted_clonotypes  20"

# ── Step 4: TCR reactivity selection ─────────────────────────────────────────
run_nb "notebooks/04_tcr_reactivity_selection.ipynb" python3 \
    "clonotype_table_path ${SCVI_RANKED}" \
    "output_selected_tcr  ${SCVI_SELECTED}"

# ── Step 5: PPV validation ────────────────────────────────────────────────────
run_nb "notebooks/05_ppv_validation.ipynb" python3 \
    "selected_tcr_path ${SCVI_SELECTED}" \
    "output_ppv_table  ${SCVI_PPV}"

log "scANVI pipeline complete. Query patient: ${QUERY_PATIENT}"
log "Executed notebooks: ${OUT_DIR}/${TIMESTAMP}_query-${QUERY_PATIENT}_scvi/"
log "Data outputs: data/scvi/*_${QUERY_PATIENT}.*"
