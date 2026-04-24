#!/usr/bin/env bash
# =============================================================================
# run_cca_pipeline.sh — CCA / Seurat Reference Projection Pipeline
# =============================================================================
# Method: Seurat CCA label transfer (R/Seurat)
# Notebooks 00–03: R / Seurat (ir kernel) — executed via Papermill + IRkernel
# Notebooks 04–05: Python / scikit-learn (python3 kernel) — shared with scVI
# Outputs: data/cca/
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG="${SCRIPT_DIR}/config/params.yaml"
NB_DIR="${SCRIPT_DIR}/notebooks"
OUT_DIR="${NB_DIR}/executed"
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")

# Colours
GREEN='\033[0;32m'; YELLOW='\033[1;33m'; RED='\033[0;31m'; CYAN='\033[0;36m'; NC='\033[0m'

log()  { echo -e "${GREEN}[$(date '+%H:%M:%S')] $*${NC}"; }
warn() { echo -e "${YELLOW}[$(date '+%H:%M:%S')] WARNING: $*${NC}"; }
fail() { echo -e "${RED}[$(date '+%H:%M:%S')] FAILED: $*${NC}"; exit 1; }
mem()  { echo -e "${CYAN}[$(date '+%H:%M:%S')] [MEM] $*${NC}"; }

# ── Memory monitor (background) ───────────────────────────────────────────────
# Prints RAM free/used every MONITOR_INTERVAL seconds to stdout + log file.
MONITOR_INTERVAL=30
MONITOR_LOG="${SCRIPT_DIR}/data/memory_monitor.log"
_MONITOR_PID=""

mem_snapshot() {
    local total_bytes free_bytes used_bytes total_gb free_gb used_gb pct_used
    total_bytes=$(sysctl -n hw.memsize 2>/dev/null || echo 0)
    # vm_stat page size and free+inactive pages (available to apps on macOS)
    local vm; vm=$(vm_stat 2>/dev/null)
    local page_size; page_size=$(echo "$vm" | awk '/page size of/ {print $NF}')
    local pages_free; pages_free=$(echo "$vm" | awk '/^Pages free:/ {gsub(/\./,"",$3); print $3}')
    local pages_inactive; pages_inactive=$(echo "$vm" | awk '/^Pages inactive:/ {gsub(/\./,"",$3); print $3}')
    local pages_wired; pages_wired=$(echo "$vm" | awk '/^Pages wired down:/ {gsub(/\./,"",$4); print $4}')
    local pages_active; pages_active=$(echo "$vm" | awk '/^Pages active:/ {gsub(/\./,"",$3); print $3}')
    page_size=${page_size:-16384}
    pages_free=${pages_free:-0}; pages_inactive=${pages_inactive:-0}
    pages_wired=${pages_wired:-0}; pages_active=${pages_active:-0}
    free_bytes=$(( (pages_free + pages_inactive) * page_size ))
    used_bytes=$(( (pages_wired + pages_active) * page_size ))
    total_gb=$(awk "BEGIN {printf \"%.1f\", $total_bytes/1073741824}")
    free_gb=$(awk  "BEGIN {printf \"%.1f\", $free_bytes/1073741824}")
    used_gb=$(awk  "BEGIN {printf \"%.1f\", $used_bytes/1073741824}")
    pct_used=$(awk "BEGIN {printf \"%.0f\", $used_bytes/$total_bytes*100}" 2>/dev/null || echo "?")
    echo "total=${total_gb}GB  used=${used_gb}GB (${pct_used}%)  free=${free_gb}GB"
}

start_memory_monitor() {
    mkdir -p "$(dirname "$MONITOR_LOG")"
    echo "timestamp,total_gb,used_gb,free_gb,pct_used" > "$MONITOR_LOG"
    (
        while true; do
            local snap; snap=$(mem_snapshot)
            local ts; ts=$(date '+%H:%M:%S')
            # Append CSV row
            local total used free pct
            total=$(echo "$snap" | grep -o 'total=[^G]*' | cut -d= -f2)
            used=$(echo  "$snap" | grep -o 'used=[^G]*'  | cut -d= -f2)
            free=$(echo  "$snap" | grep -o 'free=[^G]*'  | cut -d= -f2)
            pct=$(echo   "$snap" | grep -oE '\([0-9]+%\)' | tr -d '()%')
            echo "$ts,$total,$used,$free,$pct" >> "$MONITOR_LOG"
            # Print coloured line to the same stdout as the pipeline
            echo -e "${CYAN}[${ts}] [MEM] ${snap}${NC}"
            sleep "$MONITOR_INTERVAL"
        done
    ) &
    _MONITOR_PID=$!
}

stop_memory_monitor() {
    if [[ -n "$_MONITOR_PID" ]]; then
        kill "$_MONITOR_PID" 2>/dev/null || true
        wait "$_MONITOR_PID" 2>/dev/null || true
        _MONITOR_PID=""
        mem "Monitor stopped. Full log: ${MONITOR_LOG}"
    fi
}
# Ensure monitor is killed if the script exits for any reason
trap 'stop_memory_monitor' EXIT

# ── Pre-flight checks ──────────────────────────────────────────────────────────
command -v papermill >/dev/null 2>&1 || fail "papermill not found. Install with: pip install papermill"
command -v python    >/dev/null 2>&1 || fail "python not found"
command -v Rscript   >/dev/null 2>&1 || fail "Rscript not found. Install R and run: IRkernel::installspec()"
[ -f "$CONFIG" ] || fail "Config file not found: $CONFIG"

mkdir -p "$OUT_DIR"
log "Pipeline start — timestamp: ${TIMESTAMP}"
log "Config: ${CONFIG}"
log "Output dir: ${OUT_DIR}"

# ── Helper: read a scalar value from params.yaml ───────────────────────────────
yaml_val() {
    python -c "import yaml; d=yaml.safe_load(open('$CONFIG')); print(d.get('$1', ''))"
}

# ── Read parameters from config/params.yaml ───────────────────────────────────
# Raw input files
RAW_COUNTS_PATH=$(yaml_val raw_counts_path)
RAW_METADATA_PATH=$(yaml_val raw_metadata_path)
RAW_TCR_PATH=$(yaml_val raw_tcr_path)
CAUSHI_CD8_GZ=$(yaml_val caushi_cd8_gz)

# Derived paths
REFERENCE_PATH=$(yaml_val reference_path)
QUERY_PATH=$(yaml_val query_path)
TCR_PATH=$(yaml_val tcr_path)
CAUSHI_CLONOTYPES_PATH=$(yaml_val caushi_clonotypes_path)
REFERENCE_SPLIT=$(yaml_val reference_split_fraction)

# QC thresholds
MIN_GENES=$(yaml_val min_genes)
MAX_GENES=$(yaml_val max_genes)
MAX_PCT_MITO=$(yaml_val max_pct_mito)

# Analysis parameters
N_DIMS=$(yaml_val n_dims)
MIN_CELLS_PER_CLONOTYPE=$(yaml_val min_cells_per_clonotype)
N_EXHAUSTED_CLONOTYPES=$(yaml_val n_exhausted_clonotypes)
QUERY_PATIENT=$(yaml_val query_patient)
# Ensure a concrete patient is set (required for output file naming)
[[ -z "${QUERY_PATIENT}" ]] && { QUERY_PATIENT="su009"; warn "query_patient not set in config — defaulting to ${QUERY_PATIENT}"; }

# ── Patient-stamped data paths ─────────────────────────────────────────────────
# Every output file encodes the held-out patient ID so different LOO folds
# can coexist on disk without overwriting each other.
CCA_REF="data/cca/reference_no_${QUERY_PATIENT}.rds"
CCA_QRY="data/cca/query_${QUERY_PATIENT}.rds"
CCA_QRY_PROJ="data/cca/query_${QUERY_PATIENT}_projected.rds"
CCA_QRY_CLONO="data/cca/query_${QUERY_PATIENT}_with_clonotypes.rds"
CCA_RANKED="data/cca/ranked_clonotypes_${QUERY_PATIENT}.csv"
CCA_SELECTED="data/cca/selected_tumor_reactive_tcrs_${QUERY_PATIENT}.csv"
CCA_PPV="data/cca/ppv_summary_${QUERY_PATIENT}.csv"

log "Parameters loaded:"
log "  raw_counts=${RAW_COUNTS_PATH}"
log "  raw_metadata=${RAW_METADATA_PATH}"
log "  raw_tcr=${RAW_TCR_PATH}"
log "  caushi_cd8=${CAUSHI_CD8_GZ}"
log "  reference_path=${REFERENCE_PATH}  query_path=${QUERY_PATH}"
log "  tcr_path=${TCR_PATH}  n_dims=${N_DIMS}"
log "  min_genes=${MIN_GENES}  max_genes=${MAX_GENES}  max_pct_mito=${MAX_PCT_MITO}"
log "  min_cells_per_clonotype=${MIN_CELLS_PER_CLONOTYPE}  n_exhausted_clonotypes=${N_EXHAUSTED_CLONOTYPES}"
log "  query_patient=${QUERY_PATIENT}  (held-out → query; rest → reference)"

# ── Start memory monitor ───────────────────────────────────────────────────────
mem "Initial: $(mem_snapshot)"
start_memory_monitor
log "Memory monitor started (every ${MONITOR_INTERVAL}s → ${MONITOR_LOG})"

# ── Run pipeline ───────────────────────────────────────────────────────────────
run_notebook() {
    local nb_name="$1"
    local kernel="$2"       # 'ir' for R notebooks, 'python3' for Python
    # Optional 3rd arg: output sub-path (defaults to nb_name).
    # Pass e.g. "cca/04_tcr_reactivity_selection.ipynb" so shared notebooks
    # land inside the cca/ executed subfolder instead of the flat root.
    local out_name
    if [[ "$3" == *".ipynb" && "$3" != -* ]]; then
        out_name="$3"; shift 3
    else
        out_name="${nb_name}"; shift 2
    fi
    local extra_params=("$@")

    local input="${NB_DIR}/${nb_name}"
    local output="${OUT_DIR}/${TIMESTAMP}_query-${QUERY_PATIENT}_${out_name}"
    mkdir -p "$(dirname "$output")"

    [ -f "$input" ] || fail "Notebook not found: $input"

    mem "Before ${nb_name}: $(mem_snapshot)"
    log "Running ${nb_name} (kernel: ${kernel}) ..."

    # Build parameter string
    local param_str=""
    for p in "${extra_params[@]}"; do
        param_str="${param_str} -p ${p}"
    done

    # shellcheck disable=SC2086
    papermill \
        "$input" \
        "$output" \
        --kernel "${kernel}" \
        --report-mode \
        $param_str \
        2>&1 | tee -a "${OUT_DIR}/${TIMESTAMP}_query-${QUERY_PATIENT}_pipeline.log"

    mem "After  ${nb_name}: $(mem_snapshot)"
    log "${nb_name} complete → ${output}"
}

# ── Run pipeline ───────────────────────────────────────────────────────────────

# Step 0: Build all derived data files from raw GEO downloads
run_notebook "cca/00_data_acquisition.ipynb" ir \
    "raw_counts_path ${RAW_COUNTS_PATH}" \
    "raw_metadata_path ${RAW_METADATA_PATH}" \
    "raw_tcr_path ${RAW_TCR_PATH}" \
    "caushi_cd8_gz ${CAUSHI_CD8_GZ}" \
    "output_query ${CCA_QRY}" \
    "output_reference ${CCA_REF}" \
    "output_tcr ${TCR_PATH}" \
    "output_caushi_clonotypes ${CAUSHI_CLONOTYPES_PATH}" \
    "query_patient ${QUERY_PATIENT}" \
    "n_dims ${N_DIMS}"

# Step 1: QC and preprocessing on query cells
run_notebook "cca/01_preprocessing.ipynb" ir \
    "query_path ${CCA_QRY}" \
    "min_genes ${MIN_GENES}" \
    "max_genes ${MAX_GENES}" \
    "max_pct_mito ${MAX_PCT_MITO}" \
    "output_path ${CCA_QRY}"

# Step 2: CCA projection — maps query cells onto reference T cell atlas
run_notebook "cca/02_reference_projection.ipynb" ir \
    "reference_path ${CCA_REF}" \
    "query_path ${CCA_QRY}" \
    "output_path ${CCA_QRY_PROJ}"

# Step 3: Link TCR clonotypes to exhaustion scores
run_notebook "cca/03_clonotype_exhaustion.ipynb" ir \
    "query_path ${CCA_QRY_PROJ}" \
    "exhaustion_panel_path data/exhaustion_gene_panel.txt" \
    "output_clonotype_table ${CCA_RANKED}" \
    "output_adata ${CCA_QRY_CLONO}" \
    "tcr_path ${TCR_PATH}" \
    "min_cells_per_clonotype ${MIN_CELLS_PER_CLONOTYPE}" \
    "n_exhausted_clonotypes ${N_EXHAUSTED_CLONOTYPES}"

# Step 4: RF classifier to select tumour-reactive TCRs
# Output goes into the cca/ subfolder so it isn't mixed with scVI executed notebooks
run_notebook "04_tcr_reactivity_selection.ipynb" python3 "cca/04_tcr_reactivity_selection.ipynb" \
    "clonotype_table_path ${CCA_RANKED}" \
    "caushi_table_path ${CAUSHI_CLONOTYPES_PATH}" \
    "model_path models/rf_tumor_reactivity.pkl" \
    "output_selected_tcr ${CCA_SELECTED}"

# Step 5: PPV validation — precision/recall of selected TCRs
run_notebook "05_ppv_validation.ipynb" python3 "cca/05_ppv_validation.ipynb" \
    "selected_tcr_path ${CCA_SELECTED}" \
    "output_ppv_table ${CCA_PPV}"

# ── Done ───────────────────────────────────────────────────────────────────────
stop_memory_monitor
log "═══════════════════════════════════════════════════════"
log "Pipeline complete. Query patient: ${QUERY_PATIENT}"
log "Executed notebooks: ${OUT_DIR}/${TIMESTAMP}_query-${QUERY_PATIENT}_cca/"
log "Log: ${OUT_DIR}/${TIMESTAMP}_query-${QUERY_PATIENT}_pipeline.log"
log "Data outputs: data/cca/*_${QUERY_PATIENT}.*"

# Print PPV summary if available
PPV_FILE="${SCRIPT_DIR}/${CCA_PPV}"
if [ -f "$PPV_FILE" ]; then
    log "PPV Summary:"
    python -c "
import pandas as pd
df = pd.read_csv('${PPV_FILE}')
print(df.to_string(index=False))
"
fi
