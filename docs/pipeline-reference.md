# Pipeline Reference

Detailed per-step logic, data flow, technical stack, and extension guide.

---

## Data Flow

Both pipelines share the same raw inputs (NB00) and the same downstream steps (NB04–05). Only the cell-type projection (NB01–NB03) differs.

```
data/raw_downloads/
        │
        ▼
[CCA]  notebooks/cca/00_data_acquisition   →  data/cca/reference.rds
                                               data/cca/query.rds
        │
        ├── notebooks/cca/01_preprocessing       →  data/cca/query.rds (QC'd)
        ├── notebooks/cca/02_reference_projection →  data/cca/query_projected.rds
        └── notebooks/cca/03_clonotype_exhaustion →  data/cca/ranked_clonotypes.csv
                                                      data/cca/query_with_clonotypes.rds

[scVI] notebooks/scvi/00_convert           →  data/scvi/reference.h5ad
                                               data/scvi/query.h5ad
        │
        ├── notebooks/scvi/01_train_reference    →  models/scanvi_reference/
        ├── notebooks/scvi/02_project_query      →  data/scvi/query_projected.h5ad
        └── notebooks/scvi/03_clonotype_exhaustion → data/scvi/ranked_clonotypes.csv
                                                      data/scvi/query_with_clonotypes.h5ad

        ▼ (either pipeline feeds here)
notebooks/04_tcr_reactivity_selection  →  data/{cca,scvi}/selected_tumor_reactive_tcrs.csv
notebooks/05_ppv_validation            →  data/{cca,scvi}/ppv_summary.csv
```

---

## Step-by-step Logic

### NB00 · Data Acquisition (CCA; scVI reads its outputs)
- Downloads Yost 2019 BCC data from GEO; 70/30 random split into reference / query Seurat objects
- **Outputs:** `data/cca/reference.rds`, `data/cca/query.rds`, `data/cca/yost2019_bcc_tcr.rds`
- scVI pipeline converts these via `notebooks/scvi/00_convert.ipynb` (rpy2 → AnnData h5ad via MEX)

### CCA — NB01 · Preprocessing
- QC: `nFeature_RNA`, `percent.mt` thresholds → `NormalizeData` → `FindVariableFeatures` (2000 HVGs) → `ScaleData` → `RunPCA`
- **Output:** `data/cca/query.rds`  |  **Papermill params:** `min_genes`, `max_genes`, `max_pct_mito`

### CCA — NB02 · Reference Projection
- `FindTransferAnchors(reduction = "pcaproject", reference.reduction = "pca")` + `MapQuery`
- Transfers cell-type labels; projects query onto reference UMAP
- **Output:** `data/cca/query_projected.rds`  |  **Papermill params:** `reference_path`, `query_path`, `n_dims`

### scVI — NB01 · Train Reference
- Trains SCVI on reference for dimensionality reduction, then SCANVI for semi-supervised label transfer
- `batch_size=512`, `accelerator="mps"` (Apple Silicon) — change to `"cuda"` on GPU servers
- **Output:** `models/scanvi_reference/`

### scVI — NB02 · Project Query
- Pure inference: loads trained SCANVI, projects query cells into the reference latent space
- Predicts cell-type labels and per-cell soft probabilities
- **Output:** `data/scvi/query_projected.h5ad`

### CCA — NB03 · Clonotype × Exhaustion
- Merges TCR clonotype table; computes exhaustion module score via `AddModuleScore` (18-gene panel; see [exhaustion-gene-panel.md](exhaustion-gene-panel.md))
- Ranks clonotypes by composite score: `mean_exhaustion × log1p(clone_size)`
- **Outputs:** `data/cca/ranked_clonotypes.csv`, `data/cca/query_with_clonotypes.rds`

### scVI — NB03 · Clonotype × Exhaustion
- Mirror of CCA NB03 on AnnData; gene scoring via `sc.tl.score_genes`
- Produces identical output schema for seamless handoff to NB04
- **Outputs:** `data/scvi/ranked_clonotypes.csv`, `data/scvi/query_with_clonotypes.h5ad`

### NB04 · TCR Reactivity Selection *(placeholder)*
- **Current implementation:** Random Forest (`class_weight="balanced"`, 5-fold CV) on exhaustion score, clone frequency, CDR3 length, TRBV family, rank score
- This is one candidate; `benchmarking/benchmark_tcr_ml.ipynb` compares LR, RF, and manual ranking. A definitive method will be selected once labelled training data is available.
- **Input:** `ranked_clonotypes.csv` from either pipeline  |  **Output:** `selected_tumor_reactive_tcrs.csv`

### NB05 · PPV Validation *(placeholder)*
- Computes Positive Predictive Value against an external `is_tumor_reactive` ground-truth column
- Framework is in place; PPV figures will be meaningful once a validated reactivity dataset is supplied
- **Output:** `data/{cca,scvi}/ppv_summary.csv`

---

## Technical Stack

| Step | Method | Language | Technology |
|------|--------|----------|------------|
| Data acquisition (NB00) | CCA | R | Seurat, GEOquery |
| Preprocessing (NB01) | CCA | R | [Seurat](https://satijalab.org/seurat/) |
| Reference projection (NB02) | CCA | R | Seurat `FindTransferAnchors` + `MapQuery` |
| Clonotype × exhaustion (NB03) | CCA | R | Seurat `AddModuleScore`, scRepertoire |
| Format conversion (NB00) | scVI | Python | rpy2, anndata |
| Model training (NB01) | scVI | Python | [scvi-tools](https://scvi-tools.org) SCANVI |
| Query projection (NB02) | scVI | Python | scvi-tools |
| Clonotype × exhaustion (NB03) | scVI | Python | scanpy `sc.tl.score_genes` |
| TCR selection (NB04) | Shared | Python | [scikit-learn](https://scikit-learn.org) Random Forest |
| PPV validation (NB05) | Shared | Python | pandas / numpy |
| Orchestration | Both | bash | [Papermill](https://papermill.readthedocs.io) |

---

## Adding a Third Method (e.g. scGPT)

1. Create `notebooks/scgpt/` with method-specific notebooks `00`–`03`
2. Create `data/scgpt/` for its outputs
3. Copy `run_scvi_pipeline.sh` → `run_scgpt_pipeline.sh`; update notebook paths and output directory to `data/scgpt/`
4. NB04 and NB05 are already shared — pass `-p clonotype_table_path data/scgpt/ranked_clonotypes.csv` via the pipeline script
