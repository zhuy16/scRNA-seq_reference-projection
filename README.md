# scRNA-seq Reference Projection Pipeline

A reproducible, [Papermill](https://papermill.readthedocs.io)-parameterised pipeline for single-cell RNA-seq analysis focusing on **tumor-infiltrating lymphocytes (TIL)** and **CD8+ exhausted T cell** identification for tumor-reactive TCR selection.

Two parallel preprocessing methods are implemented and can be run independently or compared:

| Method | Language | Script |
|--------|----------|--------|
| **CCA** — Seurat Canonical Correlation Analysis | R | `run_cca_pipeline.sh` |
| **scVI/scANVI** — variational autoencoder | Python | `run_scvi_pipeline.sh` |

---

## Demo Dataset

This pipeline uses the **Yost et al. 2019 BCC dataset** ([GSE123813](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE123813), GEO open access) as a fully reproducible stand-in for clinical TIL data.

**Why this dataset:**
- Paired **scRNA-seq + TCR-seq** from 11 basal cell carcinoma (BCC) patients (su001–su012, no su011), both pre- and post-anti-PD1 immunotherapy
- Includes annotated **exhausted CD8+ T cell subtypes** (`CD8_ex`, `CD8_ex_act`) that serve as proxies for chronic antigen stimulation by tumor neoantigens
- Multi-patient cohort enables leave-one-out cross-validation to assess projection robustness across individuals
- Fully open access — no institutional data transfer agreement required

**Limitation:** Yost 2019 does not include a direct tumor-reactivity assay. Tumor reactivity is inferred from exhaustion state and clonal expansion. PPV validation (NB05) therefore cross-references predicted TCRs against the independent **Caushi et al. 2021** neoantigen reactivity ground truth (lung cancer; see [Citation / Acknowledgements](#citation--acknowledgements)).

> In a production setting, replace the raw GEO inputs with in-house TIL scRNA-seq data. All paths are configured in `config/params.yaml`.

**Reference:** Yost KE, Satpathy AT, Wells DK, et al. *Clonal replacement of tumor-specific T cells following PD-1 blockade.* Nature Medicine. 2019;25(8):1251–1259. https://doi.org/10.1038/s41591-019-0522-3

---

## Clinical Context

Adoptive cell therapy with tumor-infiltrating lymphocytes (TIL) requires identifying the small fraction of CD8+ T cells that are genuinely tumor-reactive. This pipeline automates:

1. **QC and preprocessing** of scRNA-seq data from patient biopsies
2. **Reference projection** — assigns cell-type labels to new samples without manual re-annotation using two independent methods
3. **Clonotype × exhaustion integration** — identifies expanded, antigen-experienced CD8+ T cell clones
4. **TCR reactivity selection** — selects top tumor-reactive TCRs using a trained Random Forest
5. **PPV validation** against neoantigen reactivity ground truth

---

## Repository Structure

```
scRNA-seq_reference-projection/
├── run_cca_pipeline.sh          # CCA pipeline: NB00–03 (R) + NB04–05 (Python)
├── run_scvi_pipeline.sh         # scVI pipeline: NB00–03 (Python) + NB04–05 (Python)
├── run_loo_benchmark.sh         # Leave-one-out patient benchmark (CCA + scVI)
├── Makefile                     # Shortcuts: make setup / run-cca / run-scvi / run-loo / …
├── monitor_live.py              # Live CPU/RAM/MPS dashboard (run in a second terminal)
├── Dockerfile                   # Self-contained Linux image (alternative to conda)
├── .dockerignore
├── .gitignore
├── config/
│   └── params.yaml              # Shared parameters (single source of truth)
│
├── notebooks/
│   ├── cca/                     # R/Seurat notebooks (CCA method)
│   │   ├── 00_data_acquisition.ipynb      # Download + split into reference/query
│   │   ├── 01_preprocessing.ipynb         # QC, normalisation, HVG, PCA
│   │   ├── 02_reference_projection.ipynb  # CCA label transfer + MapQuery
│   │   └── 03_clonotype_exhaustion.ipynb  # TCR × exhaustion (R/Seurat)
│   │
│   ├── scvi/                    # Python/scVI notebooks (scVI method)
│   │   ├── 00_convert.ipynb               # RDS → h5ad (reads from data/cca/)
│   │   ├── 01_train_reference.ipynb       # Train SCVI + SCANVI on reference
│   │   ├── 02_project_query.ipynb         # Project query (pure inference)
│   │   └── 03_clonotype_exhaustion.ipynb  # TCR × exhaustion (Python/scanpy)
│   │
│   ├── 04_tcr_reactivity_selection.ipynb  # Shared: Random Forest TCR selection
│   ├── 05_ppv_validation.ipynb            # Shared: PPV against neoantigen ground truth
│   └── executed/                          # Timestamped executed notebooks (gitignored)
│
├── data/
│   ├── raw_downloads/           # Raw GEO input files (gitignored)
│   ├── cca/                     # CCA pipeline outputs
│   │   ├── reference.rds                  # Reference Seurat object (NB00)
│   │   ├── query.rds                      # Query Seurat object (NB00/01)
│   │   ├── yost2019_bcc_tcr.rds           # TCR clonotype data (NB00)
│   │   ├── query_projected.rds            # After CCA label transfer (NB02)
│   │   ├── query_with_clonotypes.rds      # CD8 cells + clonotype metadata (NB03)
│   │   ├── ranked_clonotypes.csv          # Exhaustion-ranked clonotypes (NB03)
│   │   ├── selected_tumor_reactive_tcrs.csv  # Top TCRs (NB04)
│   │   └── ppv_summary.csv                # PPV results (NB05)
│   │
│   ├── scvi/                    # scVI pipeline outputs
│   │   ├── reference.h5ad                 # Reference AnnData (NB00_convert)
│   │   ├── query.h5ad                     # Query AnnData (NB00_convert)
│   │   ├── reference_trained.h5ad         # After SCVI training (NB01)
│   │   ├── query_projected.h5ad           # After SCANVI projection (NB02)
│   │   ├── query_with_clonotypes.h5ad     # CD8 cells + clonotype metadata (NB03)
│   │   ├── ranked_clonotypes.csv          # Exhaustion-ranked clonotypes (NB03)
│   │   ├── selected_tumor_reactive_tcrs.csv  # Top TCRs (NB04)
│   │   └── ppv_summary.csv                # PPV results (NB05)
│   │
│   ├── exhaustion_gene_panel.txt          # 18-gene exhaustion panel (shared)
│   ├── caushi2021_ranked_clonotypes.csv   # Ground-truth neoantigen data (gitignored)
│   ├── caushi2021_cd8_annotated.rds       # Caushi CD8 reference (gitignored)
│   └── loo/                               # LOO per-fold prediction CSVs (gitignored)
│       ├── fold_<patient>_cca.csv         #   true/predicted labels — CCA
│       └── fold_<patient>_scvi.csv        #   true/predicted labels — scVI
│
├── models/
│   ├── scanvi_reference/            # Main SCANVI model weights (gitignored)
│   └── scanvi_reference_su*/        # Per-patient LOO SCANVI models (gitignored)
│
├── benchmarking/
│   ├── benchmark_celltype.ipynb      # CCA vs scVI/scANVI cell-type accuracy
│   ├── benchmark_tcr_ml.ipynb        # Manual rank vs LR vs RF comparison
│   ├── loo_bar_charts.png            # LOO accuracy bar charts (per patient)
│   ├── loo_boxplots.png              # LOO metric distribution boxplots
│   ├── loo_confusion_matrices.png    # Per-method confusion matrices
│   ├── loo_scatter.png               # CCA vs scVI per-patient scatter
│   ├── loo_per_patient_metrics.csv   # Per-patient accuracy / F1 / MCC table
│   └── loo_summary.csv              # Aggregated LOO summary statistics
│
├── environment.yml              # Conda environment (Python + R runtime)
├── setup_r_env.R                # R package installer (writes renv.lock)
└── renv.lock                    # Exact R package versions
```

---

## Data Flow

Both pipelines share the same raw inputs (NB00) and the same downstream analysis (NB04–05). Only the projection step (NB01–03) differs.

```
data/raw_downloads/
        │
        ▼
[CCA]  notebooks/cca/00_data_acquisition   →  data/cca/reference.rds
                                               data/cca/query.rds
        │
        ├──[CCA]  notebooks/cca/01_preprocessing       →  data/cca/query.rds (QC'd)
        ├──[CCA]  notebooks/cca/02_reference_projection →  data/cca/query_projected.rds
        └──[CCA]  notebooks/cca/03_clonotype_exhaustion →  data/cca/ranked_clonotypes.csv
                                                            data/cca/query_with_clonotypes.rds

[scVI] notebooks/scvi/00_convert           →  data/scvi/reference.h5ad
                                               data/scvi/query.h5ad
        │
        ├──[scVI] notebooks/scvi/01_train_reference    →  models/scanvi_reference/
        ├──[scVI] notebooks/scvi/02_project_query      →  data/scvi/query_projected.h5ad
        └──[scVI] notebooks/scvi/03_clonotype_exhaustion → data/scvi/ranked_clonotypes.csv
                                                            data/scvi/query_with_clonotypes.h5ad

        ▼ (either pipeline feeds here)
notebooks/04_tcr_reactivity_selection  →  data/{cca,scvi}/selected_tumor_reactive_tcrs.csv
notebooks/05_ppv_validation            →  data/{cca,scvi}/ppv_summary.csv
```

---

## Pipeline Logic

### Shared — NB00 · Data Acquisition (CCA only; scVI reads its outputs)
- Downloads Yost 2019 BCC data from GEO; splits into reference / query Seurat objects
- **Outputs:** `data/cca/reference.rds`, `data/cca/query.rds`, `data/cca/yost2019_bcc_tcr.rds`
- scVI pipeline reads these via `notebooks/scvi/00_convert.ipynb`

### CCA — NB01 · Preprocessing (R / Seurat)
- QC filtering (`nFeature_RNA`, `percent.mt`) → `NormalizeData` → `FindVariableFeatures` (2000 HVGs) → `ScaleData` → `RunPCA`
- **Output:** `data/cca/query.rds`
- **Papermill params:** `min_genes`, `max_genes`, `max_pct_mito`

### CCA — NB02 · Reference Projection (R / Seurat CCA)
- `FindTransferAnchors(reduction = "pcaproject", reference.reduction = "pca")` + `MapQuery`
- Transfers cell-type labels; projects query onto reference UMAP
- **Output:** `data/cca/query_projected.rds`
- **Papermill params:** `reference_path`, `query_path`, `n_dims`

### scVI — NB00 · Convert (Python)
- Exports R Seurat objects to AnnData h5ad via MEX format
- **Inputs:** `data/cca/reference.rds`, `data/cca/query.rds`
- **Outputs:** `data/scvi/reference.h5ad`, `data/scvi/query.h5ad`

### scVI — NB01 · Train Reference (Python / scvi-tools)
- Trains SCVI on reference, then SCANVI for label transfer
- `batch_size=512`, `accelerator="mps"` (Apple Silicon) — change to `"cuda"` on GPU servers
- **Output:** `models/scanvi_reference/`

### scVI — NB02 · Project Query (Python / scvi-tools)
- Pure inference: projects query cells into the trained SCANVI latent space
- Predicts cell-type labels and soft probabilities
- **Output:** `data/scvi/query_projected.h5ad`

### Shared — NB03 · Clonotype × Exhaustion Integration
- **CCA version** (`notebooks/cca/03_clonotype_exhaustion.ipynb`): R/Seurat `AddModuleScore`
- **scVI version** (`notebooks/scvi/03_clonotype_exhaustion.ipynb`): Python `sc.tl.score_genes`
- Both merge TCR VDJ clonotype data (production: `scRepertoire`; demo: synthetic power-law)
- Rank by composite score: `mean_exhaustion × log1p(clone_size)`
- **Outputs:** `data/{cca,scvi}/ranked_clonotypes.csv`, `data/{cca,scvi}/query_with_clonotypes.{rds,h5ad}`

### Shared — NB04 · TCR Reactivity Selection (Python)
- Random Forest classifier (`class_weight="balanced"`, 5-fold CV)
- **Features:** exhaustion score, clone frequency, CDR3 length, TRBV family, rank score
- **Input:** `ranked_clonotypes.csv` from either method (set by pipeline script or papermill param)
- **Output:** `data/{cca,scvi}/selected_tumor_reactive_tcrs.csv`

### Shared — NB05 · PPV Validation (Python)
- Cross-references selected TCRs against Caushi et al. neoantigen reactivity data
- **Positive Predictive Value (PPV)** at k=20 selection; clinical target PPV ≥ 0.70
- **Output:** `data/{cca,scvi}/ppv_summary.csv`

---

## Exhaustion Gene Panel

18 genes curated from published TIL exhaustion literature (Wherry et al. 2007; Sade-Feldman et al. 2018; Oliveira et al. 2021). Stored in `data/exhaustion_gene_panel.txt` (shared by both pipelines):

```
PDCD1  HAVCR2  LAG3   TIGIT  CTLA4  TOX    NR4A1  ENTPD1  CXCL13
PRDM1  EOMES   TBX21  GZMB   PRF1   IFNG   TNF    CD38    VCAM1
```

---

## Quick Start

### 1 · Create the conda environment

```bash
conda env create -f environment.yml
conda activate scrnaseq

# Install R packages (first run only)
Rscript setup_r_env.R
```

> **Why both conda and renv?** conda manages Python + the R runtime; `renv` tracks exact R package versions in `renv.lock` the same way `requirements.txt` tracks Python packages.

**To reproduce on another machine:**
```bash
conda env create -f environment.yml
conda activate scrnaseq
Rscript -e "renv::restore()"
```

### 2 · Place raw data files

Download GEO files into `data/raw_downloads/`:

| GEO accession | File |
|--------------|------|
| GSE123813 | `GSE123813_bcc_scRNA_counts.txt.gz` |
| GSE123813 | `GSE123813_bcc_tcell_metadata.txt.gz` |
| GSE123813 | `GSE123813_bcc_tcr.txt.gz` |
| GSE176021 | `GSE176021_CD8_annotations.rds.gz` |

Paths are configured in `config/params.yaml` — no changes needed with the default folder layout.

### 3 · Run the CCA pipeline

```bash
conda run -n scrnaseq bash run_cca_pipeline.sh
```

Outputs land in `data/cca/`. Executed notebooks with full cell outputs are saved to `notebooks/executed/`.

### 4 · Run the scVI pipeline

The scVI pipeline requires the CCA NB00 outputs (`data/cca/reference.rds`, `data/cca/query.rds`) as its starting point:

```bash
# If you haven't run the CCA pipeline yet, run NB00 first:
conda run -n scrnaseq bash -c '
  papermill notebooks/cca/00_data_acquisition.ipynb /dev/null -k ir
'

# Then run the full scVI pipeline:
conda run -n scrnaseq bash run_scvi_pipeline.sh
```

Outputs land in `data/scvi/`.

### 5 · Run notebooks interactively

Open any notebook in JupyterLab or VS Code. All notebooks have a **parameters cell** (first code cell, tagged `parameters`) showing the default paths, which always point to the correct method subfolder. They can be run without the pipeline scripts.

### 6 · Run benchmarking

```bash
papermill benchmarking/benchmark_celltype.ipynb benchmarking/executed_benchmark_celltype.ipynb
papermill benchmarking/benchmark_tcr_ml.ipynb   benchmarking/executed_benchmark_tcr_ml.ipynb
```

### 7 · Run the leave-one-out patient benchmark

Holds out each patient in turn as the query, runs both CCA and scVI pipelines, and aggregates per-patient metrics:

```bash
conda run -n scrnaseq bash run_loo_benchmark.sh
# or: make run-loo
```

Outputs land in `data/loo/fold_<patient>_{cca,scvi}.csv` (gitignored — regeneratable).  
Summary figures and tables are written to `benchmarking/loo_*.png` and `benchmarking/loo_*.csv` (committed).  
Estimated runtime: ~4 hours for all 11 patients × 2 methods on an Apple Silicon MacBook.

To benchmark a subset of patients only:

```bash
conda run -n scrnaseq bash run_loo_benchmark.sh su001 su009
```

### Makefile shortcuts

All common operations have a `make` target:

```bash
make setup            # create conda env + install R packages (first run)
make restore          # recreate exact env from environment.yml + renv.lock
make run-cca          # conda run -n scrnaseq bash run_cca_pipeline.sh
make run-scvi         # conda run -n scrnaseq bash run_scvi_pipeline.sh
make run-loo          # conda run -n scrnaseq bash run_loo_benchmark.sh
make snapshot         # update environment.yml and renv.lock after adding packages
make docker-build     # build the Docker image
make docker-lab       # start JupyterLab in Docker (localhost:8888)
make docker-run-cca   # run CCA pipeline in Docker
make docker-run-scvi  # run scVI pipeline in Docker
```

### Live resource monitoring

During long training steps (scVI NB01 in particular), open a second terminal and run:

```bash
conda run -n scrnaseq python3 monitor_live.py
```

This polls CPU, RAM, and MPS/GPU memory every 3 seconds and displays a live dashboard. The pipeline scripts also write a CSV memory trace to `data/memory_monitor.log`.

---

## Docker (Alternative to Conda)

A `Dockerfile` is provided as a fully self-contained alternative to the conda + renv setup. The image packages Python 3.11, R 4.4, all Python libraries, and all R packages into a single reproducible Linux container.

> **When to prefer Docker over conda:**
> - Sharing the environment with collaborators who have no conda/R experience
> - Running on a cloud instance or HPC node without conda
> - CI/CD pipelines where you want environment parity

### Build the image

```bash
docker build -t scrna-projection .
# or: make docker-build
```

First build takes 10–20 minutes (Seurat + scvi-tools are large). Subsequent builds are fully cached unless `environment.yml` or `setup_r_env.R` change.

### Interactive JupyterLab

```bash
docker run --rm -p 8888:8888 \
    -v "$(pwd)/data":/project/data \
    -v "$(pwd)/models":/project/models \
    scrna-projection
# or: make docker-lab
```

Then open [http://localhost:8888](http://localhost:8888) in your browser. `data/` and `models/` are volume-mounted so all pipeline outputs persist on your host machine after the container exits.

### Run a pipeline non-interactively

```bash
# CCA pipeline
docker run --rm \
    -v "$(pwd)/data":/project/data \
    -v "$(pwd)/models":/project/models \
    scrna-projection bash run_cca_pipeline.sh
# or: make docker-run-cca

# scVI pipeline (requires data/cca/ outputs from CCA NB00 first)
docker run --rm \
    -v "$(pwd)/data":/project/data \
    -v "$(pwd)/models":/project/models \
    scrna-projection bash run_scvi_pipeline.sh
# or: make docker-run-scvi
```

### GPU / accelerator notes

| Platform | What to do |
|----------|-----------|
| macOS Apple Silicon (host) | Use conda — MPS is not available inside Docker |
| Linux + NVIDIA GPU | Replace the base image with a CUDA-enabled one and set `accelerator="cuda"` in `notebooks/scvi/01_train_reference.ipynb` |
| CPU-only (any platform) | No changes needed — scvi-tools falls back to CPU automatically |

### Supply raw data files inside the container

The `data/raw_downloads/` directory is volume-mounted together with `data/`, so place GEO files there before running:

```bash
mkdir -p data/raw_downloads
# copy GSE123813 and GSE176021 files into data/raw_downloads/
docker run --rm \
    -v "$(pwd)/data":/project/data \
    -v "$(pwd)/models":/project/models \
    scrna-projection bash run_cca_pipeline.sh
```

---

## Key Parameters (`config/params.yaml`)

Used by `run_cca_pipeline.sh`. The scVI pipeline reads the same file for raw input paths and QC thresholds.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `raw_counts_path` | `data/raw_downloads/GSE123813_bcc_scRNA_counts.txt.gz` | Raw count matrix |
| `reference_path` | `data/cca/reference.rds` | CCA reference Seurat object |
| `query_path` | `data/cca/query.rds` | CCA query Seurat object |
| `tcr_path` | `data/cca/yost2019_bcc_tcr.rds` | TCR clonotype table |
| `min_genes` | 200 | Minimum genes per cell (QC) |
| `max_genes` | 4000 | Maximum genes per cell (doublet proxy) |
| `max_pct_mito` | 25 | Maximum % mitochondrial reads |
| `n_dims` | 30 | PCA dimensions for CCA / UMAP |
| `min_cells_per_clonotype` | 2 | Minimum cells to include a clonotype |
| `n_exhausted_clonotypes` | 20 | Top-N exhausted clonotypes to carry forward |

---

## Adding a Third Method (e.g. scGPT)

The project is structured to make adding a new projection method straightforward:

1. Create `notebooks/scgpt/` with method-specific notebooks `00`–`03`
2. Create `data/scgpt/` for its outputs
3. Copy `run_scvi_pipeline.sh` → `run_scgpt_pipeline.sh` and update notebook paths and output paths to `data/scgpt/`
4. NB04 and NB05 are already shared — just pass `-p clonotype_table_path data/scgpt/ranked_clonotypes.csv` etc. via the pipeline script

---

## Technical Stack

| Component | Method | Language | Technology |
|-----------|--------|----------|------------|
| Data acquisition (NB00) | CCA | R | Seurat, GEOquery |
| Preprocessing (NB01) | CCA | R | [Seurat](https://satijalab.org/seurat/) |
| Reference projection (NB02) | CCA | R | Seurat `FindTransferAnchors` + `MapQuery` |
| Clonotype exhaustion (NB03) | CCA | R | Seurat `AddModuleScore`, scRepertoire |
| Format conversion (NB00) | scVI | Python | rpy2, anndata |
| Model training (NB01) | scVI | Python | [scvi-tools](https://scvi-tools.org) SCANVI |
| Query projection (NB02) | scVI | Python | scvi-tools |
| Clonotype exhaustion (NB03) | scVI | Python | scanpy `sc.tl.score_genes` |
| TCR selection (NB04) | Shared | Python | [scikit-learn](https://scikit-learn.org) Random Forest |
| PPV validation (NB05) | Shared | Python | pandas / numpy |
| Benchmarking | — | Python | scvi-tools, scGPT (optional) |
| Orchestration | Both | bash | [Papermill](https://papermill.readthedocs.io) |

---

## Notes on Demo vs. Production Data

| Aspect | Demo (this repo) | Production |
|--------|-----------------|------------|
| Input data | Yost 2019 BCC (GEO open access) | TIL scRNA-seq from patient biopsies |
| TCR data | Synthetic clonotypes (power-law) | 10x VDJ-seq `filtered_contig_annotations.csv` |
| Exhaustion scores | Simulated when panel genes absent from HVGs | Real expression from TIL samples |
| Reactivity labels | Caushi 2021 CDR3b overlap | Co-culture assays (IFN-γ ELISPOT / ICS) |
| Reference atlas | Yost 2019 70% split | Pan-cancer TIL atlas (Zheng et al. / internal) |

---

## Citation / Acknowledgements

Pipeline structure inspired by:
- Zheng et al. (2021) *Pan-cancer single-cell landscape of tumor-infiltrating T cells* — Science
- Sade-Feldman et al. (2018) *Defining T cell states associated with response to checkpoint immunotherapy* — Cell
- Oliveira et al. (2021) *Phenotype, specificity and avidity of antitumour CD8+ T cells* — Nature
- Caushi et al. (2021) *Transcriptional programs of neoantigen-specific TIL in lung cancer* — Nature
