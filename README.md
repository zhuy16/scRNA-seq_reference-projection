# scRNA-seq Reference Projection Pipeline

A reproducible, [Papermill](https://papermill.readthedocs.io)-parameterised pipeline for **tumor-infiltrating lymphocyte (TIL)** analysis with a primary focus on **cell type annotation**: QC → reference projection → clonotype × exhaustion integration.

NB04 and NB05 are currently under development as **placeholder downstream prediction modules** (TCR reactivity selection and PPV validation) and are not the main focus of this repository.

Two projection methods run in parallel and can be compared head-to-head:

| Method | Language | Script |
|--------|----------|--------|
| **CCA** — Seurat Canonical Correlation Analysis | R | `run_cca_pipeline.sh` |
| **scVI/scANVI** — variational autoencoder | Python | `run_scvi_pipeline.sh` |

---

## Demo Dataset

**Yost et al. 2019 BCC** ([GSE123813](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE123813), GEO open access) — paired scRNA-seq + TCR-seq from 11 basal cell carcinoma patients (pre/post anti-PD1). Used as a fully reproducible stand-in for clinical TIL data; replace GEO paths in `config/params.yaml` to run on in-house data.

> Yost KE et al. *Clonal replacement of tumor-specific T cells following PD-1 blockade.* Nat Med. 2019. https://doi.org/10.1038/s41591-019-0522-3

---

## Pipeline Overview

```
Raw GEO data
    │
    ├─[CCA]   NB00 acquire → NB01 QC/preprocess → NB02 CCA projection ─┐
    │                                                                     ├─ NB03 TCR × exhaustion scoring
    └─[scVI]  NB00 convert → NB01 train SCANVI  → NB02 project query  ─┘
                                                                          │
                                                          NB04 TCR selection (placeholder)
                                                          NB05 PPV validation (placeholder)
```

NB03 has parallel R (Seurat `AddModuleScore`) and Python (scanpy `sc.tl.score_genes`) implementations — both produce identical output schemas. NB04–05 are placeholder implementations; the optimal ML approach is under benchmarking pending labelled training data. See [docs/pipeline-reference.md](docs/pipeline-reference.md) for per-step details.

---

## Repository Structure

```
├── run_cca_pipeline.sh / run_scvi_pipeline.sh / run_loo_benchmark.sh
├── Makefile                   # make setup / run-cca / run-scvi / run-loo / …
├── Dockerfile                 # self-contained alternative to conda (→ docs/docker.md)
├── config/params.yaml         # single source of truth for all paths + QC params
│
├── notebooks/
│   ├── cca/                   # NB00–03 in R/Seurat
│   ├── scvi/                  # NB00–03 in Python/scvi-tools
│   ├── benchmarking/          # benchmark notebooks + LOO summary outputs
│   ├── 04_tcr_reactivity_selection.ipynb
│   ├── 05_ppv_validation.ipynb
│   └── examples/              # selected executed notebooks with outputs (see below)
│
├── data/exhaustion_gene_panel.txt  # 18-gene panel (→ docs/exhaustion-gene-panel.md)
├── environment.yml / renv.lock     # exact Python + R environments
└── docs/                           # detailed references (see below)
```

---

## Benchmarking: CCA vs scVI/scANVI

Leave-one-out benchmark across 11 patients. Each patient held out as query; remaining 10 as reference. Full results: `notebooks/benchmarking/loo_summary.csv`, plots: `notebooks/benchmarking/loo_*.png`.

| Metric | CCA | scANVI |
|--------|-----|--------|
| Overall accuracy | 86.2 ± 4.8% | 85.9 ± 5.0% |
| Macro F1 | 73.7 ± 7.5% | 73.7 ± 4.3% |
| CD8\_ex recall | 56.3 ± 42.9% | **88.3 ± 13.5%** |
| CD8\_ex F1 | 51.2 ± 42.5% | **68.3 ± 20.9%** |

Overall accuracy is equivalent. The key difference is **exhausted CD8 recovery**: scANVI achieves higher recall (88% vs 56%) and F1 (68% vs 51%) with much lower patient-to-patient variance — the critical subtype for TIL selection. CCA is a fast, GPU-free baseline suited for development and interpretability.

---

## Example Notebooks

Selected executed notebooks with full outputs are committed to illustrate code style, logic flow, and visualisations without running the full pipeline:

| Notebook | What to look for |
|----------|-----------------|
| [explore_yost2019_bcc](notebooks/examples/explore_yost2019_bcc.ipynb) | Dataset structure, cell-type composition, TCR overlap overview |
| [benchmark_celltype](notebooks/examples/benchmark_celltype.ipynb) | CCA vs scANVI accuracy, per-patient LOO results, confusion matrices |
| [loo_su001_cca/02_reference_projection](notebooks/executed/loo_su001_cca/02_reference_projection.ipynb) | LOO reference projection behavior with held-out patient |
| [loo_su001_scvi/02_project_query](notebooks/executed/loo_su001_scvi/02_project_query.ipynb) | LOO SCANVI inference on held-out patient |

---

## Quick Start

```bash
# 1. Create environment
conda env create -f environment.yml && conda activate scrnaseq
Rscript setup_r_env.R          # first run only — installs R packages via renv

# 2. Download GEO data (GSE123813) into data/raw_downloads/

# 3. Run CCA pipeline (R)
conda run -n scrnaseq bash run_cca_pipeline.sh

# 4. Run scVI pipeline (Python) — requires CCA NB00 outputs
conda run -n scrnaseq bash run_scvi_pipeline.sh

# 5. LOO benchmark across all 11 patients (~4 h)
conda run -n scrnaseq bash run_loo_benchmark.sh
```

All steps have `make` shortcuts (`make setup`, `make run-cca`, `make run-scvi`, `make run-loo`). For Docker, see [docs/docker.md](docs/docker.md).

---

## Documentation

| File | Contents |
|------|----------|
| [docs/pipeline-reference.md](docs/pipeline-reference.md) | Data flow diagram, per-step logic, technical stack, adding new methods |
| [docs/configuration.md](docs/configuration.md) | `params.yaml` reference, demo vs production comparison |
| [docs/exhaustion-gene-panel.md](docs/exhaustion-gene-panel.md) | 18-gene panel, per-gene rationale, source references |
| [docs/docker.md](docs/docker.md) | Docker build, JupyterLab, GPU notes |

---

## Citation / Acknowledgements

- Yost KE et al. (2019) *Clonal replacement of tumor-specific T cells following PD-1 blockade* — Nature Medicine
- Zheng et al. (2021) *Pan-cancer single-cell landscape of tumor-infiltrating T cells* — Science
- Sade-Feldman et al. (2018) *Defining T cell states associated with response to checkpoint immunotherapy* — Cell
- Oliveira et al. (2021) *Phenotype, specificity and avidity of antitumour CD8+ T cells* — Nature
- Caushi et al. (2021) *Transcriptional programs of neoantigen-specific TIL in lung cancer* — Nature

