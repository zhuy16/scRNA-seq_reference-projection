# Configuration Reference

## `config/params.yaml` — key parameters

Used by `run_cca_pipeline.sh` and `run_scvi_pipeline.sh`. All notebook defaults point here via Papermill injection.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `raw_counts_path` | `data/raw_downloads/GSE123813_bcc_scRNA_counts.txt.gz` | Raw count matrix |
| `reference_path` | `data/cca/reference.rds` | CCA reference Seurat object |
| `query_path` | `data/cca/query.rds` | CCA query Seurat object |
| `tcr_path` | `data/cca/yost2019_bcc_tcr.rds` | TCR clonotype table |
| `min_genes` | 200 | Minimum genes per cell (QC lower bound) |
| `max_genes` | 4000 | Maximum genes per cell (doublet proxy) |
| `max_pct_mito` | 25 | Maximum % mitochondrial reads |
| `n_dims` | 30 | PCA dimensions for CCA / UMAP |
| `min_cells_per_clonotype` | 2 | Minimum cells to include a clonotype |
| `n_exhausted_clonotypes` | 20 | Top-N exhausted clonotypes to carry forward |

---

## Demo vs. Production

| Aspect | Demo (this repo) | Production |
|--------|-----------------|------------|
| Input data | Yost 2019 BCC (GEO open access) | TIL scRNA-seq from patient biopsies |
| TCR data | Synthetic clonotypes (power-law) | 10x VDJ-seq `filtered_contig_annotations.csv` |
| Exhaustion scores | Simulated when panel genes absent from HVGs | Real expression from patient TIL samples |
| Reactivity labels | None — NB04/NB05 are placeholders | Co-culture assays (IFN-γ ELISPOT / ICS) |
| Reference atlas | Yost 2019 70% split | Pan-cancer TIL atlas (Zheng et al. / internal) |

To adapt to production data, update the paths in `config/params.yaml`. No notebook code changes are needed for standard TIL scRNA-seq data in 10x Genomics format.
