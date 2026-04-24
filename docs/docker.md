# Docker — Reproducible Container Environment

A `Dockerfile` is provided as a fully self-contained alternative to the conda + renv setup. The image packages Python 3.11, R 4.4, all Python libraries, and all R packages into a single reproducible Linux container.

**When to prefer Docker over conda:**
- Sharing with collaborators who have no conda/R experience
- Running on a cloud instance or HPC node without conda
- CI/CD pipelines where you want environment parity

## Build the image

```bash
docker build -t scrna-projection .
# or: make docker-build
```

First build takes 10–20 minutes (Seurat + scvi-tools are large). Subsequent builds are fully cached unless `environment.yml` or `setup_r_env.R` change.

## Interactive JupyterLab

```bash
docker run --rm -p 8888:8888 \
    -v "$(pwd)/data":/project/data \
    -v "$(pwd)/models":/project/models \
    scrna-projection
# or: make docker-lab
```

Then open http://localhost:8888 in your browser. `data/` and `models/` are volume-mounted so all pipeline outputs persist on your host machine after the container exits.

## Run a pipeline non-interactively

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

## GPU / accelerator notes

| Platform | What to do |
|----------|-----------|
| macOS Apple Silicon (host) | Use conda — MPS is not available inside Docker |
| Linux + NVIDIA GPU | Replace the base image with a CUDA-enabled one and set `accelerator="cuda"` in `notebooks/scvi/01_train_reference.ipynb` |
| CPU-only (any platform) | No changes needed — scvi-tools falls back to CPU automatically |

## Supplying raw data files

The `data/raw_downloads/` directory is volume-mounted together with `data/`, so place GEO files there before running:

```bash
mkdir -p data/raw_downloads
# copy GSE123813 files into data/raw_downloads/
docker run --rm \
    -v "$(pwd)/data":/project/data \
    -v "$(pwd)/models":/project/models \
    scrna-projection bash run_cca_pipeline.sh
```
