# =============================================================================
# Dockerfile — scRNA-seq Reference Projection Pipeline
# =============================================================================
# Builds a self-contained Linux image with Python 3.11, R 4.4, and all
# dependencies required to run both the CCA and scVI pipelines.
#
# Build:
#   docker build -t scrna-projection .
#
# Interactive JupyterLab:
#   docker run --rm -p 8888:8888 \
#     -v "$(pwd)/data":/project/data \
#     -v "$(pwd)/models":/project/models \
#     scrna-projection
#
# Run a pipeline (non-interactive):
#   docker run --rm \
#     -v "$(pwd)/data":/project/data \
#     -v "$(pwd)/models":/project/models \
#     scrna-projection bash run_cca_pipeline.sh
#
# Notes:
#   - MPS (Apple Silicon GPU) is not available inside Docker.
#     The scVI pipeline automatically falls back to CPU.
#   - For CUDA GPU support, replace the base image with a CUDA-enabled one
#     (e.g. condaforge/miniforge3 on top of nvidia/cuda) and set
#     accelerator="cuda" in notebooks/scvi/01_train_reference.ipynb.
# =============================================================================

FROM condaforge/miniforge3:latest

LABEL maintainer="scRNA-seq Reference Projection Pipeline"
LABEL description="CCA (Seurat) + scVI/scANVI dual-method TIL reference projection"

# ── System dependencies required by R packages ───────────────────────────────
# libgsl-dev      : used by Seurat / scattermore
# libhdf5-dev     : HDF5 for anndata / h5seurat
# libcurl4-*      : R httr / GEOquery downloads
# libxml2-dev     : XML / GEOquery
# libfont*/harfbuzz/fribidi : ggplot2 text rendering
# libfreetype/png/tiff/jpeg : image output (ggsave etc.)
RUN apt-get update && apt-get install -y --no-install-recommends \
        libgsl-dev \
        libhdf5-dev \
        libcurl4-openssl-dev \
        libssl-dev \
        libxml2-dev \
        libfontconfig1-dev \
        libharfbuzz-dev \
        libfribidi-dev \
        libfreetype6-dev \
        libpng-dev \
        libtiff5-dev \
        libjpeg-dev \
        build-essential \
        git \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /project

# ── Python + R conda environment ─────────────────────────────────────────────
# Copy only the spec first so this layer is re-used when only source code changes.
COPY environment.yml .
RUN mamba env create -f environment.yml \
    && mamba clean --all --yes

# All subsequent RUN steps execute inside the conda env
SHELL ["conda", "run", "--no-capture-output", "-n", "scrnaseq", "/bin/bash", "-c"]

# ── R packages (Seurat ecosystem) ────────────────────────────────────────────
# setup_r_env.R installs CRAN + Bioc packages and registers the ir kernel.
COPY setup_r_env.R .
RUN Rscript setup_r_env.R

# ── Project source code ───────────────────────────────────────────────────────
# .dockerignore excludes data/, models/, notebooks/executed/, renv/library/ etc.
COPY . .

# Ensure scripts are executable
RUN chmod +x run_cca_pipeline.sh run_scvi_pipeline.sh

# ── Runtime configuration ─────────────────────────────────────────────────────
# Disable MPS (not available in Linux containers).
# scvi-tools will fall back to CPU automatically when accelerator="mps" is set
# but MPS is unavailable; setting this variable makes the fallback explicit.
ENV PYTORCH_ENABLE_MPS_FALLBACK=1
# Prevent OMP conflicts in Linux (no DYLD_INSERT_LIBRARIES needed here)
ENV KMP_DUPLICATE_LIB_OK=TRUE
ENV OMP_NUM_THREADS=1
ENV OPENBLAS_NUM_THREADS=1

# JupyterLab listens on 8888
EXPOSE 8888

# Default entrypoint: run commands inside the conda env
ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "scrnaseq"]

# Default command: interactive JupyterLab
# Override with e.g. `docker run ... scrna-projection bash run_cca_pipeline.sh`
CMD ["jupyter", "lab", \
     "--ip=0.0.0.0", \
     "--port=8888", \
     "--no-browser", \
     "--allow-root", \
     "--NotebookApp.token=''", \
     "--NotebookApp.password=''"]
