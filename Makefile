IMAGE_NAME ?= scrna-projection
IMAGE_TAG  ?= latest

.PHONY: snapshot setup restore run-cca run-scvi run-loo \
        docker-build docker-lab docker-run-cca docker-run-scvi

# ── Snapshot both Python and R environments ───────────────────────────────────
# Run this after installing any new package (R or Python) to keep lockfiles current.
#   make snapshot
snapshot:
	@echo "Snapshotting conda/pip environment → environment.yml"
	conda env export --no-builds | grep -v "^prefix" > environment.yml
	@echo "Snapshotting R packages → renv.lock"
	Rscript -e "renv::snapshot(prompt = FALSE)"
	@echo "Done. Commit environment.yml and renv.lock."

# ── First-time environment setup ──────────────────────────────────────────────
# Creates the conda env (Python + R) and installs all R packages via renv.
#   make setup
setup:
	conda env create -f environment.yml
	conda run -n scrnaseq Rscript setup_r_env.R

# ── Restore exact environment on a new machine ────────────────────────────────
# Re-creates the conda env and restores pinned R packages from renv.lock.
#   make restore
restore:
	conda env create -f environment.yml
	conda run -n scrnaseq Rscript -e "renv::restore(prompt = FALSE)"

# ── Run the CCA / Seurat pipeline (NB00–03 R, NB04–05 Python) ────────────────
#   make run-cca
run-cca:
	conda run -n scrnaseq bash run_cca_pipeline.sh

# ── Run the scVI / scANVI pipeline (NB00–03 Python, NB04–05 Python) ──────────
# Requires data/cca/reference.rds and data/cca/query.rds from run-cca NB00.
#   make run-scvi
run-scvi:
	conda run -n scrnaseq bash run_scvi_pipeline.sh

# ── Run the leave-one-out patient benchmark (CCA + scVI) ─────────────────────
# Holds out each patient in turn, runs both pipelines, aggregates metrics.
# Outputs: data/loo/ (gitignored) + benchmarking/loo_*.{png,csv} (committed).
# Estimated runtime: ~4 h for all 11 patients × 2 methods.
#   make run-loo
run-loo:
	conda run -n scrnaseq bash run_loo_benchmark.sh

# ── Docker ────────────────────────────────────────────────────────────────────

# Build the Docker image (~10–20 min first time; cached on subsequent builds).
#   make docker-build
docker-build:
	docker build -t $(IMAGE_NAME):$(IMAGE_TAG) .

# Start an interactive JupyterLab session inside Docker.
# data/ and models/ are mounted as volumes so outputs persist on the host.
#   make docker-lab
docker-lab:
	docker run --rm -p 8888:8888 \
		-v "$(CURDIR)/data":/project/data \
		-v "$(CURDIR)/models":/project/models \
		$(IMAGE_NAME):$(IMAGE_TAG)

# Run the CCA pipeline non-interactively inside Docker.
#   make docker-run-cca
docker-run-cca:
	docker run --rm \
		-v "$(CURDIR)/data":/project/data \
		-v "$(CURDIR)/models":/project/models \
		$(IMAGE_NAME):$(IMAGE_TAG) bash run_cca_pipeline.sh

# Run the scVI pipeline non-interactively inside Docker.
# Requires data/cca/ outputs from docker-run-cca (or run-cca) first.
#   make docker-run-scvi
docker-run-scvi:
	docker run --rm \
		-v "$(CURDIR)/data":/project/data \
		-v "$(CURDIR)/models":/project/models \
		$(IMAGE_NAME):$(IMAGE_TAG) bash run_scvi_pipeline.sh
