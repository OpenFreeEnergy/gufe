# syntax=docker/dockerfile:1
#
# Self-contained dev image for the gufe framejs visualization work.
# Builds from gufe's own environment.yml (native arm64, fast). Everything the
# `.view()` integration targets (LigandNetwork/SmallMoleculeComponent/
# ProteinComponent/LigandAtomMapping) lives in gufe, so this single light image
# covers all notebook-visualization development — no openfe needed.
#
# BuildKit cache mounts keep the conda + pip downloads across builds, so an
# interrupted build RESUMES instead of re-downloading.
FROM mambaorg/micromamba:latest

# git + certs: env files may include `pip: git+https://...` lines.
USER root
RUN apt-get update \
 && apt-get install -y --no-install-recommends git ca-certificates \
 && rm -rf /var/lib/apt/lists/*
USER $MAMBA_USER

# --- conda env solve (cached layer; downloads cached via the pkgs mount) -------
ARG CONDA_ENV_FILE=environment.yml
COPY --chown=$MAMBA_USER:$MAMBA_USER ${CONDA_ENV_FILE} /tmp/environment.yml
RUN --mount=type=cache,target=/opt/conda/pkgs,uid=57439,gid=57439,sharing=locked \
    micromamba install -y -n base -f /tmp/environment.yml

# activate base for the remaining RUN/CMD steps
ARG MAMBA_DOCKERFILE_ACTIVATE=1

# Notebook tooling, the widget (PyPI), and editable-install build backends.
# Backends must be present for `pip install -e --no-build-isolation`:
#   setuptools/setuptools-scm → gufe;  hatchling/editables → metaframe-widget
#   (only needed for `just dev-widget`, harmless otherwise).
RUN --mount=type=cache,target=/home/mambauser/.cache/pip,uid=57439,gid=57439,sharing=locked \
    python -m pip install \
      jupyterlab marimo ipykernel \
      "metaframe-widget==0.3.11" \
      "setuptools>=77.0.3" "setuptools-scm>=8" wheel \
      hatchling editables

COPY --chown=$MAMBA_USER:$MAMBA_USER visualization-demo/dev-start.sh /usr/local/bin/dev-start.sh
RUN chmod +x /usr/local/bin/dev-start.sh

WORKDIR /workspace
# Put the base env on PATH so `docker compose exec <svc> python|pytest ...`
# works even though exec bypasses micromamba's activation entrypoint.
ENV PATH=/opt/conda/bin:$PATH \
    PYTHONDONTWRITEBYTECODE=1 \
    SETUPTOOLS_SCM_PRETEND_VERSION=0.0.0+docker

EXPOSE 8888 2718
# ENTRYPOINT stays micromamba's _entrypoint.sh (env activation).
# The service command (dev-start.sh + server) is set per-service in compose.
