#!/usr/bin/env bash
# Container preamble: editable-install the mounted gufe source over the env, then
# exec the service command. Fast: build backends are pre-installed and we skip
# build isolation + dependency resolution (deps come from the conda env).
set -uo pipefail

export SETUPTOOLS_SCM_PRETEND_VERSION="${SETUPTOOLS_SCM_PRETEND_VERSION:-0.0.0+docker}"

editable() {
  local d="$1"
  if [ -d "$d" ] && [ -f "$d/pyproject.toml" ]; then
    echo "[dev-start] pip install -e $d"
    python -m pip install -e "$d" --no-build-isolation --no-deps -q \
      || echo "[dev-start] WARN: editable install failed for $d (continuing)"
  fi
}

# gufe (the viz classes live here). metaframe-widget only if a checkout is
# mounted for widget dev; otherwise the PyPI build baked into the image is used.
editable /src/gufe
editable /src/metaframe-widget

echo "[dev-start] launching: $*"
exec "$@"
