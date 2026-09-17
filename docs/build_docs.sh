#!/bin/bash
# Builds the mkdocs documentation site (source: docs/, config: ../mkdocs.yml)
# into docs/site/. Uses mkdocs from the "pybio" micromamba environment
# directly by path, so this works whether or not that environment is
# currently activated. Run from anywhere -- it cd's to the repo root itself.
#
# Usage:
#   docs/build_docs.sh          # build the static site into docs/site/
#   docs/build_docs.sh serve    # serve locally with live-reload for editing

set -e

MKDOCS=/home/gregor/micromamba/envs/pybio/bin/mkdocs
# this script lives in docs/, but mkdocs.yml lives one level up at the repo root
cd "$(dirname "$0")/.."

if [ "$1" == "serve" ]; then
    "$MKDOCS" serve
else
    "$MKDOCS" build --strict
    echo "Built docs/site/ -- open docs/site/index.html or run 'docs/build_docs.sh serve' to preview locally."
fi
