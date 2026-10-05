#!/usr/bin/env bash
# Run the full Comorbidity-Cluster KB pipeline end to end.
# Neo4j connection comes from NEO4J_URI / NEO4J_USER / NEO4J_PASSWORD
# (URI/user default to bolt://localhost:7687 and neo4j; the password has no default). If Neo4j is not configured or unreachable,
# step 4/5 fall back to SQLite + NetworkX automatically.
set -euo pipefail
cd "$(dirname "$0")"
PY="${PYTHON:-.venv/bin/python}"
export TMPDIR="$PWD/data/tmp"; mkdir -p "$TMPDIR"
export PYTHONWARNINGS="ignore"

"$PY" src/s01_opentargets.py   # disease IDs + disease-gene associations (Open Targets)
"$PY" src/s02_kegg.py          # KEGG gene->pathway links
"$PY" src/s03_geo_de.py        # GEO download + differential expression
"$PY" src/s04_load.py          # load Neo4j (+ SQLite mirror, GraphML)
"$PY" src/s05_queries.py       # hub queries, Cypher/SQL cross-check, sensitivity
"$PY" src/s06_validation.py    # INSR/IRS1/IGF1 + statistical validation
"$PY" src/s07_figures.py       # figures
echo "Done. See results/."
