#!/usr/bin/env bash
# Launch the Streamlit frontend. Neo4j settings come from NEO4J_URI / NEO4J_USER / NEO4J_PASSWORD
# or .streamlit/secrets.toml (see secrets.toml.example); without them it runs in demo mode on the bundled SQLite.
set -euo pipefail
cd "$(dirname "$0")"
PY="${PYTHON:-.venv/bin/python}"
if [ ! -f results/comorbidity_kb.sqlite ]; then
  echo "Knowledge base not built yet - running ./run_all.sh first"; ./run_all.sh
fi
exec "$PY" -m streamlit run app/app.py --server.port "${PORT:-8501}" "$@"
