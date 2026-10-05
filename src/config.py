"""Shared configuration for the Comorbidity-Cluster KB pipeline."""
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "data" / "raw"
PROC = ROOT / "data" / "processed"
RESULTS = ROOT / "results"
FIGURES = RESULTS / "figures"
QUERIES = ROOT / "queries"
SQLITE_PATH = RESULTS / "comorbidity_kb.sqlite"

for p in (RAW / "geo", RAW / "kegg", RAW / "opentargets", PROC, RESULTS, FIGURES):
    p.mkdir(parents=True, exist_ok=True)

# GEOparse writes temp files to TMPDIR; keep them inside the project.
_tmp = ROOT / "data" / "tmp"
_tmp.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("TMPDIR", str(_tmp))

OPENTARGETS_URL = "https://api.platform.opentargets.org/api/v4/graphql"
KEGG_URL = "https://rest.kegg.jp"

# Disease short name -> search string used to resolve the Open Targets ontology ID.
DISEASES = {
    "PCOS": "polycystic ovary syndrome",
    "T2D": "type 2 diabetes mellitus",
    "NAFLD": "non-alcoholic fatty liver disease",
    "Obesity": "obesity",
}

# Expected resolutions (verified 2026-10-05); used as a guard against search drift.
EXPECTED_IDS = {
    "PCOS": "MONDO_0008487",
    "T2D": "MONDO_0005148",
    "NAFLD": "EFO_1001248",
    "Obesity": "MONDO_0011122",
}

# Minimum Open Targets overall association score for a Disease-[ASSOCIATED_WITH]->Gene edge.
MIN_ASSOC_SCORE = float(os.environ.get("MIN_ASSOC_SCORE", 0.1))
# A datatype (genetic, literature, clinical, ...) counts as supporting an edge if its score >= this.
EVIDENCE_MIN = 0.1
# Thresholds used for the sensitivity analysis of hub calls.
SENSITIVITY_SCORES = [0.1, 0.2, 0.3]

# Significance rule for a DEG (written onto EXPRESSED_IN.significant).
DEG_P = 0.05          # nominal moderated-t p-value
DEG_ABS_LFC = 0.0     # no fold-change floor (see ASSUMPTIONS.md)

# Pathway hub: disease gene set enriched in pathway (hypergeometric, BH-adjusted).
PATHWAY_FDR = 0.05

NEO4J_URI = os.environ.get("NEO4J_URI", "bolt://localhost:7687")
NEO4J_USER = os.environ.get("NEO4J_USER", "neo4j")
# No default password: set NEO4J_PASSWORD in the environment. Without it, steps 4-5 fall back to SQLite.
NEO4J_PASSWORD = os.environ.get("NEO4J_PASSWORD")

VALIDATION_GENES = ["INSR", "IRS1", "IGF1"]
