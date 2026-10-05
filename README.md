# Comorbidity-Cluster Knowledge Base (PCOS · T2D · NAFLD · Obesity)

A knowledge graph of the genes and pathways shared by four commonly co-occurring metabolic diseases.
It is built from **Open Targets** (disease→gene), **KEGG** (gene→pathway) and **GEO** (gene expression
in patients vs controls), stored in **Neo4j** and mirrored to SQLite. Implements
`Comorbidity_Cluster_KB_Writeup-2.docx`.

- Plain-language walkthrough → [EXPLAINED.md](EXPLAINED.md)
- Every judgment call → [ASSUMPTIONS.md](ASSUMPTIONS.md)
- Validation (INSR / IRS1 / IGF1 + statistical checks) → [VALIDATION.md](VALIDATION.md)

## Quick start

Requirements: Python 3.10+ (tested 3.12), internet access, and optionally Docker for Neo4j.

```bash
# 1. Python environment (requirements.txt = app only; the pipeline needs the extra build packages)
python3 -m venv .venv
.venv/bin/pip install -r requirements-pipeline.txt

# 2. Neo4j 5 on bolt://localhost:7687 (skip if you already run one, e.g. the neo4j-kb container).
#    Credentials come only from the environment; nothing is stored in the repo.
export NEO4J_PASSWORD=<your password>          # optional: NEO4J_URI, NEO4J_USER
docker compose up -d

# 3. Run everything (about 5-10 min, mostly the GEO downloads; GSE48452 and GSE55200 are 27-39 MB)
./run_all.sh
```

Without Neo4j the pipeline logs a warning and continues on SQLite + NetworkX. All queries and results still
work (`results/backend.txt` records which backend was used).

Explore the graph in the Neo4j Browser at <http://localhost:7474>:

```cypher
:param min_score => 0.1
:param exclude_drug_only => false
// then paste any query from queries/hubs.cypher, e.g. the INSR neighbourhood:
MATCH (d:Disease)-[a:ASSOCIATED_WITH]->(g:Gene {symbol:'INSR'})-[:PARTICIPATES_IN]->(p:Pathway)
RETURN d, a, g, p LIMIT 50
```

Or use SQLite: `sqlite3 results/comorbidity_kb.sqlite` and run queries from `queries/hubs.sql`
(replace `:min_score` with `0.1` and `:exclude_drug_only` with `0`).

## Web app

An interactive Streamlit frontend in `app/` reads from the Neo4j database built above. When Neo4j is not
configured or not reachable, it runs in **demo mode** on the bundled read-only SQLite copy
(`results/comorbidity_kb.sqlite`, committed to the repo). The sidebar says which mode is active.

```bash
pip install -r requirements.txt   # app only
streamlit run app/app.py          # or ./run_app.sh
# open http://localhost:8501
```

To use the live Neo4j database, set `NEO4J_URI`, `NEO4J_USER` and `NEO4J_PASSWORD` as environment
variables, or copy `.streamlit/secrets.toml.example` to `.streamlit/secrets.toml` (git-ignored) and fill it
in. `KB_BACKEND=sqlite` forces demo mode.

| Page | What it shows |
|---|---|
| Home | Summary, headline numbers, schema diagram, data sources |
| Hub Explorer | Filterable hub-gene table (min diseases, min score, drug-only toggle, consistent-only, search) + CSV download |
| Gene Profile | One gene's associations, evidence types, log2FC/p per disease, pathways, interactive neighbourhood graph |
| Disease Overlap | UpSet-style intersection chart for 2-4 diseases, pairwise overlap, shared-gene table |
| Pathways | Over-representation statistics per disease, pathway-to-diseases view |
| Network | Interactive graph of top hubs, diseases and shared pathways (pyvis) |
| Expression | log2FC heatmap of top hubs across GSE34526 / GSE16415 / GSE48452 / GSE55200 |
| Validation & Limitations | INSR / IRS1 / IGF1 results and caveats, pulled from VALIDATION.md and results/ |
| Query Console | Read-only Cypher with example queries (Neo4j mode), or read-only SQL on the SQLite copy (demo mode) |

Checks: `.venv/bin/python app/tests/smoke_test.py` renders every page headlessly and fails on any
exception. With the app running, `.venv/bin/python app/tests/screenshots.py` opens each page in
headless Chrome, fails on error elements and saves `results/screenshots/*.png`.

## Deploy to Streamlit Community Cloud

The repo is ready to deploy as-is. The app runs entirely from the committed SQLite database, so no
Neo4j or secrets are needed.

1. Push the repo to GitHub (public, or private with Streamlit access).
2. At <https://share.streamlit.io> choose **Create app**, then select the repo and branch `main`, and set
   **Main file path** to `app/app.py`.
3. Under **Advanced settings**, pick **Python 3.12**. No secrets are required. To point a deployment at a
   reachable Neo4j instead, add `NEO4J_URI`, `NEO4J_USER` and `NEO4J_PASSWORD` under **Secrets**.

Streamlit installs `requirements.txt` (app dependencies only) and picks up `.streamlit/config.toml`.

## Pipeline

| Step | Script | Output |
|---|---|---|
| 1 | `src/s01_opentargets.py`: resolve disease IDs, fetch associations (score ≥ 0.1) with per-datatype evidence | `data/processed/diseases.csv`, `associations.csv` |
| 2 | `src/s02_kegg.py`: `rest.kegg.jp/link/pathway/hsa` + pathway names + symbol→Entrez | `gene_pathway.csv`, `pathways.csv` |
| 3 | `src/s03_geo_de.py` (+ `src/limma_py.py`): download GEO SOFT, verify groups, probe→gene, moderated t-test | `de_<disease>.csv`, `expression.csv`, `results/geo_datasets.json` |
| 4 | `src/s04_load.py`: load Neo4j, SQLite mirror, GraphML | Neo4j DB, `results/comorbidity_kb.sqlite`, `results/comorbidity_kb.graphml` |
| 5 | `src/s05_queries.py`: hub queries in Cypher **and** SQL, cross-checked; pathway enrichment; sensitivity | `results/hub_*.csv`, `sensitivity.csv`, `backend_crosscheck.txt` |
| 6 | `src/s06_validation.py`: INSR/IRS1/IGF1, insulin-resistance enrichment, permutation and background nulls | `results/validation*.{json,csv}` |
| 7 | `src/s07_figures.py`: figures | `results/figures/*.png` |

Shared settings (thresholds, Neo4j credentials) are in `src/config.py`.

## Schema (Neo4j)

```
(:Disease {id, name, full_name, source_id, efo_xref, ot_total_associations, kegg_universe_size})
(:Gene    {id (Ensembl), symbol, entrez_id, name})
(:Pathway {id, kegg_id, name, category, size})

(:Disease)-[:ASSOCIATED_WITH {score, evidence_types, drug_only, literature_only}]->(:Gene)
(:Gene)-[:PARTICIPATES_IN]->(:Pathway)
(:Gene)-[:EXPRESSED_IN {log2FC, direction, p_value, adj_p, significant}]->(:Disease)
```

| | count |
|---|---|
| Disease / Gene / Pathway nodes | 4 / 3,303 / 366 |
| ASSOCIATED_WITH / PARTICIPATES_IN / EXPRESSED_IN | 4,512 / 11,356 / 3,275 |

## Data

| Disease | Open Targets ID | genes (score ≥ 0.1) | GEO | tissue | case vs control |
|---|---|---|---|---|---|
| PCOS | MONDO_0008487 | 287 | GSE34526 (GPL570) | granulosa cells | 7 vs 3 |
| T2D | MONDO_0005148 | 2,244 | GSE16415 (GPL2986) | visceral adipose | 5 vs 5 |
| NAFLD | EFO_1001248 | 422 | GSE48452 (GPL11532) | liver | 32 (steatosis+NASH) vs 14 |
| Obesity | MONDO_0011122 | 1,559 | GSE55200 (GPL17692) | subcutaneous adipose | 16 (MHO+MUO) vs 7 |

## Key results (score ≥ 0.1)

- **855 hub genes** (associated with ≥ 2 diseases): 99 shared by all four, 255 by three or more.
  About half of the 4-disease hubs are drug-only artifacts (mostly metformin's complex I targets).
  With those edges removed, 48 remain.
- **41 expression-consistent hubs**, e.g. CD4 ↑ (PCOS, T2D, Obesity), MMP9 ↑, TGFB1 ↑, CCL2 ↑,
  PPARG ↓, ADIPOQ ↓, INSR ↓, RBP4 ↓, TXNIP ↓. 9 are discordant (e.g. GHR, CD36, FABP4).
- **143 hub pathways** (enriched in ≥ 2 diseases, FDR < 0.05), including Insulin resistance, AMPK,
  adipocytokine, FoxO, Type II diabetes mellitus and NAFLD.
- **Validation:** INSR, IRS1 and IGF1 are all hubs. Only INSR is expression-consistent, and not in the
  PCOS/T2D datasets. IGF1's hub status is fragile. See [VALIDATION.md](VALIDATION.md).

## Figures

| | |
|---|---|
| `results/figures/fig1_overlap_graph.png` | Disease–gene–pathway overlap graph (top robust hubs and shared pathways) |
| `results/figures/fig2_hub_gene_table.png` | Hub gene table: association dot-matrix, GEO direction, evidence |
| `results/figures/fig3_expression_heatmap.png` | log2FC heatmap of the top 40 hub genes across the four datasets |

## Repository layout

```
app/            Streamlit frontend (app.py, kb.py data layer, ui.py, views/, tests/)
src/            pipeline steps s01–s07, config.py, limma_py.py
queries/        hubs.cypher (Neo4j) and hubs.sql (SQLite), same queries in both
data/raw/       downloaded GEO/KEGG/Open Targets files (git-ignored, about 90 MB, re-created by the pipeline)
data/processed/ tidy CSVs fed into the graph
results/        query results, validation outputs, SQLite DB, GraphML, figures/, screenshots/
```

## Caveats

Data reflect the live APIs on 2026-10-05. The PCOS and T2D expression datasets are small, and the four
datasets come from different tissues. Read the limitations in [EXPLAINED.md](EXPLAINED.md#limitations)
before drawing biological conclusions.
