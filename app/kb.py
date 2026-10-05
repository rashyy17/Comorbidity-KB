"""Data access for the Streamlit app.

The KB is small (about 20k edges), so each table is pulled once from Neo4j (or the SQLite mirror if Neo4j
is unreachable) and cached. Hub, overlap and enrichment calculations then run in pandas, using the same
rules as src/s05_queries.py.
"""
import os
import re
import sqlite3
from pathlib import Path

import numpy as np
import pandas as pd
import streamlit as st
from scipy import stats

ROOT = Path(__file__).resolve().parent.parent
SQLITE_PATH = ROOT / "results" / "comorbidity_kb.sqlite"
RESULTS = ROOT / "results"



def _setting(name, default=None):
    """Read a setting from the environment first, then st.secrets (no credentials live in the repo)."""
    if os.environ.get(name):
        return os.environ[name]
    try:
        return st.secrets.get(name, default)
    except Exception:  # no secrets.toml present
        return default


NEO4J_URI = _setting("NEO4J_URI", "bolt://localhost:7687")
NEO4J_USER = _setting("NEO4J_USER", "neo4j")
NEO4J_PASSWORD = _setting("NEO4J_PASSWORD")
# "auto" (default): use Neo4j when configured and reachable, else the bundled SQLite. "sqlite": always demo mode.
KB_BACKEND = str(_setting("KB_BACKEND", "auto")).lower()
SQLITE_URI = SQLITE_PATH.as_uri() + "?mode=ro"

DISEASES = ["PCOS", "T2D", "NAFLD", "Obesity"]
GEO = {"PCOS": ("GSE34526", "granulosa cells"), "T2D": ("GSE16415", "visceral adipose"),
       "NAFLD": ("GSE48452", "liver"), "Obesity": ("GSE55200", "subcutaneous adipose")}
DEG_P = 0.05
PATHWAY_FDR = 0.05


def bh_fdr(p):
    """Benjamini-Hochberg adjusted p-values (same as statsmodels multipletests(method='fdr_bh'))."""
    p = np.asarray(p, dtype=float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order] * n / np.arange(1, n + 1)
    adj = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.minimum(adj, 1.0)
    return out


# ---------------------------------------------------------------- connection

def neo4j_status():
    """'disabled' (demo mode forced or no credentials), 'connected', or 'unreachable'."""
    if KB_BACKEND == "sqlite" or not NEO4J_PASSWORD:
        return "disabled"
    return "connected" if get_driver() is not None else "unreachable"


@st.cache_resource(show_spinner=False)
def get_driver():
    """Return a verified Neo4j driver, or None if Neo4j is disabled or unreachable."""
    if KB_BACKEND == "sqlite" or not NEO4J_PASSWORD:
        return None
    try:
        from neo4j import GraphDatabase
        driver = GraphDatabase.driver(NEO4J_URI, auth=(NEO4J_USER, NEO4J_PASSWORD),
                                      connection_timeout=3)
        driver.verify_connectivity()
        with driver.session() as s:
            if s.run("MATCH (d:Disease) RETURN count(d) AS n").single()["n"] == 0:
                driver.close()
                return None  # reachable but empty: the SQLite copy is the better source
        return driver
    except Exception:
        return None


def backend():
    return "neo4j" if get_driver() is not None else "sqlite"


def _neo(query, **params):
    with get_driver().session(default_access_mode="READ") as s:
        return pd.DataFrame(s.run(query, **params).data())


def _sql(query, params=()):
    con = sqlite3.connect(SQLITE_URI, uri=True)
    try:
        return pd.read_sql_query(query, con, params=params)
    finally:
        con.close()


CYPHER = {
    "diseases": """MATCH (d:Disease) RETURN d.id AS id, d.full_name AS full_name, d.source_id AS source_id,
                   d.efo_xref AS efo_xref, d.ot_total_associations AS ot_total_associations,
                   d.kegg_universe_size AS kegg_universe_size ORDER BY id""",
    "genes": "MATCH (g:Gene) RETURN g.id AS gene_id, g.symbol AS symbol, g.entrez_id AS entrez_id, g.name AS name",
    "pathways": """MATCH (p:Pathway) RETURN p.kegg_id AS kegg_id, p.name AS pathway, p.category AS category,
                   p.size AS size""",
    "assoc": """MATCH (d:Disease)-[a:ASSOCIATED_WITH]->(g:Gene)
                RETURN d.id AS disease, g.id AS gene_id, g.symbol AS symbol, a.score AS score,
                       a.drug_only AS drug_only, a.literature_only AS literature_only,
                       a.evidence_types AS evidence_types""",
    "gene_pathway": """MATCH (g:Gene)-[:PARTICIPATES_IN]->(p:Pathway)
                       RETURN g.symbol AS symbol, p.kegg_id AS kegg_id""",
    "expr": """MATCH (g:Gene)-[e:EXPRESSED_IN]->(d:Disease)
               RETURN g.symbol AS symbol, d.id AS disease, e.log2FC AS log2FC, e.direction AS direction,
                      e.p_value AS p_value, e.adj_p AS adj_p, e.significant AS significant""",
}
SQL = {
    "diseases": """SELECT id, full_name, source_id, efo_xref, ot_total_associations, kegg_universe_size
                   FROM disease ORDER BY id""",
    "genes": "SELECT id AS gene_id, symbol, entrez_id, name FROM gene",
    "pathways": "SELECT kegg_id, name AS pathway, category, size FROM pathway",
    "assoc": """SELECT a.disease_id AS disease, a.gene_id, g.symbol, a.score, a.drug_only, a.literature_only,
                       a.evidence_types FROM associated_with a JOIN gene g ON g.id = a.gene_id""",
    "gene_pathway": """SELECT g.symbol, pi.pathway_id AS kegg_id
                       FROM participates_in pi JOIN gene g ON g.id = pi.gene_id""",
    "expr": """SELECT g.symbol, e.disease_id AS disease, e.log2FC, e.direction, e.p_value, e.adj_p, e.significant
               FROM expressed_in e JOIN gene g ON g.id = e.gene_id""",
}


@st.cache_data(show_spinner="Loading knowledge base...")
def load_tables(source: str) -> dict:
    """Pull every node/edge table once. `source` is part of the cache key."""
    if source == "neo4j":
        t = {k: _neo(q) for k, q in CYPHER.items()}
        t["assoc"]["evidence_types"] = t["assoc"].evidence_types.apply(
            lambda x: ";".join(x) if isinstance(x, list) else (x or ""))
    else:
        t = {k: _sql(q) for k, q in SQL.items()}
        t["assoc"]["evidence_types"] = t["assoc"].evidence_types.fillna("")
    t["assoc"]["drug_only"] = t["assoc"].drug_only.astype(bool)
    t["assoc"]["literature_only"] = t["assoc"].literature_only.astype(bool)
    t["expr"]["significant"] = t["expr"].significant.astype(bool)
    return t


def tables():
    return load_tables(backend())


# ---------------------------------------------------------------- derived views

def filtered_assoc(min_score, exclude_drug_only):
    a = tables()["assoc"]
    keep = a.score >= min_score
    if exclude_drug_only:
        keep &= ~a.drug_only
    return a[keep]


@st.cache_data(show_spinner="Computing hub genes...")
def hub_genes(min_score: float, exclude_drug_only: bool, source: str = "") -> pd.DataFrame:
    """One row per gene associated with >= 2 diseases, with per-disease score / log2FC / significance."""
    a = filtered_assoc(min_score, exclude_drug_only)
    expr = tables()["expr"]
    genes = tables()["genes"].drop_duplicates("symbol").set_index("symbol")

    g = a.groupby("symbol").agg(n_diseases=("disease", "nunique"),
                                diseases=("disease", lambda s: ", ".join(d for d in DISEASES if d in set(s))),
                                mean_score=("score", "mean"), max_score=("score", "max"),
                                drug_only_links=("drug_only", "sum"))
    g = g[g.n_diseases >= 2].copy()
    score = a.pivot_table(index="symbol", columns="disease", values="score")
    linked = a[["symbol", "disease"]].assign(linked=True)
    e = expr.merge(linked, on=["symbol", "disease"])  # expression only counts in associated diseases
    lfc = e.pivot_table(index="symbol", columns="disease", values="log2FC")
    sig = e[e.significant]
    for d in DISEASES:
        g[f"score_{d}"] = score.get(d, pd.Series(dtype=float)).reindex(g.index)
        g[f"log2FC_{d}"] = lfc.get(d, pd.Series(dtype=float)).reindex(g.index)
        g[f"p_{d}"] = e[e.disease == d].set_index("symbol").p_value.reindex(g.index)
        g[f"sig_{d}"] = g.index.isin(sig[sig.disease == d].symbol)

    sig_hub = sig[sig.symbol.isin(g.index)]
    dirs = sig_hub.groupby("symbol").direction.agg(["nunique", "count", "first"])
    g["n_deg"] = dirs["count"].reindex(g.index).fillna(0).astype(int)
    multi = dirs[dirs["count"] >= 2]
    cons = multi[multi["nunique"] == 1]
    g["expression_class"] = "insufficient data"
    g.loc[g.index.isin(multi.index), "expression_class"] = "discordant"
    g.loc[g.index.isin(cons.index), "expression_class"] = "consistent"
    g["direction"] = cons["first"].reindex(g.index)
    g["name"] = genes.name.reindex(g.index)
    g = g.reset_index()
    order = {"consistent": 0, "discordant": 1, "insufficient data": 2}
    g["_o"] = g.expression_class.map(order)
    g = g.sort_values(["n_diseases", "_o", "mean_score"], ascending=[False, True, False])
    return g.drop(columns="_o").reset_index(drop=True)


@st.cache_data(show_spinner="Running pathway over-representation tests...")
def pathway_enrichment(min_score: float, exclude_drug_only: bool, source: str = "") -> pd.DataFrame:
    """Hypergeometric test of each disease's genes in each non-overview KEGG pathway (BH FDR per disease)."""
    t = tables()
    a = filtered_assoc(min_score, exclude_drug_only)
    gp = t["gene_pathway"]
    pw = t["pathways"]
    universe = int(t["diseases"].kegg_universe_size.iloc[0])
    testable = pw[pw.category != "overview"]
    gp = gp[gp.kegg_id.isin(testable.kegg_id)]
    in_kegg = set(t["gene_pathway"].symbol)

    rows = []
    for d in DISEASES:
        dg = set(a[a.disease == d].symbol) & in_kegg
        k = gp[gp.symbol.isin(dg)].groupby("kegg_id").symbol.nunique()
        df = testable[["kegg_id", "pathway", "category", "size"]].copy()
        df["disease"] = d
        df["k"] = df.kegg_id.map(k).fillna(0).astype(int)
        df["n_disease_genes"] = len(dg)
        df["expected"] = df.n_disease_genes * df["size"] / universe
        df["p_value"] = np.where(df.k > 0, stats.hypergeom.sf(df.k - 1, universe, df["size"], len(dg)), 1.0)
        df["fdr"] = bh_fdr(df.p_value)
        rows.append(df)
    out = pd.concat(rows, ignore_index=True)
    out["fold_enrichment"] = np.where(out.expected > 0, out.k / out.expected, np.nan)
    out["enriched"] = (out.fdr < PATHWAY_FDR) & (out.k > 0)
    return out


@st.cache_data(show_spinner=False)
def pathway_summary(min_score: float, exclude_drug_only: bool, source: str = "") -> pd.DataFrame:
    enr = pathway_enrichment(min_score, exclude_drug_only, source)
    wide = enr.pivot_table(index=["kegg_id", "pathway", "category", "size"], columns="disease",
                           values=["k", "fdr"], aggfunc="first")
    wide.columns = [f"{a}_{b}" for a, b in wide.columns]
    wide = wide.reset_index()
    sig = enr[enr.enriched].groupby("kegg_id").disease.apply(
        lambda s: ", ".join(d for d in DISEASES if d in set(s)))
    wide["enriched_in"] = wide.kegg_id.map(sig).fillna("")
    wide["n_enriched"] = wide.enriched_in.apply(lambda s: len(s.split(", ")) if s else 0)
    wide["min_fdr"] = wide[[f"fdr_{d}" for d in DISEASES]].min(axis=1)
    return wide.sort_values(["n_enriched", "min_fdr"], ascending=[False, True]).reset_index(drop=True)


def disease_gene_sets(min_score, exclude_drug_only):
    a = filtered_assoc(min_score, exclude_drug_only)
    return {d: set(a[a.disease == d].symbol) for d in DISEASES}


# ---------------------------------------------------------------- query console

WRITE_KEYWORDS = ["CREATE", "MERGE", "DELETE", "DETACH", "SET", "REMOVE", "DROP", "FOREACH", "LOAD"]
SQL_WRITE_KEYWORDS = ["INSERT", "UPDATE", "DELETE", "DROP", "CREATE", "ALTER", "REPLACE", "ATTACH", "DETACH",
                      "PRAGMA", "VACUUM", "REINDEX"]


def blocked_keywords(query: str, keywords=WRITE_KEYWORDS) -> list:
    """Write keywords present in the query, ignoring string literals and comments."""
    stripped = re.sub(r"(//|--)[^\n]*", " ", query)
    stripped = re.sub(r"'(?:[^'\\]|\\.)*'|\"(?:[^\"\\]|\\.)*\"|`[^`]*`", " ", stripped)
    return [k for k in keywords if re.search(rf"\b{k}\b", stripped, flags=re.IGNORECASE)]


def _plain(v):
    """Convert neo4j graph objects to JSON-friendly values for display."""
    from neo4j.graph import Node, Path as GPath, Relationship
    if isinstance(v, Node):
        return {"labels": sorted(v.labels), **dict(v)}
    if isinstance(v, Relationship):
        return {"type": v.type, **dict(v)}
    if isinstance(v, GPath):
        return [_plain(n) for n in v.nodes]
    if isinstance(v, list):
        return [_plain(x) for x in v]
    return v


@st.cache_data(show_spinner=False, ttl=300)
def run_readonly_cypher(query: str, max_rows: int = 1000):
    """Run a query inside a READ transaction (Neo4j also rejects writes there). Returns (df, truncated)."""
    def work(tx):
        res = tx.run(query)
        rows, truncated = [], False
        for i, rec in enumerate(res):
            if i >= max_rows:
                truncated = True
                break
            rows.append({k: _plain(v) for k, v in rec.items()})
        return rows, truncated

    with get_driver().session(default_access_mode="READ") as s:
        rows, truncated = s.execute_read(work)
    df = pd.DataFrame(rows)
    for c in df.columns:  # dicts/lists -> strings so st.dataframe renders them
        if df[c].map(lambda x: isinstance(x, (dict, list))).any():
            df[c] = df[c].astype(str)
    return df, truncated


@st.cache_data(show_spinner=False, ttl=300)
def run_readonly_sql(query: str, max_rows: int = 1000):
    """Run SQL on a read-only connection to the bundled SQLite copy. Returns (df, truncated)."""
    con = sqlite3.connect(SQLITE_URI, uri=True)
    try:
        cur = con.execute(query)
        cols = [c[0] for c in cur.description or []]
        rows = cur.fetchmany(max_rows + 1)
    finally:
        con.close()
    return pd.DataFrame(rows[:max_rows], columns=cols), len(rows) > max_rows
