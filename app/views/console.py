import time

import streamlit as st

import kb
import ui

EXAMPLES = {
    "Hub genes shared by all four diseases (no drug-only links)": """MATCH (d:Disease)-[a:ASSOCIATED_WITH]->(g:Gene)
WHERE a.score >= 0.1 AND NOT a.drug_only
WITH g, collect(d.id) AS diseases, avg(a.score) AS mean_score
WHERE size(diseases) = 4
RETURN g.symbol AS gene, round(mean_score, 3) AS mean_score
ORDER BY mean_score DESC
LIMIT 25""",
    "Expression-consistent hub genes": """MATCH (d:Disease)-[a:ASSOCIATED_WITH]->(g:Gene)
WHERE a.score >= 0.1
WITH g, collect(d.id) AS diseases WHERE size(diseases) >= 2
MATCH (g)-[e:EXPRESSED_IN]->(d2:Disease)
WHERE d2.id IN diseases AND e.significant
WITH g, collect(d2.id) AS deg_diseases, collect(DISTINCT e.direction) AS dirs
WHERE size(deg_diseases) >= 2 AND size(dirs) = 1
RETURN g.symbol AS gene, deg_diseases, dirs[0] AS direction
ORDER BY size(deg_diseases) DESC, gene""",
    "Everything about INSR": """MATCH (d:Disease)-[a:ASSOCIATED_WITH]->(g:Gene {symbol: 'INSR'})
OPTIONAL MATCH (g)-[e:EXPRESSED_IN]->(d)
RETURN d.id AS disease, round(a.score, 3) AS score, a.evidence_types AS evidence,
       round(e.log2FC, 3) AS log2FC, e.direction AS direction, e.p_value AS p_value, e.significant AS significant
ORDER BY disease""",
    "Pathways with the most PCOS + T2D genes": """MATCH (p:Pathway)<-[:PARTICIPATES_IN]-(g:Gene)
WHERE p.category <> 'overview'
  AND EXISTS { (:Disease {id: 'PCOS'})-[:ASSOCIATED_WITH]->(g) }
  AND EXISTS { (:Disease {id: 'T2D'})-[:ASSOCIATED_WITH]->(g) }
RETURN p.kegg_id AS kegg_id, p.name AS pathway, count(g) AS shared_genes, collect(g.symbol)[..10] AS examples
ORDER BY shared_genes DESC
LIMIT 15""",
    "Pairwise gene overlap between diseases": """MATCH (d1:Disease)-[a1:ASSOCIATED_WITH]->(g:Gene)<-[a2:ASSOCIATED_WITH]-(d2:Disease)
WHERE d1.id < d2.id
RETURN d1.id AS disease_a, d2.id AS disease_b, count(g) AS shared_genes
ORDER BY shared_genes DESC""",
    "Drug-only links (the metformin effect)": """MATCH (d:Disease)-[a:ASSOCIATED_WITH]->(g:Gene)
WHERE a.drug_only
WITH g, collect(d.id) AS diseases
WHERE size(diseases) >= 3
RETURN g.symbol AS gene, diseases
ORDER BY size(diseases) DESC, gene
LIMIT 30""",
}


SQL_EXAMPLES = {
    "Hub genes shared by all four diseases (no drug-only links)": """SELECT g.symbol AS gene, ROUND(AVG(a.score), 3) AS mean_score
FROM associated_with a JOIN gene g ON g.id = a.gene_id
WHERE a.score >= 0.1 AND a.drug_only = 0
GROUP BY g.id
HAVING COUNT(DISTINCT a.disease_id) = 4
ORDER BY mean_score DESC
LIMIT 25""",
    "Expression-consistent hub genes": """WITH hubs AS (
  SELECT gene_id FROM associated_with WHERE score >= 0.1
  GROUP BY gene_id HAVING COUNT(*) >= 2
), degs AS (
  SELECT e.gene_id, e.disease_id, e.direction FROM expressed_in e
  JOIN associated_with a ON a.gene_id = e.gene_id AND a.disease_id = e.disease_id
  WHERE e.significant = 1
)
SELECT g.symbol AS gene, GROUP_CONCAT(d.disease_id) AS deg_diseases, MIN(d.direction) AS direction
FROM hubs h JOIN degs d ON d.gene_id = h.gene_id JOIN gene g ON g.id = h.gene_id
GROUP BY h.gene_id
HAVING COUNT(*) >= 2 AND COUNT(DISTINCT d.direction) = 1
ORDER BY COUNT(*) DESC, gene""",
    "Everything about INSR": """SELECT a.disease_id AS disease, ROUND(a.score, 3) AS score, a.evidence_types AS evidence,
       ROUND(e.log2FC, 3) AS log2FC, e.direction, e.p_value, e.significant
FROM gene g
JOIN associated_with a ON a.gene_id = g.id
LEFT JOIN expressed_in e ON e.gene_id = g.id AND e.disease_id = a.disease_id
WHERE g.symbol = 'INSR'
ORDER BY disease""",
    "Pathways with the most PCOS + T2D genes": """SELECT p.kegg_id, p.name AS pathway, COUNT(*) AS shared_genes
FROM pathway p
JOIN participates_in pi ON pi.pathway_id = p.id
WHERE p.category <> 'overview'
  AND pi.gene_id IN (SELECT gene_id FROM associated_with WHERE disease_id = 'PCOS')
  AND pi.gene_id IN (SELECT gene_id FROM associated_with WHERE disease_id = 'T2D')
GROUP BY p.id
ORDER BY shared_genes DESC
LIMIT 15""",
    "Pairwise gene overlap between diseases": """SELECT a1.disease_id AS disease_a, a2.disease_id AS disease_b, COUNT(*) AS shared_genes
FROM associated_with a1
JOIN associated_with a2 ON a1.gene_id = a2.gene_id AND a1.disease_id < a2.disease_id
GROUP BY 1, 2
ORDER BY shared_genes DESC""",
    "Drug-only links (the metformin effect)": """SELECT g.symbol AS gene, GROUP_CONCAT(a.disease_id) AS diseases
FROM associated_with a JOIN gene g ON g.id = a.gene_id
WHERE a.drug_only = 1
GROUP BY g.id
HAVING COUNT(*) >= 3
ORDER BY COUNT(*) DESC, gene
LIMIT 30""",
}


def show_result(df, truncated, ms):
    if df.empty:
        ui.empty_state("The query returned no rows.")
        return
    st.success(f"{len(df):,} row(s){' (truncated at 1,000)' if truncated else ''} in {ms:.0f} ms",
               icon=":material/check_circle:")
    st.dataframe(df, width="stretch", hide_index=True)
    st.download_button("Download as CSV", df.to_csv(index=False).encode(), "query_result.csv", "text/csv",
                       icon=":material/download:")


def render_sql_console():
    """Demo mode: read-only SQL against the bundled SQLite copy (same schema, as tables)."""
    st.info("**Live Cypher is only available when the app runs locally with Neo4j.** In demo mode you can "
            "query the bundled SQLite copy of the knowledge base with read-only SQL instead.",
            icon=":material/info:")
    st.caption("Tables: disease(id, full_name, source_id), gene(id, symbol, entrez_id, name), "
               "pathway(id, kegg_id, name, category, size), associated_with(disease_id, gene_id, score, "
               "evidence_types, drug_only), participates_in(gene_id, pathway_id), expressed_in(gene_id, "
               "disease_id, log2FC, direction, p_value, adj_p, significant). The database is opened read-only; "
               "only SELECT / WITH queries are accepted. Results are capped at 1,000 rows.")
    example = st.selectbox("Example queries", list(SQL_EXAMPLES), key="sq_example")
    if st.session_state.get("_sq_last_example") != example:
        st.session_state["sq_text"] = SQL_EXAMPLES[example]
        st.session_state["_sq_last_example"] = example
    query = st.text_area("SQL", height=220, key="sq_text")
    run = st.button("Run query", type="primary", icon=":material/play_arrow:")
    if not (run or st.session_state.get("_sq_autorun", True)):
        return
    st.session_state["_sq_autorun"] = False
    q = query.strip().rstrip(";")
    if not q:
        ui.empty_state("Enter a query.")
        return
    bad = kb.blocked_keywords(q, kb.SQL_WRITE_KEYWORDS)
    if bad or not q.lower().startswith(("select", "with")) or ";" in q:
        why = f"contains {', '.join(bad)}" if bad else "must be a single SELECT or WITH statement"
        st.error(f"Blocked: this console is read-only and the query {why}.", icon=":material/block:")
        return
    try:
        with st.spinner("Running query..."):
            t0 = time.perf_counter()
            df, truncated = kb.run_readonly_sql(q)
            ms = (time.perf_counter() - t0) * 1000
    except Exception as e:
        st.error(f"Query failed: {e}", icon=":material/error:")
        return
    show_result(df, truncated, ms)


def render():
    ui.page_header("Query Console", "Run read-only queries against the knowledge base")
    if kb.backend() != "neo4j":
        render_sql_console()
        return

    st.caption("Queries run inside a READ transaction and are also screened for write keywords "
               f"({', '.join(kb.WRITE_KEYWORDS)}). Results are capped at 1,000 rows. "
               "Schema: (:Disease)-[:ASSOCIATED_WITH {score, drug_only, evidence_types}]->(:Gene)"
               "-[:PARTICIPATES_IN]->(:Pathway); (:Gene)-[:EXPRESSED_IN {log2FC, direction, p_value, "
               "significant}]->(:Disease).")
    example = st.selectbox("Example queries", list(EXAMPLES), key="cq_example")
    if st.session_state.get("_cq_last_example") != example:
        st.session_state["cq_text"] = EXAMPLES[example]
        st.session_state["_cq_last_example"] = example
    query = st.text_area("Cypher", height=220, key="cq_text")
    run = st.button("Run query", type="primary", icon=":material/play_arrow:")

    if run or st.session_state.get("_cq_autorun", True):
        st.session_state["_cq_autorun"] = False
        if not query.strip():
            ui.empty_state("Enter a query.")
            return
        bad = kb.blocked_keywords(query)
        if bad:
            st.error(f"Blocked: this console is read-only and the query contains {', '.join(bad)}.",
                     icon=":material/block:")
            return
        try:
            with st.spinner("Running query..."):
                t0 = time.perf_counter()
                df, truncated = kb.run_readonly_cypher(query)
                ms = (time.perf_counter() - t0) * 1000
        except Exception as e:
            msg = getattr(e, "message", None) or str(e)
            st.error(f"Query failed: {msg}", icon=":material/error:")
            return
        show_result(df, truncated, ms)
