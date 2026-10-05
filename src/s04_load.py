"""Step 4: load Disease/Gene/Pathway nodes and ASSOCIATED_WITH / PARTICIPATES_IN / EXPRESSED_IN edges.

Primary store: Neo4j (bolt). A SQLite mirror with the same schema is always written so the KB is
usable without a running Neo4j; if Neo4j is unreachable the pipeline continues on SQLite + NetworkX.
"""
import sqlite3
import sys

import networkx as nx
import pandas as pd

from config import NEO4J_PASSWORD, NEO4J_URI, NEO4J_USER, PROC, RESULTS, SQLITE_PATH


def load_tables():
    diseases = pd.read_csv(PROC / "diseases.csv")
    assoc = pd.read_csv(PROC / "associations.csv")
    entrez = pd.read_csv(PROC / "gene_entrez.csv", dtype={"entrez_id": str})
    gp = pd.read_csv(PROC / "gene_pathway.csv", dtype={"entrez_id": str})
    pathways = pd.read_csv(PROC / "pathways.csv")
    universe = pd.read_csv(PROC / "kegg_universe.csv", dtype={"entrez_id": str})
    expr = pd.read_csv(PROC / "expression.csv")

    genes = (assoc[["ensembl_id", "symbol", "gene_name"]].drop_duplicates("ensembl_id")
             .merge(entrez[["symbol", "entrez_id"]], on="symbol", how="left")
             .rename(columns={"ensembl_id": "id", "gene_name": "name"}))
    # Pathway size = number of human genes in the pathway (all of KEGG, not just KB genes); used for enrichment.
    links = pd.read_csv(next((PROC.parent / "raw" / "kegg").glob("link_pathway_hsa.tsv")), sep="\t",
                        header=None, names=["g", "p"])
    links["kegg_id"] = links.p.str.replace("path:", "", regex=False)
    sizes = links.groupby("kegg_id").size().rename("size")
    pathways = pathways.merge(sizes, left_on="kegg_id", right_index=True, how="left").fillna({"size": 0})
    pathways["id"] = pathways.kegg_id
    pathways = pathways[pathways.kegg_id.isin(gp.kegg_id)]  # only pathways touched by a KB gene
    diseases["kegg_universe_size"] = len(universe)

    gp = gp.merge(genes[["id", "symbol"]], on="symbol")[["id", "kegg_id"]].rename(columns={"id": "gene_id"})
    assoc = assoc.rename(columns={"ensembl_id": "gene_id"})[
        ["disease", "gene_id", "score", "evidence_types", "drug_only", "literature_only"]]
    sym2id = genes.drop_duplicates("symbol").set_index("symbol").id
    expr["gene_id"] = expr.symbol.map(sym2id)
    expr = expr[["gene_id", "disease", "log2FC", "direction", "p_value", "adj_p", "significant"]]
    return diseases, genes, pathways, assoc, gp, expr


def load_sqlite(diseases, genes, pathways, assoc, gp, expr):
    SQLITE_PATH.unlink(missing_ok=True)
    con = sqlite3.connect(SQLITE_PATH)
    con.executescript("""
    CREATE TABLE disease (id TEXT PRIMARY KEY, name TEXT, full_name TEXT, source_id TEXT, efo_xref TEXT,
                          ot_total_associations INT, kegg_universe_size INT);
    CREATE TABLE gene (id TEXT PRIMARY KEY, symbol TEXT, entrez_id TEXT, name TEXT);
    CREATE TABLE pathway (id TEXT PRIMARY KEY, kegg_id TEXT, name TEXT, category TEXT, size INT);
    CREATE TABLE associated_with (disease_id TEXT REFERENCES disease, gene_id TEXT REFERENCES gene, score REAL,
                                  evidence_types TEXT, drug_only INT, literature_only INT,
                                  PRIMARY KEY (disease_id, gene_id));
    CREATE TABLE participates_in (gene_id TEXT REFERENCES gene, pathway_id TEXT REFERENCES pathway,
                                  PRIMARY KEY (gene_id, pathway_id));
    CREATE TABLE expressed_in (gene_id TEXT REFERENCES gene, disease_id TEXT REFERENCES disease, log2FC REAL,
                               direction TEXT, p_value REAL, adj_p REAL, significant INT,
                               PRIMARY KEY (gene_id, disease_id));
    """)
    diseases[["id", "name", "full_name", "source_id", "efo_xref", "ot_total_associations",
              "kegg_universe_size"]].to_sql("disease", con, if_exists="append", index=False)
    genes[["id", "symbol", "entrez_id", "name"]].to_sql("gene", con, if_exists="append", index=False)
    pathways[["id", "kegg_id", "name", "category", "size"]].to_sql("pathway", con, if_exists="append", index=False)
    assoc.rename(columns={"disease": "disease_id"}).astype({"drug_only": int, "literature_only": int}).to_sql("associated_with", con, if_exists="append", index=False)
    gp.rename(columns={"kegg_id": "pathway_id"}).to_sql("participates_in", con, if_exists="append", index=False)
    expr.rename(columns={"disease": "disease_id"}).astype({"significant": int}).to_sql(
        "expressed_in", con, if_exists="append", index=False)
    con.commit()
    con.close()
    print(f"SQLite written: {SQLITE_PATH}")


def load_neo4j(diseases, genes, pathways, assoc, gp, expr):
    from neo4j import GraphDatabase

    if not NEO4J_PASSWORD:
        raise RuntimeError("NEO4J_PASSWORD is not set")
    driver = GraphDatabase.driver(NEO4J_URI, auth=(NEO4J_USER, NEO4J_PASSWORD))
    driver.verify_connectivity()

    def rows(df):
        return df.astype(object).where(df.notna(), None).to_dict("records")

    with driver.session() as s:
        # Only remove this KB's own labels; never wipe unrelated data in the database.
        for label in ("Disease", "Gene", "Pathway"):
            s.run(f"MATCH (n:{label}) CALL (n) {{ DETACH DELETE n }} IN TRANSACTIONS OF 5000 ROWS").consume()
        for label in ("Disease", "Gene", "Pathway"):
            s.run(f"CREATE CONSTRAINT {label.lower()}_id IF NOT EXISTS FOR (n:{label}) REQUIRE n.id IS UNIQUE")
        s.run("CREATE INDEX gene_symbol IF NOT EXISTS FOR (g:Gene) ON (g.symbol)")

        s.run("UNWIND $rows AS r MERGE (d:Disease {id: r.id}) SET d += r", rows=rows(diseases))
        s.run("UNWIND $rows AS r MERGE (g:Gene {id: r.id}) SET g += r", rows=rows(genes))
        s.run("UNWIND $rows AS r MERGE (p:Pathway {id: r.id}) SET p += r",
              rows=rows(pathways[["id", "kegg_id", "name", "category", "size"]]))
        s.run("""UNWIND $rows AS r MATCH (d:Disease {id: r.disease}), (g:Gene {id: r.gene_id})
                 MERGE (d)-[a:ASSOCIATED_WITH]->(g)
                 SET a.score = r.score, a.evidence_types = split(r.evidence_types, ';'),
                     a.drug_only = r.drug_only, a.literature_only = r.literature_only""", rows=rows(assoc))
        s.run("""UNWIND $rows AS r MATCH (g:Gene {id: r.gene_id}), (p:Pathway {id: r.kegg_id})
                 MERGE (g)-[:PARTICIPATES_IN]->(p)""", rows=rows(gp))
        s.run("""UNWIND $rows AS r MATCH (g:Gene {id: r.gene_id}), (d:Disease {id: r.disease})
                 MERGE (g)-[e:EXPRESSED_IN]->(d)
                 SET e.log2FC = r.log2FC, e.direction = r.direction, e.p_value = r.p_value,
                     e.adj_p = r.adj_p, e.significant = r.significant""", rows=rows(expr))
        counts = s.run("""CALL () { MATCH (n) RETURN labels(n)[0] AS k, count(*) AS c
                                    UNION ALL MATCH ()-[r]->() RETURN type(r) AS k, count(*) AS c }
                          RETURN k, c""").data()
    driver.close()
    print("Neo4j loaded:", {c["k"]: c["c"] for c in counts})


def export_networkx(diseases, genes, pathways, assoc, gp, expr):
    G = nx.MultiDiGraph()
    for r in diseases.itertuples():
        G.add_node(f"D:{r.id}", kind="Disease", name=r.id)
    for r in genes.itertuples():
        G.add_node(r.id, kind="Gene", name=r.symbol)
    for r in pathways.itertuples():
        G.add_node(r.id, kind="Pathway", name=r.name)
    for r in assoc.itertuples():
        G.add_edge(f"D:{r.disease}", r.gene_id, type="ASSOCIATED_WITH", score=r.score, drug_only=bool(r.drug_only))
    for r in gp.itertuples():
        G.add_edge(r.gene_id, r.kegg_id, type="PARTICIPATES_IN")
    for r in expr.itertuples():
        G.add_edge(r.gene_id, f"D:{r.disease}", type="EXPRESSED_IN", log2FC=r.log2FC, direction=r.direction,
                   p_value=r.p_value, significant=bool(r.significant))
    nx.write_graphml(G, RESULTS / "comorbidity_kb.graphml")
    print(f"NetworkX graph: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")


def main():
    tables = load_tables()
    load_sqlite(*tables)
    export_networkx(*tables)
    try:
        load_neo4j(*tables)
        (RESULTS / "backend.txt").write_text("neo4j\n")
    except Exception as e:
        print(f"Neo4j unavailable ({e.__class__.__name__}: {e}); falling back to SQLite + NetworkX.",
              file=sys.stderr)
        (RESULTS / "backend.txt").write_text("sqlite\n")


if __name__ == "__main__":
    main()
