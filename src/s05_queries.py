"""Step 5: run the hub queries (Cypher on Neo4j, SQL on SQLite), cross-check them, write results/."""
import re
import sqlite3
import sys

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

from config import (MIN_ASSOC_SCORE, NEO4J_PASSWORD, NEO4J_URI, NEO4J_USER, PATHWAY_FDR, PROC, QUERIES,
                    RESULTS, SENSITIVITY_SCORES, SQLITE_PATH, VALIDATION_GENES)

DISEASES = ["NAFLD", "Obesity", "PCOS", "T2D"]
# (min_score, exclude_drug_only): primary run, score-threshold sensitivity, and drug-evidence robustness.
SCENARIOS = [(s, False) for s in SENSITIVITY_SCORES] + [(MIN_ASSOC_SCORE, True)]


def parse(path, marker):
    """Split a query file into {name: query} using '<marker> name: X' header lines."""
    out, name, buf = {}, None, []
    for line in path.read_text().splitlines():
        m = re.match(rf"^{re.escape(marker)}\s*name:\s*(\w+)", line)
        if m:
            if name:
                out[name] = "\n".join(buf).strip().rstrip(";")
            name, buf = m.group(1), []
        elif name:
            buf.append(line)
    if name:
        out[name] = "\n".join(buf).strip().rstrip(";")
    return out


class SQLBackend:
    def __init__(self):
        self.con = sqlite3.connect(SQLITE_PATH)
        self.q = parse(QUERIES / "hubs.sql", "--")

    def run(self, name, **params):
        sql = self.q[name]
        used = {k: int(v) if isinstance(v, bool) else v for k, v in params.items() if f":{k}" in sql}
        return pd.read_sql_query(sql, self.con, params=used)


class Neo4jBackend:
    def __init__(self):
        from neo4j import GraphDatabase
        if not NEO4J_PASSWORD:
            raise RuntimeError("NEO4J_PASSWORD is not set")
        self.driver = GraphDatabase.driver(NEO4J_URI, auth=(NEO4J_USER, NEO4J_PASSWORD))
        self.driver.verify_connectivity()
        self.q = parse(QUERIES / "hubs.cypher", "//")

    def run(self, name, **params):
        with self.driver.session() as s:
            cy = self.q[name]
            return pd.DataFrame(s.run(cy, **{k: v for k, v in params.items() if f"${k}" in cy}).data())


def listify(x):
    return ",".join(map(str, x)) if isinstance(x, list) else x


def crosscheck(neo, sql, name, key_cols, **params):
    """Assert the Cypher and SQL versions of a query return the same rows (on key columns)."""
    a = neo.run(name, **params).map(listify)[key_cols]
    b = sql.run(name, **params)[key_cols]
    norm = lambda df: df.astype(str).sort_values(key_cols).reset_index(drop=True)
    same = norm(a).equals(norm(b))
    return f"{name}({params}): neo4j={len(a)} rows, sqlite={len(b)} rows, " \
           f"identical={same}", same


def pathway_enrichment(db, min_score, exclude_drug_only=False):
    """Hypergeometric over-representation of each disease's genes in each KEGG pathway."""
    kw = dict(min_score=min_score, exclude_drug_only=exclude_drug_only)
    counts = db.run("pathway_disease_counts", **kw)
    denom = db.run("disease_kegg_gene_counts", **kw).set_index("disease")
    universe = int(denom.universe.iloc[0])
    counts["n_disease_genes"] = counts.disease.map(denom.n_genes_in_kegg)
    # P(X >= k), X ~ Hypergeom(N=universe, K=pathway size, n=disease genes in KEGG)
    counts["p_value"] = stats.hypergeom.sf(counts.k - 1, universe, counts.pathway_size, counts.n_disease_genes)
    counts["expected"] = counts.n_disease_genes * counts.pathway_size / universe
    counts["fold_enrichment"] = counts.k / counts.expected
    # BH within each disease over every non-overview pathway (untested pathways have k=0 -> p=1).
    n_tests = db.run("pathway_disease_counts", min_score=0.0, exclude_drug_only=False).kegg_id.nunique()
    counts["fdr"] = np.nan
    for d, idx in counts.groupby("disease").groups.items():
        p = counts.loc[idx, "p_value"].to_numpy()
        padded = np.concatenate([p, np.ones(max(n_tests - len(p), 0))])
        counts.loc[idx, "fdr"] = multipletests(padded, method="fdr_bh")[1][:len(p)]
    counts["enriched"] = counts.fdr < PATHWAY_FDR
    return counts


def hub_pathways(enr):
    wide = enr.pivot_table(index=["kegg_id", "pathway", "category", "pathway_size"], columns="disease",
                           values=["k", "fdr"], aggfunc="first")
    wide.columns = [f"{a}_{b}" for a, b in wide.columns]
    wide = wide.reset_index()
    for d in DISEASES:
        wide[f"k_{d}"] = wide.get(f"k_{d}", pd.Series(0, index=wide.index)).fillna(0).astype(int)
        wide[f"fdr_{d}"] = wide.get(f"fdr_{d}", pd.Series(1.0, index=wide.index)).fillna(1.0)
    enriched = enr[enr.enriched].groupby("kegg_id").disease.apply(lambda s: ",".join(sorted(s)))
    wide["enriched_in"] = wide.kegg_id.map(enriched).fillna("")
    wide["n_enriched"] = wide.enriched_in.apply(lambda s: len(s.split(",")) if s else 0)
    wide["min_fdr"] = wide[[f"fdr_{d}" for d in DISEASES]].min(axis=1)
    wide = wide.sort_values(["n_enriched", "min_fdr"], ascending=[False, True])
    return wide


def main():
    sql = SQLBackend()
    try:
        neo = Neo4jBackend()
        db, backend = neo, "neo4j"
    except Exception as e:
        print(f"Neo4j unavailable ({e}); using SQLite backend only.", file=sys.stderr)
        neo, db, backend = None, sql, "sqlite"

    log = [f"primary backend: {backend}"]
    if neo:
        checks = [
            ("hub_genes", ["symbol", "n_diseases", "diseases"]),
            ("pathway_disease_counts", ["kegg_id", "disease", "k"]),
            ("disease_kegg_gene_counts", ["disease", "n_genes_in_kegg"]),
            ("hub_pathways_naive", ["kegg_id", "n_diseases"]),
            ("expression_consistent_hubs", ["symbol", "n_diseases", "n_deg", "deg_diseases", "direction"]),
            ("expression_discordant_hubs", ["symbol", "n_diseases", "n_deg"]),
            ("shared_genes_pairwise", ["disease_a", "disease_b", "shared_genes"]),
        ]
        ok = True
        for score, excl in SCENARIOS:
            for name, cols in checks:
                msg, same = crosscheck(neo, sql, name, cols, min_score=score, exclude_drug_only=excl)
                ok &= same
                log.append(msg)
        log.append(f"ALL CYPHER/SQL CROSS-CHECKS PASSED: {ok}")
    (RESULTS / "backend_crosscheck.txt").write_text("\n".join(log) + "\n")
    print("\n".join(log[-3:]))

    expr = pd.read_csv(PROC / "expression.csv")
    lfc = expr.pivot(index="symbol", columns="disease", values="log2FC")
    pv = expr.pivot(index="symbol", columns="disease", values="p_value")
    sig = expr.pivot(index="symbol", columns="disease", values="significant")

    sens = []
    for score, excl in SCENARIOS:
        kw = dict(min_score=score, exclude_drug_only=excl)
        hubs = db.run("hub_genes", **kw).map(listify)
        cons = db.run("expression_consistent_hubs", **kw).map(listify)
        disc = db.run("expression_discordant_hubs", **kw).map(listify)
        enr = pathway_enrichment(db, score, excl)
        hp = hub_pathways(enr)
        row = {"min_score": score, "exclude_drug_only": excl, "hub_genes": len(hubs),
               "hub_genes_4_diseases": int((hubs.n_diseases == 4).sum()),
               "hub_genes_3plus": int((hubs.n_diseases >= 3).sum()),
               "expression_consistent_hubs": len(cons), "expression_discordant_hubs": len(disc),
               "hub_pathways_enriched_2plus": int((hp.n_enriched >= 2).sum()),
               "hub_pathways_naive_2plus": len(db.run("hub_pathways_naive", **kw))}
        for g in VALIDATION_GENES:
            h = hubs[hubs.symbol == g]
            row[f"{g}_n_diseases"] = int(h.n_diseases.iloc[0]) if len(h) else 0
            row[f"{g}_expr_consistent"] = g in set(cons.symbol)
        sens.append(row)

        if excl:
            hubs.to_csv(RESULTS / "hub_genes_excluding_drug_only.csv", index=False)
            cons.to_csv(RESULTS / "expression_consistent_hubs_excluding_drug_only.csv", index=False)
            hp[hp.n_enriched >= 2].to_csv(RESULTS / "hub_pathways_excluding_drug_only.csv", index=False)
        elif score == MIN_ASSOC_SCORE:
            for d in DISEASES:
                hubs[f"log2FC_{d}"] = hubs.symbol.map(lfc.get(d, pd.Series(dtype=float)))
                hubs[f"p_{d}"] = hubs.symbol.map(pv.get(d, pd.Series(dtype=float)))
                hubs[f"sig_{d}"] = hubs.symbol.map(sig.get(d, pd.Series(dtype=object)))
            assoc = pd.read_csv(PROC / "associations.csv")
            drug = assoc[assoc.score >= score].groupby("symbol").drug_only.agg(["sum", "count"])
            hubs["n_drug_only_edges"] = hubs.symbol.map(drug["sum"]).fillna(0).astype(int)
            hubs["hub_without_drug_only"] = (hubs.n_diseases - hubs.n_drug_only_edges) >= 2
            hubs["expression_class"] = np.select(
                [hubs.symbol.isin(cons.symbol), hubs.symbol.isin(disc.symbol)],
                ["consistent", "discordant"], default="insufficient")
            hubs.to_csv(RESULTS / "hub_genes.csv", index=False)
            cons.to_csv(RESULTS / "expression_consistent_hubs.csv", index=False)
            disc.to_csv(RESULTS / "expression_discordant_hubs.csv", index=False)
            enr.to_csv(RESULTS / "pathway_enrichment.csv", index=False)
            hp.to_csv(RESULTS / "pathways_by_disease_enrichment.csv", index=False)
            hp[hp.n_enriched >= 2].to_csv(RESULTS / "hub_pathways.csv", index=False)
            db.run("shared_genes_pairwise", **kw).to_csv(RESULTS / "pairwise_overlap.csv", index=False)
            print(f"min_score={score}: {len(hubs)} hub genes, {len(cons)} expression-consistent, "
                  f"{len(disc)} discordant, {(hp.n_enriched >= 2).sum()} hub pathways")
    pd.DataFrame(sens).to_csv(RESULTS / "sensitivity.csv", index=False)
    print(pd.DataFrame(sens).T.to_string())


if __name__ == "__main__":
    main()
