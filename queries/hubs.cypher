// Comorbidity-Cluster KB: example Cypher queries.
// Each query starts with a "// name:" line; parameters use $param syntax.
// Run in Neo4j Browser (http://localhost:7474) after
//   :param min_score => 0.1
//   :param exclude_drug_only => false
// or via src/s05_queries.py. exclude_drug_only=true drops edges supported only by drug (ChEMBL) evidence.

// name: hub_genes
// Hub gene = associated (score >= $min_score) with 2+ of the 4 diseases (degree centrality).
MATCH (d:Disease)-[a:ASSOCIATED_WITH]->(g:Gene)
WHERE a.score >= $min_score AND NOT (a.drug_only AND $exclude_drug_only)
WITH g, collect(d.id) AS diseases, avg(a.score) AS mean_score, max(a.score) AS max_score
WHERE size(diseases) >= 2
WITH g, diseases, mean_score, max_score ORDER BY g.symbol
UNWIND diseases AS dd WITH g, dd, mean_score, max_score ORDER BY dd
WITH g, collect(dd) AS diseases, mean_score, max_score
RETURN g.symbol AS symbol, g.entrez_id AS entrez_id, size(diseases) AS n_diseases,
       diseases, round(mean_score, 4) AS mean_score, round(max_score, 4) AS max_score
ORDER BY n_diseases DESC, mean_score DESC, symbol;

// name: pathway_disease_counts
// Per pathway x disease: how many of the disease's associated genes fall in the pathway.
// Hub pathways are derived from these counts in s05_queries.py (hypergeometric enrichment, BH FDR).
MATCH (d:Disease)-[a:ASSOCIATED_WITH]->(g:Gene)-[:PARTICIPATES_IN]->(p:Pathway)
WHERE a.score >= $min_score AND NOT (a.drug_only AND $exclude_drug_only) AND p.category <> 'overview'
RETURN p.kegg_id AS kegg_id, p.name AS pathway, p.category AS category, p.size AS pathway_size,
       d.id AS disease, count(DISTINCT g) AS k
ORDER BY kegg_id, disease;

// name: disease_kegg_gene_counts
// Number of each disease's associated genes that have at least one KEGG pathway (enrichment denominator).
MATCH (d:Disease)-[a:ASSOCIATED_WITH]->(g:Gene)
WHERE a.score >= $min_score AND NOT (a.drug_only AND $exclude_drug_only) AND EXISTS { (g)-[:PARTICIPATES_IN]->(:Pathway) }
RETURN d.id AS disease, count(DISTINCT g) AS n_genes_in_kegg, d.kegg_universe_size AS universe
ORDER BY disease;

// name: hub_pathways_naive
// Naive pathway degree: number of diseases with >= 1 associated gene in the pathway (no enrichment).
MATCH (d:Disease)-[a:ASSOCIATED_WITH]->(:Gene)-[:PARTICIPATES_IN]->(p:Pathway)
WHERE a.score >= $min_score AND NOT (a.drug_only AND $exclude_drug_only)
WITH p, collect(DISTINCT d.id) AS diseases
WHERE size(diseases) >= 2
RETURN p.kegg_id AS kegg_id, p.name AS pathway, size(diseases) AS n_diseases, diseases
ORDER BY n_diseases DESC, kegg_id;

// name: expression_consistent_hubs
// Hub genes that are significant DEGs in 2+ of their associated diseases, all in the same direction.
MATCH (d:Disease)-[a:ASSOCIATED_WITH]->(g:Gene)
WHERE a.score >= $min_score AND NOT (a.drug_only AND $exclude_drug_only)
WITH g, collect(d.id) AS diseases
WHERE size(diseases) >= 2
MATCH (g)-[e:EXPRESSED_IN]->(d2:Disease)
WHERE d2.id IN diseases AND e.significant
WITH g, diseases, d2, e ORDER BY d2.id
WITH g, diseases, collect(d2.id) AS deg_diseases, collect(DISTINCT e.direction) AS dirs,
     collect(round(e.log2FC, 3)) AS log2FCs
WHERE size(deg_diseases) >= 2 AND size(dirs) = 1
RETURN g.symbol AS symbol, size(diseases) AS n_diseases, diseases, deg_diseases,
       size(deg_diseases) AS n_deg, dirs[0] AS direction, log2FCs
ORDER BY n_deg DESC, n_diseases DESC, symbol;

// name: expression_discordant_hubs
// Hub genes significant in 2+ associated diseases but moving in opposite directions.
MATCH (d:Disease)-[a:ASSOCIATED_WITH]->(g:Gene)
WHERE a.score >= $min_score AND NOT (a.drug_only AND $exclude_drug_only)
WITH g, collect(d.id) AS diseases
WHERE size(diseases) >= 2
MATCH (g)-[e:EXPRESSED_IN]->(d2:Disease)
WHERE d2.id IN diseases AND e.significant
WITH g, diseases, d2, e ORDER BY d2.id
WITH g, diseases, collect(d2.id) AS deg_diseases, collect(DISTINCT e.direction) AS dirs,
     collect(e.direction) AS directions
WHERE size(deg_diseases) >= 2 AND size(dirs) > 1
RETURN g.symbol AS symbol, size(diseases) AS n_diseases, diseases, deg_diseases,
       size(deg_diseases) AS n_deg, directions
ORDER BY n_deg DESC, n_diseases DESC, symbol;

// name: gene_profile
// Everything the KB knows about one gene (used for validation of INSR / IRS1 / IGF1).
MATCH (g:Gene {symbol: $symbol})
OPTIONAL MATCH (d:Disease)-[a:ASSOCIATED_WITH]->(g)
OPTIONAL MATCH (g)-[e:EXPRESSED_IN]->(d)
RETURN g.symbol AS symbol, d.id AS disease, round(a.score, 4) AS score,
       round(e.log2FC, 3) AS log2FC, e.direction AS direction, e.p_value AS p_value, e.significant AS significant
ORDER BY disease;

// name: shared_genes_pairwise
// Size of the gene overlap between every pair of diseases (the comorbidity overlap, computed at query time).
MATCH (d1:Disease)-[a1:ASSOCIATED_WITH]->(g:Gene)<-[a2:ASSOCIATED_WITH]-(d2:Disease)
WHERE d1.id < d2.id AND a1.score >= $min_score AND a2.score >= $min_score
RETURN d1.id AS disease_a, d2.id AS disease_b, count(g) AS shared_genes
ORDER BY disease_a, disease_b;
