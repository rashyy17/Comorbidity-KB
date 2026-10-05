-- Comorbidity-Cluster KB: SQL equivalents of queries/hubs.cypher (SQLite backend, results/comorbidity_kb.sqlite).
-- Each query starts with a "-- name:" line; :min_score / :exclude_drug_only (0/1) / :symbol are named parameters.
-- Lists are returned as comma-separated strings sorted alphabetically.

-- name: hub_genes
SELECT g.symbol, g.entrez_id, COUNT(*) AS n_diseases,
       (SELECT GROUP_CONCAT(disease_id, ',') FROM (SELECT a2.disease_id FROM associated_with a2
         WHERE a2.gene_id = g.id AND a2.score >= :min_score AND NOT (a2.drug_only = 1 AND :exclude_drug_only = 1) ORDER BY a2.disease_id)) AS diseases,
       ROUND(AVG(a.score), 4) AS mean_score, ROUND(MAX(a.score), 4) AS max_score
FROM associated_with a JOIN gene g ON g.id = a.gene_id
WHERE a.score >= :min_score AND NOT (a.drug_only = 1 AND :exclude_drug_only = 1)
GROUP BY g.id
HAVING COUNT(*) >= 2
ORDER BY n_diseases DESC, AVG(a.score) DESC, g.symbol;

-- name: pathway_disease_counts
SELECT p.kegg_id, p.name AS pathway, p.category, p.size AS pathway_size, a.disease_id AS disease,
       COUNT(DISTINCT a.gene_id) AS k
FROM associated_with a
JOIN participates_in pi ON pi.gene_id = a.gene_id
JOIN pathway p ON p.id = pi.pathway_id
WHERE a.score >= :min_score AND NOT (a.drug_only = 1 AND :exclude_drug_only = 1) AND p.category <> 'overview'
GROUP BY p.kegg_id, a.disease_id
ORDER BY p.kegg_id, a.disease_id;

-- name: disease_kegg_gene_counts
SELECT a.disease_id AS disease, COUNT(DISTINCT a.gene_id) AS n_genes_in_kegg, d.kegg_universe_size AS universe
FROM associated_with a JOIN disease d ON d.id = a.disease_id
WHERE a.score >= :min_score AND NOT (a.drug_only = 1 AND :exclude_drug_only = 1) AND a.gene_id IN (SELECT gene_id FROM participates_in)
GROUP BY a.disease_id
ORDER BY a.disease_id;

-- name: hub_pathways_naive
SELECT p.kegg_id, p.name AS pathway, COUNT(DISTINCT a.disease_id) AS n_diseases
FROM associated_with a
JOIN participates_in pi ON pi.gene_id = a.gene_id
JOIN pathway p ON p.id = pi.pathway_id
WHERE a.score >= :min_score AND NOT (a.drug_only = 1 AND :exclude_drug_only = 1)
GROUP BY p.kegg_id
HAVING COUNT(DISTINCT a.disease_id) >= 2
ORDER BY n_diseases DESC, p.kegg_id;

-- name: expression_consistent_hubs
WITH hubs AS (
  SELECT gene_id, COUNT(*) AS n_diseases FROM associated_with
  WHERE score >= :min_score AND NOT (drug_only = 1 AND :exclude_drug_only = 1) GROUP BY gene_id HAVING COUNT(*) >= 2
), degs AS (
  SELECT e.gene_id, e.disease_id, e.direction, e.log2FC
  FROM expressed_in e
  JOIN associated_with a ON a.gene_id = e.gene_id AND a.disease_id = e.disease_id AND a.score >= :min_score AND NOT (a.drug_only = 1 AND :exclude_drug_only = 1)
  WHERE e.significant = 1
)
SELECT g.symbol, h.n_diseases, COUNT(*) AS n_deg,
       (SELECT GROUP_CONCAT(disease_id, ',') FROM (SELECT disease_id FROM degs d2
         WHERE d2.gene_id = g.id ORDER BY disease_id)) AS deg_diseases,
       MIN(d.direction) AS direction
FROM hubs h JOIN degs d ON d.gene_id = h.gene_id JOIN gene g ON g.id = h.gene_id
GROUP BY h.gene_id
HAVING COUNT(*) >= 2 AND COUNT(DISTINCT d.direction) = 1
ORDER BY n_deg DESC, h.n_diseases DESC, g.symbol;

-- name: expression_discordant_hubs
WITH hubs AS (
  SELECT gene_id, COUNT(*) AS n_diseases FROM associated_with
  WHERE score >= :min_score AND NOT (drug_only = 1 AND :exclude_drug_only = 1) GROUP BY gene_id HAVING COUNT(*) >= 2
), degs AS (
  SELECT e.gene_id, e.disease_id, e.direction
  FROM expressed_in e
  JOIN associated_with a ON a.gene_id = e.gene_id AND a.disease_id = e.disease_id AND a.score >= :min_score AND NOT (a.drug_only = 1 AND :exclude_drug_only = 1)
  WHERE e.significant = 1
)
SELECT g.symbol, h.n_diseases, COUNT(*) AS n_deg
FROM hubs h JOIN degs d ON d.gene_id = h.gene_id JOIN gene g ON g.id = h.gene_id
GROUP BY h.gene_id
HAVING COUNT(*) >= 2 AND COUNT(DISTINCT d.direction) > 1
ORDER BY n_deg DESC, h.n_diseases DESC, g.symbol;

-- name: gene_profile
SELECT g.symbol, a.disease_id AS disease, ROUND(a.score, 4) AS score, ROUND(e.log2FC, 3) AS log2FC,
       e.direction, e.p_value, e.significant
FROM gene g
JOIN associated_with a ON a.gene_id = g.id
LEFT JOIN expressed_in e ON e.gene_id = g.id AND e.disease_id = a.disease_id
WHERE g.symbol = :symbol
ORDER BY a.disease_id;

-- name: shared_genes_pairwise
SELECT a1.disease_id AS disease_a, a2.disease_id AS disease_b, COUNT(*) AS shared_genes
FROM associated_with a1 JOIN associated_with a2 ON a1.gene_id = a2.gene_id AND a1.disease_id < a2.disease_id
WHERE a1.score >= :min_score AND a2.score >= :min_score
GROUP BY a1.disease_id, a2.disease_id
ORDER BY disease_a, disease_b;
