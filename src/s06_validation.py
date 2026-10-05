"""Step 6: validation.

1. Literature genes (INSR, IRS1, IGF1): are they hubs? expression-consistent? on what evidence?
2. Unbiased reference set: are hub genes enriched for KEGG hsa04931 (Insulin resistance) members?
3. Is expression consistency among hubs better than chance? (gene-label permutation within each dataset)
4. Genome-wide log2FC correlation between datasets (do the two adipose datasets agree?).
"""
import json

import numpy as np
import pandas as pd
from scipy import stats

from config import DEG_P, PROC, RESULTS, SENSITIVITY_SCORES, VALIDATION_GENES

DISEASES = ["NAFLD", "Obesity", "PCOS", "T2D"]
RNG = np.random.default_rng(20261005)
N_PERM = 2000


def validation_genes(assoc, expr, sens):
    rows = []
    for g in VALIDATION_GENES:
        for d in DISEASES:
            a = assoc[(assoc.symbol == g) & (assoc.disease == d)]
            e = expr[(expr.symbol == g) & (expr.disease == d)]
            de = pd.read_csv(PROC / f"de_{d}.csv", index_col=0)
            rows.append({
                "gene": g, "disease": d,
                "ot_score": round(float(a.score.iloc[0]), 3) if len(a) else None,
                "evidence_types": a.evidence_types.iloc[0] if len(a) else None,
                "measured_in_geo": g in de.index,
                "log2FC": round(float(de.loc[g, "log2FC"]), 3) if g in de.index else None,
                "p_value": float(de.loc[g, "p_value"]) if g in de.index else None,
                "significant": bool(e.significant.iloc[0]) if len(e) else (
                    bool(de.loc[g, "p_value"] < DEG_P) if g in de.index else None),
            })
    table = pd.DataFrame(rows)
    table.to_csv(RESULTS / "validation_genes.csv", index=False)
    summary = {}
    for g in VALIDATION_GENES:
        t = table[table.gene == g]
        summary[g] = {
            "degree_by_threshold": {f"{s}": int(((t.ot_score.fillna(0) >= s)).sum()) for s in SENSITIVITY_SCORES},
            "hub_at_primary_threshold": bool((t.ot_score.fillna(0) >= 0.1).sum() >= 2),
            "expression_consistent": bool(sens.loc[sens.min_score.eq(0.1) & ~sens.exclude_drug_only,
                                                   f"{g}_expr_consistent"].iloc[0]),
        }
    return table, summary


def insulin_resistance_enrichment(hubs, assoc):
    gp = pd.read_csv(PROC / "gene_pathway.csv")
    kb_genes = set(assoc.symbol) & set(gp.symbol)  # KB genes that are in KEGG at all
    ir = set(gp[gp.kegg_id == "hsa04931"].symbol)
    hub = set(hubs.symbol) & kb_genes
    a = len(hub & ir); b = len(hub - ir)
    c = len((kb_genes - hub) & ir); d = len((kb_genes - hub) - ir)
    odds, p = stats.fisher_exact([[a, b], [c, d]], alternative="greater")
    return {"reference_set": "KEGG hsa04931 Insulin resistance (genes in the KB)",
            "hub_in_set": a, "hub_not_in_set": b, "nonhub_in_set": c, "nonhub_not_in_set": d,
            "odds_ratio": round(odds, 2), "p_value": float(p),
            "hub_members": sorted(hub & ir)}


def consistency_counts(hub_diseases, sig, dirn):
    """hub_diseases: {gene: [diseases]}; sig/dirn: {disease: {gene: bool/str}}."""
    cons = disc = 0
    for g, ds in hub_diseases.items():
        dirs = [dirn[d][g] for d in ds if sig[d].get(g, False)]
        if len(dirs) >= 2:
            if len(set(dirs)) == 1:
                cons += 1
            else:
                disc += 1
    return cons, disc


def permutation_test(hubs, assoc):
    hub_diseases = {r.symbol: r.diseases.split(",") for r in hubs.itertuples()}
    de = {d: pd.read_csv(PROC / f"de_{d}.csv", index_col=0) for d in DISEASES}
    sig = {d: (de[d].p_value < DEG_P).to_dict() for d in DISEASES}
    dirn = {d: np.where(de[d].log2FC > 0, "up", "down") for d in DISEASES}
    dirn = {d: dict(zip(de[d].index, dirn[d])) for d in DISEASES}
    obs_c, obs_d = consistency_counts(hub_diseases, sig, dirn)

    null_c, null_rate = [], []
    for _ in range(N_PERM):
        s2, d2 = {}, {}
        for d in DISEASES:
            idx = de[d].index.to_numpy()
            perm = RNG.permutation(len(idx))
            s2[d] = dict(zip(idx, np.array(list(sig[d].values()))[perm]))
            d2[d] = dict(zip(idx, np.array(list(dirn[d].values()))[perm]))
        c, dd = consistency_counts(hub_diseases, s2, d2)
        null_c.append(c)
        null_rate.append(c / (c + dd) if c + dd else np.nan)
    null_c, null_rate = np.array(null_c), np.array(null_rate)
    obs_rate = obs_c / (obs_c + obs_d)
    return {
        "observed_consistent": obs_c, "observed_discordant": obs_d, "observed_concordance": round(obs_rate, 3),
        "null_consistent_mean": round(float(null_c.mean()), 1),
        "null_consistent_95pct": [int(np.percentile(null_c, 2.5)), int(np.percentile(null_c, 97.5))],
        "p_consistent_count": float((np.sum(null_c >= obs_c) + 1) / (N_PERM + 1)),
        "null_concordance_mean": round(float(np.nanmean(null_rate)), 3),
        "p_concordance": float((np.sum(null_rate >= obs_rate) + 1) / (N_PERM + 1)),
        "n_permutations": N_PERM,
    }


def background_concordance(perm_result):
    """Concordance among ALL genes significant in 2+ datasets (no Open Targets filter).

    The permutation null above breaks the genome-wide correlation between datasets; this baseline keeps it,
    so it is the fairer reference for "is 82% concordance special to hub genes?".
    """
    de = {d: pd.read_csv(PROC / f"de_{d}.csv", index_col=0) for d in DISEASES}
    sig = pd.concat({d: de[d].p_value < DEG_P for d in DISEASES}, axis=1).fillna(False)
    up = pd.concat({d: de[d].log2FC > 0 for d in DISEASES}, axis=1)
    multi = sig.sum(axis=1) >= 2
    cons = disc = 0
    for g in sig.index[multi]:
        dirs = set(up.loc[g, sig.loc[g]].tolist())
        cons += len(dirs) == 1
        disc += len(dirs) > 1
    rate = cons / (cons + disc)
    k_obs = perm_result["observed_consistent"]
    n_obs = k_obs + perm_result["observed_discordant"]
    p = stats.binomtest(k_obs, n_obs, rate, alternative="greater").pvalue
    return {"all_genes_sig_in_2plus": int(cons + disc), "consistent": int(cons), "discordant": int(disc),
            "background_concordance": round(rate, 3), "hub_concordance": round(k_obs / n_obs, 3),
            "binomial_p_hub_vs_background": float(p)}


def dataset_correlations():
    lfc = pd.concat({d: pd.read_csv(PROC / f"de_{d}.csv", index_col=0).log2FC for d in DISEASES}, axis=1)
    out = []
    for i, a in enumerate(DISEASES):
        for b in DISEASES[i + 1:]:
            both = lfc[[a, b]].dropna()
            rho, p = stats.spearmanr(both[a], both[b])
            out.append({"a": a, "b": b, "n_genes": len(both), "spearman_rho": round(rho, 3), "p_value": float(p)})
    pd.DataFrame(out).to_csv(RESULTS / "dataset_log2fc_correlation.csv", index=False)
    return out


def main():
    assoc = pd.read_csv(PROC / "associations.csv")
    expr = pd.read_csv(PROC / "expression.csv")
    hubs = pd.read_csv(RESULTS / "hub_genes.csv")
    sens = pd.read_csv(RESULTS / "sensitivity.csv")

    table, gene_summary = validation_genes(assoc, expr, sens)
    result = {
        "validation_genes": gene_summary,
        "insulin_resistance_enrichment": insulin_resistance_enrichment(hubs, assoc),
        "expression_consistency_permutation": (perm := permutation_test(hubs, assoc)),
        "expression_consistency_background": background_concordance(perm),
        "dataset_log2fc_correlation": dataset_correlations(),
    }
    (RESULTS / "validation.json").write_text(json.dumps(result, indent=2))
    print(table.to_string(index=False))
    print(json.dumps(result, indent=1))


if __name__ == "__main__":
    main()
