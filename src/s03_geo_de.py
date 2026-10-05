"""Step 3: download GEO series, verify groups, map probes to genes, run moderated-t DE.

Writes data/processed/de_<DISEASE>.csv (all measured genes) and expression.csv
(Open Targets disease-gene pairs that were measured, with log2FC/direction/p).
"""
import json
import logging
import time

import GEOparse
import numpy as np
import pandas as pd
import requests

from config import DEG_ABS_LFC, DEG_P, PROC, RAW, RESULTS
from limma_py import moderated_t

logging.getLogger("GEOparse").setLevel(logging.ERROR)
GEO_DIR = RAW / "geo"

# For each disease: series, the characteristic that encodes the group, and which values are case/control.
DATASETS = {
    "PCOS": dict(gse="GSE34526", field="disease", tissue="ovarian granulosa cells",
                 case=["Polycystic ovary Syndrome"], control=["normal"], scale="linear"),
    "T2D": dict(gse="GSE16415", field="disease state", tissue="visceral (omental) adipose tissue",
                case=["type-2 diabetes mellitus"], control=["control"], scale="linear"),
    "NAFLD": dict(gse="GSE48452", field="group", tissue="liver biopsy",
                  case=["Steatosis", "Nash"], control=["Control"], scale="log2",
                  sensitivity_control=["Healthy obese"]),
    "Obesity": dict(gse="GSE55200", field="disease state", tissue="subcutaneous adipose tissue",
                    case=["metabolically healthy obese MHO", "metabolically unhealthy obese MUO"],
                    control=["lean healthy LH"], scale="log2"),
}

LOW_EXPR_QUANTILE = 0.20  # drop probes whose mean expression is in the bottom 20%


def download(gse):
    """GEOparse's FTP downloader is flaky (size-mismatch errors); fetch the SOFT file over HTTPS."""
    path = GEO_DIR / f"{gse}_family.soft.gz"
    if path.exists():
        return path
    url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{gse[:-3]}nnn/{gse}/soft/{gse}_family.soft.gz"
    for attempt in range(6):
        try:
            with requests.get(url, stream=True, timeout=300, headers={"User-Agent": "Mozilla/5.0"}) as r:
                r.raise_for_status()
                tmp = path.with_suffix(".part")
                with open(tmp, "wb") as fh:
                    for chunk in r.iter_content(1 << 20):
                        fh.write(chunk)
            tmp.rename(path)
            return path
        except Exception as e:
            print(f"  {gse} download retry {attempt + 1}: {e}")
            time.sleep(15)
    raise RuntimeError(f"could not download {gse}")


def characteristic(gsm, field):
    for c in gsm.metadata.get("characteristics_ch1", []):
        k, _, v = c.partition(":")
        if k.strip().lower() == field.lower():
            return v.strip()
    return None


def probe_to_symbol(gpl):
    t = gpl.table.copy()
    t["ID"] = t["ID"].astype(str)
    if "Gene Symbol" in t.columns:  # GPL570 (Affymetrix) / GPL2986 (ABI)
        sym = t["Gene Symbol"].astype(str)
        # Drop probes annotated to more than one gene.
        sym = sym.where(~sym.str.contains("///") & (sym != "nan") & (sym.str.strip() != ""))
    else:  # Affymetrix Gene ST arrays: "NM_x // SYMBOL // desc // loc // entrez /// ..."
        ga = t["gene_assignment"].astype(str)
        sym = ga.str.split("///").str[0].str.split("//").str[1].str.strip()
        sym = sym.where(~ga.isin(["---", "nan"]))
        if "category" in t.columns:  # keep only main design probesets (drops controls)
            sym = sym.where(t["category"].astype(str).eq("main"))
    return pd.Series(sym.values, index=t["ID"]).dropna()


def quantile_normalize(m):
    ranks = m.rank(method="first")
    mean_sorted = np.sort(m.to_numpy(), axis=0).mean(axis=1)
    return ranks.apply(lambda col: pd.Series(mean_sorted[col.astype(int).to_numpy() - 1], index=col.index))


def prepare(gse_obj, cfg):
    expr = gse_obj.pivot_samples("VALUE").astype(float)
    expr.index = expr.index.astype(str)
    if cfg["scale"] == "linear":
        floor = np.nanpercentile(expr.to_numpy()[expr.to_numpy() > 0], 1)
        expr = np.log2(expr.clip(lower=floor))
        expr = quantile_normalize(expr)
    expr = expr.dropna()
    keep = expr.mean(axis=1) > expr.mean(axis=1).quantile(LOW_EXPR_QUANTILE)
    expr = expr[keep]
    sym = probe_to_symbol(next(iter(gse_obj.gpls.values())))
    expr = expr.loc[expr.index.intersection(sym.index)]
    expr.insert(0, "symbol", sym.loc[expr.index].values)
    # Collapse multiple probes per gene: keep the probe with the highest mean expression.
    expr["_mean"] = expr.drop(columns="symbol").mean(axis=1)
    expr = expr.sort_values("_mean", ascending=False).drop_duplicates("symbol").drop(columns="_mean")
    return expr.set_index("symbol")


def run(disease, cfg):
    path = download(cfg["gse"])
    g = GEOparse.get_GEO(filepath=str(path), silent=True)
    groups = {name: characteristic(s, cfg["field"]) for name, s in g.gsms.items()}
    case = [n for n, v in groups.items() if v in cfg["case"]]
    control = [n for n, v in groups.items() if v in cfg["control"]]
    expr = prepare(g, cfg)
    res = moderated_t(expr, case, control)
    res["direction"] = np.where(res.log2FC > 0, "up", "down")
    res["significant"] = (res.p_value < DEG_P) & (res.log2FC.abs() >= DEG_ABS_LFC)
    res.index.name = "symbol"
    res.to_csv(PROC / f"de_{disease}.csv")

    summary = {
        "disease": disease, "gse": cfg["gse"], "platform": next(iter(g.gpls)), "tissue": cfg["tissue"],
        "group_counts": pd.Series(groups).value_counts().to_dict(),
        "case_groups": cfg["case"], "control_groups": cfg["control"],
        "n_case": len(case), "n_control": len(control), "genes_tested": len(res),
        "n_p05": int((res.p_value < 0.05).sum()), "n_fdr05": int((res.adj_p < 0.05).sum()),
        "prior_df_d0": round(res.attrs["d0"], 2), "input_scale": cfg["scale"],
    }

    if "sensitivity_control" in cfg:
        alt = [n for n, v in groups.items() if v in cfg["sensitivity_control"]]
        alt_res = moderated_t(expr, case, alt)
        alt_res.index.name = "symbol"
        alt_res.to_csv(PROC / f"de_{disease}_vs_{'_'.join(cfg['sensitivity_control']).replace(' ', '')}.csv")
        both = res[["log2FC"]].join(alt_res[["log2FC"]], rsuffix="_alt")
        summary["sensitivity"] = {
            "control_groups": cfg["sensitivity_control"], "n_control": len(alt),
            "n_p05": int((alt_res.p_value < 0.05).sum()), "n_fdr05": int((alt_res.adj_p < 0.05).sum()),
            "log2FC_spearman_vs_primary": round(both.corr(method="spearman").iloc[0, 1], 3),
        }
        # Report BMI per group to document the obesity confound in the primary comparison.
        bmi = {}
        for n, s in g.gsms.items():
            v = characteristic(s, "bmi")
            try:
                bmi.setdefault(groups[n], []).append(float(v))
            except (TypeError, ValueError):
                pass
        summary["median_bmi_by_group"] = {k: round(float(np.median(v)), 1) for k, v in bmi.items()}
    return res, summary


def main():
    assoc = pd.read_csv(PROC / "associations.csv")
    summaries, rows = [], []
    for disease, cfg in DATASETS.items():
        res, summary = run(disease, cfg)
        summaries.append(summary)
        pairs = assoc[assoc.disease == disease].merge(res, left_on="symbol", right_index=True)
        summary["ot_genes_measured"] = len(pairs)
        summary["ot_genes_significant"] = int(pairs.significant.sum())
        rows.append(pairs[["disease", "symbol", "log2FC", "direction", "p_value", "adj_p", "significant"]])
        print(json.dumps(summary, indent=1))
    pd.concat(rows).to_csv(PROC / "expression.csv", index=False)
    (RESULTS / "geo_datasets.json").write_text(json.dumps(summaries, indent=2))


if __name__ == "__main__":
    main()
