import json
import re

import pandas as pd
import streamlit as st

import kb
import ui


@st.cache_data(show_spinner=False)
def load_validation():
    md = (kb.ROOT / "VALIDATION.md").read_text()
    sections = {}
    for block in re.split(r"(?m)^## ", md)[1:]:
        title, _, body = block.partition("\n")
        sections[title.strip()] = body.strip()
    js = json.loads((kb.RESULTS / "validation.json").read_text())
    genes = pd.read_csv(kb.RESULTS / "validation_genes.csv")
    sens = pd.read_csv(kb.RESULTS / "sensitivity.csv")
    geo = json.loads((kb.RESULTS / "geo_datasets.json").read_text())
    return sections, js, genes, sens, geo


def section(sections, prefix):
    return next((v for k, v in sections.items() if k.startswith(prefix)), "")


def render():
    ui.page_header("Validation & Limitations",
                   "Does the pipeline recover known insulin-resistance genes, and how far can the results be trusted?")
    sections, js, genes, sens, geo = load_validation()
    if not sections:
        ui.empty_state("VALIDATION.md not found. Run ./run_all.sh first.")
        return

    st.subheader("INSR / IRS1 / IGF1")
    cols = st.columns(3)
    verdict = {
        "INSR": ("Hub in all 4 diseases; expression-consistent (down)", "but not significant in PCOS or T2D"),
        "IRS1": ("Hub in 3-4 diseases", "not expression-consistent (significant only in Obesity)"),
        "IGF1": ("Hub only at score >= 0.1", "literature-only links scoring 0.11-0.13; drops out at 0.2"),
    }
    for col, g in zip(cols, ["INSR", "IRS1", "IGF1"]):
        v = js["validation_genes"][g]
        deg = v["degree_by_threshold"]
        col.metric(g, f"{deg['0.1']} diseases", f"{deg['0.2']} at score >= 0.2, {deg['0.3']} at >= 0.3",
                   delta_color="off", border=True)
        col.markdown(f"**{verdict[g][0]}**, {verdict[g][1]}.")

    tbl = genes.copy()
    tbl["significant"] = tbl.significant.map({True: "yes", False: "no"})
    st.dataframe(tbl, hide_index=True, width="stretch", height=36 * (len(tbl) + 1) + 4, column_config={
        "gene": "Gene", "disease": "Disease",
        "ot_score": st.column_config.ProgressColumn("OT score", min_value=0, max_value=1, format="%.3f"),
        "evidence_types": "Evidence types", "measured_in_geo": "Measured",
        "log2FC": st.column_config.NumberColumn("log2FC", format="%+.3f"),
        "p_value": st.column_config.NumberColumn("p", format="%.2e"), "significant": "Significant"})

    st.subheader("Bottom line")
    st.markdown(section(sections, "Bottom line"))

    st.subheader("Caveats, stated plainly")
    bg = js["expression_consistency_background"]
    sizes = {g["disease"]: (g["n_case"], g["n_control"]) for g in geo}
    no_drug = sens[sens.exclude_drug_only & (sens.min_score == 0.1)].iloc[0]
    primary = sens[~sens.exclude_drug_only & (sens.min_score == 0.1)].iloc[0]
    st.warning(
        f"**Drug-only (metformin) effect.** Many of the most widely shared genes are linked to the diseases "
        f"only through drug evidence. Metformin acts on mitochondrial complex I and is used or trialled in all "
        f"four conditions, so Open Targets connects its target genes (NDUF\\*, MT-ND\\*) to every one of them. "
        f"Removing drug-only links cuts the genes shared by all four diseases from "
        f"{int(primary.hub_genes_4_diseases)} to {int(no_drug.hub_genes_4_diseases)}.",
        icon=":material/medication:")
    st.warning(
        f"**Same-direction agreement is not specific to hubs.** Across *all* genes significant in two or more "
        f"datasets, {bg['background_concordance']:.1%} already change in the same direction. Hub genes reach "
        f"{bg['hub_concordance']:.0%}, which is not significantly higher (p = {bg['binomial_p_hub_vs_background']:.2f}). "
        f"Much of the consistency reflects a shared genome-wide disease signature, not hub-specific biology.",
        icon=":material/compare_arrows:")
    st.warning(
        f"**Small PCOS and T2D studies.** PCOS: {sizes['PCOS'][0]} patients vs {sizes['PCOS'][1]} controls. "
        f"T2D: {sizes['T2D'][0]} vs {sizes['T2D'][1]}. With so few samples, real changes are easily missed. "
        f"None of INSR, IRS1 or IGF1 is significant in both the PCOS and the T2D dataset.",
        icon=":material/group:")
    st.warning(
        "**Uncorrected p-value cutoff.** A gene counts as changed at p < 0.05 *without* correction for "
        "testing about 20,000 genes. This lenient cutoff is what lets the small PCOS and T2D datasets contribute "
        "at all (with FDR < 0.05 they yield 3 and 0 genes), but it also lets through false positives.",
        icon=":material/functions:")
    st.info("Other limitations: the four datasets come from different tissues (ovary, two fat depots, liver); "
            "T2D and obesity have far more known genes than PCOS; association is not causation; the NAFLD "
            "comparison is partly an obesity comparison (patients median BMI ~48 vs controls ~26).",
            icon=":material/info:")

    st.subheader("Supporting analyses (from VALIDATION.md)")
    for prefix in ["1.", "2.", "3.", "4.", "5."]:
        title = next((k for k in sections if k.startswith(prefix)), None)
        if title:
            with st.expander(title, expanded=prefix == "4."):
                st.markdown(sections[title])
