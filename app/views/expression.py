import numpy as np
import plotly.graph_objects as go
import streamlit as st

import kb
import ui


def render():
    ui.page_header("Expression", "How the top hub genes change in patient tissue, per GEO dataset")
    min_score, exclude, min_d = ui.evidence_controls("ex", show_min_diseases=True)
    c = st.columns([1, 1.6, 1])
    n = c[0].slider("Genes", 10, 80, 40, 5, key="ex_n")
    order = c[1].radio("Rank by", ["Expression evidence first", "Association score"], horizontal=True,
                       key="ex_order")
    lim = c[2].select_slider("Colour scale limit (|log2FC|)", options=[1.0, 1.5, 2.0, 3.0], value=2.0, key="ex_lim")

    with st.spinner("Loading expression data..."):
        hubs = kb.hub_genes(min_score, exclude, kb.backend())
    hubs = hubs[hubs.n_diseases >= min_d]
    measured = hubs[[f"log2FC_{d}" for d in kb.DISEASES]].notna().sum(axis=1)
    hubs = hubs[measured >= 2]  # need at least two datasets to compare
    if order == "Association score":
        hubs = hubs.sort_values("mean_score", ascending=False)
    else:
        hubs = hubs.assign(_o=hubs.expression_class.map({"consistent": 0, "discordant": 1})).sort_values(
            ["_o", "n_deg", "mean_score"], ascending=[True, False, False], na_position="last")
    hubs = hubs.head(n)
    if hubs.empty:
        ui.empty_state("No hub genes with expression measured in two or more of their diseases.")
        return

    z = hubs[[f"log2FC_{d}" for d in kb.DISEASES]].to_numpy(float)
    sig = hubs[[f"sig_{d}" for d in kb.DISEASES]].to_numpy(bool)
    pv = hubs[[f"p_{d}" for d in kb.DISEASES]].to_numpy(float)
    text = np.where(sig & ~np.isnan(z), "*", "")
    xlabels = [f"<b>{d}</b><br>{kb.GEO[d][0]}<br>{kb.GEO[d][1]}" for d in kb.DISEASES]
    hover = np.empty(z.shape, dtype=object)
    for i, g in enumerate(hubs.symbol):
        for j, d in enumerate(kb.DISEASES):
            hover[i, j] = (f"{g} in {d}: not associated / not measured" if np.isnan(z[i, j]) else
                           f"{g} in {d}<br>log2FC {z[i, j]:+.3f}<br>p = {pv[i, j]:.2e}"
                           f"<br>{'significant' if sig[i, j] else 'not significant'}")
    fig = go.Figure(go.Heatmap(
        z=np.clip(z, -lim, lim), x=xlabels, y=hubs.symbol, colorscale=ui.DIVERGING, zmid=0, zmin=-lim, zmax=lim,
        text=text, texttemplate="%{text}", textfont=dict(size=16), xgap=2, ygap=2,
        hovertext=hover, hovertemplate="%{hovertext}<extra></extra>",
        colorbar=dict(title=dict(text="log2FC", side="right"), thickness=12, len=0.5)))
    fig.update_layout(height=max(380, 22 * len(hubs) + 140), margin=dict(l=10, r=10, t=10, b=10),
                      xaxis=dict(side="top"), yaxis=dict(autorange="reversed"),
                      plot_bgcolor="#ebeae6", paper_bgcolor="rgba(0,0,0,0)")
    left, right = st.columns([1.6, 1])
    with left:
        st.plotly_chart(fig, width="stretch", config={"displayModeBar": False})
    with right:
        st.markdown("**How to read this**")
        st.markdown(
            "- Red = higher in patients than controls, blue = lower.\n"
            "- `*` = significant (moderated t, p < 0.05, **not** corrected for multiple testing).\n"
            "- Grey cells: the gene is not associated with that disease, so no EXPRESSED_IN edge exists.\n"
            "- The four datasets come from **different tissues** and the PCOS (3 vs 7) and T2D (5 vs 5) "
            "studies are very small. Treat single cells with caution.")
        counts = hubs.expression_class.value_counts()
        st.dataframe(counts.rename_axis("Expression class").reset_index(name="Genes shown"), hide_index=True,
                     width="stretch")
        st.caption("Only hubs measured in 2+ of their associated diseases are eligible.")
