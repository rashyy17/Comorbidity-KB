from itertools import combinations

import pandas as pd
import plotly.graph_objects as go
import streamlit as st
from plotly.subplots import make_subplots

import kb
import ui


def upset_figure(sets, chosen):
    """UpSet-style chart: exclusive intersection sizes (bars) over a dot matrix of which diseases take part."""
    combos = []
    for r in range(len(chosen), 0, -1):
        for c in combinations(chosen, r):
            inside = set.intersection(*(sets[d] for d in c))
            outside = set().union(*(sets[d] for d in chosen if d not in c))
            combos.append((c, len(inside - outside)))
    combos = [c for c in combos if c[1] > 0]
    combos.sort(key=lambda x: (-len(x[0]), -x[1]))
    if not combos:
        return None, combos

    labels = [" ∩ ".join(c) for c, _ in combos]
    fig = make_subplots(rows=2, cols=1, shared_xaxes=True, row_heights=[0.68, 0.32], vertical_spacing=0.03)
    colors = ["#52514e" if len(c) > 1 else ui.DCOLOR[c[0]] for c, _ in combos]
    fig.add_bar(x=labels, y=[n for _, n in combos], marker_color=colors, text=[n for _, n in combos],
                textposition="outside", cliponaxis=False, row=1, col=1,
                hovertemplate="%{x}<br>%{y} genes (exclusive)<extra></extra>")
    for i, d in enumerate(chosen[::-1]):
        member = [d in c for c, _ in combos]
        fig.add_scatter(x=labels, y=[d] * len(labels), mode="markers", row=2, col=1, showlegend=False,
                        marker=dict(size=13, color=[ui.DCOLOR[d] if m else "#e4e3df" for m in member],
                                    line=dict(width=0)), hoverinfo="skip")
    for c, _ in combos:  # connector lines between member dots
        if len(c) > 1:
            fig.add_scatter(x=[" ∩ ".join(c)] * len(c), y=list(c), mode="lines", row=2, col=1,
                            line=dict(color="#52514e", width=2), showlegend=False, hoverinfo="skip")
    fig.update_layout(height=470, margin=dict(l=10, r=10, t=20, b=10), showlegend=False,
                      plot_bgcolor="rgba(0,0,0,0)", paper_bgcolor="rgba(0,0,0,0)", bargap=0.35)
    fig.update_xaxes(showticklabels=False, row=1, col=1)
    fig.update_xaxes(showticklabels=False, row=2, col=1)
    fig.update_yaxes(title_text="genes", gridcolor="#e4e3df", row=1, col=1)
    fig.update_yaxes(categoryorder="array", categoryarray=chosen[::-1], row=2, col=1)
    # draw member dots above connector lines
    fig.data = tuple(sorted(fig.data, key=lambda tr: 0 if getattr(tr, "mode", None) == "lines" else 1))
    return fig, combos


def render():
    ui.page_header("Disease Overlap", "Which genes do the selected diseases share?")
    chosen = st.multiselect("Diseases (2-4)", kb.DISEASES, default=kb.DISEASES, key="ov_dis",
                            max_selections=4)
    min_score, exclude, _ = ui.evidence_controls("ov")
    if len(chosen) < 2:
        ui.empty_state("Pick at least two diseases.")
        return
    chosen = [d for d in kb.DISEASES if d in chosen]
    sets = kb.disease_gene_sets(min_score, exclude)

    m = st.columns(len(chosen) + 1)
    for i, d in enumerate(chosen):
        m[i].metric(f"{d} genes", f"{len(sets[d]):,}", border=True)
    shared_all = set.intersection(*(sets[d] for d in chosen))
    m[-1].metric("Shared by all selected", f"{len(shared_all):,}", border=True)

    with st.spinner("Computing intersections..."):
        fig, combos = upset_figure(sets, chosen)
    st.subheader("Intersections")
    if fig is None:
        ui.empty_state("No genes at this threshold.")
        return
    st.plotly_chart(fig, width="stretch", config={"displayModeBar": False})
    st.caption("UpSet-style chart: each bar counts genes found in exactly that combination of diseases "
               "(dots below) and none of the other selected ones. Coloured bars are disease-specific genes.")

    pair_rows = [{"Pair": f"{a} ∩ {b}", "Shared genes": len(sets[a] & sets[b]),
                  "Jaccard": len(sets[a] & sets[b]) / max(1, len(sets[a] | sets[b]))}
                 for a, b in combinations(chosen, 2)]

    left, right = st.columns([2, 1])
    with right:
        st.subheader("Pairwise")
        st.dataframe(pd.DataFrame(pair_rows), hide_index=True, width="stretch",
                     column_config={"Jaccard": st.column_config.NumberColumn(format="%.3f")})
        st.caption("Jaccard = shared / union. T2D and obesity have far more known genes, so they share "
                   "more almost by default.")
    with left:
        st.subheader("Shared genes")
        options = [" ∩ ".join(c) for c, _ in combos if len(c) > 1]
        mode = st.radio("Show", ["Shared by all selected", "A specific intersection"], horizontal=True,
                        key="ov_mode")
        if mode == "Shared by all selected":
            genes = shared_all
        else:
            if not options:
                ui.empty_state("No multi-disease intersections.")
                return
            pick = st.selectbox("Intersection", options, key="ov_pick")
            c = pick.split(" ∩ ")
            genes = set.intersection(*(sets[d] for d in c)) - set().union(
                *(sets[d] for d in chosen if d not in c))
        if not genes:
            ui.empty_state("No genes in this intersection.")
            return
        a = kb.filtered_assoc(min_score, exclude)
        tbl = a[a.symbol.isin(genes) & a.disease.isin(chosen)].pivot_table(
            index="symbol", columns="disease", values="score").reindex(columns=chosen)
        tbl["mean_score"] = tbl.mean(axis=1)
        tbl = tbl.sort_values("mean_score", ascending=False).reset_index()
        st.dataframe(ui.styled(tbl, {d: "{:.3f}" for d in chosen}), hide_index=True, width="stretch", height=380, column_config={
            "symbol": "Gene",
            "mean_score": st.column_config.ProgressColumn("Mean score", min_value=0, max_value=1, format="%.2f")})
        st.download_button("Download as CSV", tbl.to_csv(index=False).encode(), "shared_genes.csv", "text/csv",
                           icon=":material/download:")
