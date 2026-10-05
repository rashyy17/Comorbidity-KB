import numpy as np
import plotly.graph_objects as go
import streamlit as st

import kb
import ui


def render():
    ui.page_header("Pathways", "KEGG pathways over-represented among each disease's genes")
    min_score, exclude, _ = ui.evidence_controls("pw", default_exclude=False)
    c = st.columns([1, 1.4, 1.6])
    min_en = c[0].select_slider("Enriched in at least", options=[1, 2, 3, 4], value=2, key="pw_min",
                                format_func=lambda x: f"{x} disease{'s' if x > 1 else ''}")
    cats = c[1].multiselect("KEGG category", ["pathway", "disease"], default=["pathway", "disease"], key="pw_cat",
                            help="'disease' = KEGG human-disease maps (hsa05xxx). Global overview maps are "
                                 "excluded from testing.")
    search = c[2].text_input("Search pathway name", "", key="pw_q", placeholder="e.g. insulin").strip().lower()

    with st.spinner("Running over-representation tests..."):
        summ = kb.pathway_summary(min_score, exclude, kb.backend())
        enr = kb.pathway_enrichment(min_score, exclude, kb.backend())
    view = summ[(summ.n_enriched >= min_en) & summ.category.isin(cats)]
    if search:
        view = view[view.pathway.str.lower().str.contains(search, regex=False)]

    st.caption(f"{len(view)} pathways shown. A pathway counts as linked to a disease when that disease's genes "
               "are over-represented in it (hypergeometric test, Benjamini-Hochberg FDR < 0.05 within each "
               "disease). Without this test, almost every pathway would count as shared.")
    if view.empty:
        ui.empty_state("No pathways match these filters.")
    else:
        cols = ["kegg_id", "pathway", "category", "size", "n_enriched", "enriched_in", "min_fdr"] + \
               [f"k_{d}" for d in kb.DISEASES] + [f"fdr_{d}" for d in kb.DISEASES]
        st.dataframe(view[cols], hide_index=True, width="stretch", height=420, column_config={
            "kegg_id": st.column_config.TextColumn("KEGG ID", width="small"),
            "pathway": st.column_config.TextColumn("Pathway", width="medium"),
            "size": st.column_config.NumberColumn("Genes in pathway"),
            "n_enriched": st.column_config.NumberColumn("# diseases enriched"),
            "enriched_in": "Enriched in",
            "min_fdr": st.column_config.NumberColumn("Best FDR", format="%.1e"),
            **{f"k_{d}": st.column_config.NumberColumn(f"{d} genes", help=f"{d} genes in the pathway")
               for d in kb.DISEASES},
            **{f"fdr_{d}": st.column_config.NumberColumn(f"{d} FDR", format="%.1e") for d in kb.DISEASES},
        })
        st.download_button("Download as CSV", view.to_csv(index=False).encode(), "hub_pathways.csv", "text/csv",
                           icon=":material/download:")

    st.divider()
    st.subheader("Pathway to diseases")
    opts = (view if not view.empty else summ)
    labels = {f"{r.pathway} ({r.kegg_id})": r.kegg_id for r in opts.itertuples()}
    default = next((k for k, v in labels.items() if v == "hsa04931"), next(iter(labels)))
    pick = st.selectbox("Pathway", list(labels), index=list(labels).index(default), key="pw_pick")
    kid = labels[pick]
    sub = enr[enr.kegg_id == kid].set_index("disease").reindex(kb.DISEASES)

    left, right = st.columns([1.1, 1])
    with left:
        y = -np.log10(sub.fdr.clip(lower=1e-300))
        fig = go.Figure(go.Bar(
            x=kb.DISEASES, y=y, marker_color=[ui.DCOLOR[d] for d in kb.DISEASES],
            text=[f"{int(k)} genes<br>{fe:.1f}x" if k else "0" for k, fe in
                  zip(sub.k, sub.fold_enrichment.fillna(0))],
            textposition="outside", cliponaxis=False,
            customdata=np.stack([sub.k, sub.expected, sub.fdr], axis=1),
            hovertemplate="%{x}<br>%{customdata[0]} genes (expected %{customdata[1]:.1f})"
                          "<br>FDR %{customdata[2]:.2e}<extra></extra>"))
        fig.add_hline(y=-np.log10(kb.PATHWAY_FDR), line_dash="dot", line_color="#52514e",
                      annotation_text="FDR 0.05", annotation_position="top right")
        fig.update_layout(height=340, margin=dict(l=10, r=10, t=30, b=10), yaxis_title="-log10 FDR",
                          plot_bgcolor="rgba(0,0,0,0)", paper_bgcolor="rgba(0,0,0,0)",
                          yaxis=dict(gridcolor="#e4e3df", range=[0, max(2.0, float(y.max()) * 1.25)]))
        st.plotly_chart(fig, width="stretch", config={"displayModeBar": False})
        st.caption("Bar height = significance of over-representation; label = genes in pathway and fold "
                   "enrichment over chance.")
    with right:
        t = kb.tables()
        members = set(t["gene_pathway"][t["gene_pathway"].kegg_id == kid].symbol)
        a = kb.filtered_assoc(min_score, exclude)
        a = a[a.symbol.isin(members)]
        if a.empty:
            ui.empty_state("No disease genes in this pathway at this threshold.")
        else:
            tbl = a.pivot_table(index="symbol", columns="disease", values="score").reindex(columns=kb.DISEASES)
            tbl.insert(0, "# diseases", tbl.notna().sum(axis=1))
            tbl = tbl.sort_values(["# diseases"], ascending=False).reset_index()
            st.markdown(f"**Disease genes in {pick}** ({len(tbl)} of {len(members)} KB genes in the pathway)")
            st.dataframe(ui.styled(tbl), hide_index=True, width="stretch", height=300, column_config={
                "symbol": "Gene"})
        st.markdown(f"[Open {kid} in KEGG](https://www.kegg.jp/pathway/{kid})")
