import numpy as np
import pandas as pd
import plotly.graph_objects as go
import streamlit as st

import kb
import ui


def render():
    ui.page_header("Gene Profile", "Everything the knowledge base holds about one gene")
    t = kb.tables()
    symbols = sorted(t["genes"].symbol.dropna().unique())
    gene = st.selectbox("Gene", symbols, index=symbols.index("INSR") if "INSR" in symbols else 0,
                        key="gene_pick", help="Type to search")

    info = t["genes"][t["genes"].symbol == gene].iloc[0]
    a = t["assoc"][t["assoc"].symbol == gene].set_index("disease")
    e = t["expr"][t["expr"].symbol == gene].set_index("disease")
    pws = t["gene_pathway"][t["gene_pathway"].symbol == gene].merge(t["pathways"], on="kegg_id")

    st.markdown(f"### {gene}  \n{info['name'] or ''}")
    links = [f"[Ensembl {info.gene_id}](https://www.ensembl.org/Homo_sapiens/Gene/Summary?g={info.gene_id})",
             f"[Open Targets](https://platform.opentargets.org/target/{info.gene_id})"]
    if pd.notna(info.entrez_id) and info.entrez_id:
        links.append(f"[NCBI Gene {info.entrez_id}](https://www.ncbi.nlm.nih.gov/gene/{info.entrez_id})")
    st.markdown(" · ".join(links))

    m = st.columns(4)
    m[0].metric("Diseases associated", len(a), border=True)
    m[1].metric("Max OT score", f"{a.score.max():.2f}" if len(a) else "-", border=True)
    m[2].metric("KEGG pathways", len(pws), border=True)
    m[3].metric("Significant in", f"{int(e.significant.sum())} dataset(s)" if len(e) else "-", border=True)

    rows = []
    for d in kb.DISEASES:
        rows.append({
            "Disease": d,
            "OT score": a.score.get(d, np.nan),
            "Evidence types": a.evidence_types.get(d, "") if d in a.index else "not associated",
            "Drug-only": bool(a.drug_only.get(d, False)) if d in a.index else None,
            "GEO": f"{kb.GEO[d][0]} ({kb.GEO[d][1]})",
            "log2FC": e.log2FC.get(d, np.nan),
            "Direction": e.direction.get(d, "") if d in e.index else "",
            "p-value": e.p_value.get(d, np.nan),
            "Significant": bool(e.significant.get(d)) if d in e.index else None,
        })
    prof = pd.DataFrame(rows)[["Disease", "OT score", "log2FC", "Direction", "p-value", "Significant",
                                "Drug-only", "Evidence types", "GEO"]]

    left, right = st.columns([1.5, 1])
    with left:
        st.subheader("Disease associations and expression")
        st.dataframe(ui.styled(prof, {"OT score": "{:.3f}", "log2FC": "{:+.3f}", "p-value": "{:.2e}"},
                               keep_numeric=("OT score",)),
                     hide_index=True, width="stretch", column_config={
            "OT score": st.column_config.ProgressColumn(min_value=0, max_value=1, format="%.3f"),
        })
        st.caption("Expression values appear only where the gene is associated with that disease "
                   "(EXPRESSED_IN edges follow Open Targets pairs). Significant = moderated t p < 0.05, "
                   "uncorrected.")
    with right:
        st.subheader("log2 fold change")
        if prof.log2FC.notna().any():
            p = prof.dropna(subset=["log2FC"])
            fig = go.Figure(go.Bar(
                x=p.log2FC, y=p.Disease, orientation="h",
                marker=dict(color=[ui.DCOLOR[d] for d in p.Disease],
                            line=dict(color=["#0b0b0b" if s else "rgba(0,0,0,0)" for s in p.Significant],
                                      width=2)),
                text=[f"{v:+.2f}{' *' if s else ''}" for v, s in zip(p.log2FC, p.Significant)],
                textposition="outside", cliponaxis=False,
                hovertemplate="%{y}: log2FC %{x:.3f}<extra></extra>"))
            lim = max(1.0, float(p.log2FC.abs().max()) * 1.35)
            fig.update_layout(height=260, margin=dict(l=10, r=10, t=10, b=30), xaxis_range=[-lim, lim],
                              yaxis=dict(categoryorder="array", categoryarray=kb.DISEASES[::-1]),
                              plot_bgcolor="rgba(0,0,0,0)", paper_bgcolor="rgba(0,0,0,0)",
                              xaxis=dict(zeroline=True, zerolinecolor="#8a8984", gridcolor="#e4e3df"))
            st.plotly_chart(fig, width="stretch", config={"displayModeBar": False})
            st.caption("Outlined bar and * = significant. Positive = higher in patients.")
        else:
            ui.empty_state("No expression measurement for this gene in its associated diseases.")

    st.subheader("Neighbourhood")
    max_p = st.slider("Pathways to show", 0, 40, min(15, len(pws)), key="gene_np") if len(pws) else 0
    nodes = [{"id": f"G:{gene}", "label": gene, "color": "#0b0b0b", "size": 26, "font": 20,
              "title": info["name"] or gene}]
    edges = []
    for d in a.index:
        nodes.append(ui.disease_node(d, 28))
        fc = e.log2FC.get(d)
        tip = f"score {a.score[d]:.3f}" + (f" | log2FC {fc:+.2f}" if pd.notna(fc) else "")
        edges.append({"src": f"D:{d}", "dst": f"G:{gene}", "color": ui.DCOLOR[d], "width": 1 + 6 * a.score[d],
                      "title": tip})
    for r in pws.sort_values("size").head(max_p).itertuples():
        nodes.append(ui.pathway_node(r.kegg_id, r.pathway))
        edges.append({"src": f"G:{gene}", "dst": f"P:{r.kegg_id}", "color": "#c9c8c3", "title": "PARTICIPATES_IN"})
    ui.render_network(nodes, edges, height=460)
    st.caption("Edge width = Open Targets score. Smallest (most specific) pathways are shown first. "
               "Drag nodes to rearrange; hover for details.")

    if len(pws):
        with st.expander(f"All {len(pws)} KEGG pathways"):
            st.dataframe(pws[["kegg_id", "pathway", "category", "size"]].sort_values("kegg_id"), hide_index=True,
                         width="stretch", column_config={
                             "kegg_id": st.column_config.TextColumn("KEGG ID"),
                             "size": st.column_config.NumberColumn("Genes in pathway")})
