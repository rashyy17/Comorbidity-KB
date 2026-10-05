import streamlit as st

import kb
import ui


def render():
    ui.page_header("Network", "Top hub genes, the diseases they connect, and the pathways they share")
    min_score, exclude, min_d = ui.evidence_controls("net", show_min_diseases=True)
    c = st.columns(3)
    n_hubs = c[0].slider("Hub genes to show", 5, 80, 30, 5, key="net_n")
    n_pw = c[1].slider("Shared pathways to show", 0, 25, 8, key="net_p",
                       help="Pathways enriched in 3+ diseases that contain the most displayed hubs")
    physics = c[2].toggle("Live layout (physics)", value=True, key="net_phys")

    with st.spinner("Building network..."):
        hubs = kb.hub_genes(min_score, exclude, kb.backend())
        hubs = hubs[hubs.n_diseases >= min_d].head(n_hubs)
        if hubs.empty:
            ui.empty_state("No hub genes at these settings.")
            return
        a = kb.filtered_assoc(min_score, exclude)
        a = a[a.symbol.isin(hubs.symbol)]
        t = kb.tables()
        summ = kb.pathway_summary(min_score, exclude, kb.backend())
        cand = summ[(summ.n_enriched >= 3) & (summ.category == "pathway")]
        gp = t["gene_pathway"][t["gene_pathway"].symbol.isin(hubs.symbol) & t["gene_pathway"].kegg_id.isin(cand.kegg_id)]
        top_pw = gp.kegg_id.value_counts().head(n_pw).index
        gp = gp[gp.kegg_id.isin(top_pw)]
        names = t["pathways"].set_index("kegg_id").pathway

        nodes = [ui.disease_node(d) for d in kb.DISEASES]
        for r in hubs.itertuples():
            tip = (f"{r.symbol}: {r.diseases}\nmean score {r.mean_score:.2f}\nexpression: {r.expression_class}"
                   + (f" ({r.direction})" if isinstance(r.direction, str) else ""))
            nodes.append(ui.gene_node(r.symbol, r.expression_class, r.direction, size=8 + 4 * r.n_diseases,
                                      title=tip))
        for k in top_pw:
            nodes.append(ui.pathway_node(k, names.get(k, k), size=12))
        edges = [{"src": f"D:{r.disease}", "dst": f"G:{r.symbol}", "color": ui.DCOLOR[r.disease],
                  "width": 1 + 4 * r.score, "title": f"{r.disease}-{r.symbol}: score {r.score:.2f}"}
                 for r in a.itertuples()]
        edges += [{"src": f"G:{r.symbol}", "dst": f"P:{r.kegg_id}", "color": "#d5d4cf", "width": 0.8}
                  for r in gp.itertuples()]

    ui.network_legend()
    ui.render_network(nodes, edges, height=680, physics=physics)
    st.caption(f"{len(hubs)} hub genes, {len(top_pw)} pathways, {len(edges)} edges. Node size = number of "
               "diseases; edge width = Open Targets score. Drag to rearrange, scroll to zoom, hover for details.")
