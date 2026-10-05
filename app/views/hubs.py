import streamlit as st

import kb
import ui


def render():
    ui.page_header("Hub Explorer", "Genes associated with two or more diseases in the cluster")
    min_score, exclude, min_d = ui.evidence_controls("hubs", show_min_diseases=True)
    c1, c2 = st.columns([1, 2.4])
    consistent_only = c1.toggle("Expression-consistent only", value=False, key="hubs_cons")
    search = c2.text_input("Search gene symbol", "", key="hubs_search", placeholder="e.g. INSR, PPARG, IL",
                           label_visibility="collapsed").strip().upper()

    with st.spinner("Querying hubs..."):
        hubs = kb.hub_genes(min_score, exclude, kb.backend())
    view = hubs[hubs.n_diseases >= min_d]
    if consistent_only:
        view = view[view.expression_class == "consistent"]
    if search:
        view = view[view.symbol.str.upper().str.contains(search, regex=False)]

    m = st.columns(4)
    m[0].metric("Hub genes shown", f"{len(view):,}", border=True)
    m[1].metric("In all 4 diseases", int((view.n_diseases == 4).sum()), border=True)
    m[2].metric("Expression-consistent", int((view.expression_class == "consistent").sum()), border=True)
    m[3].metric("Discordant", int((view.expression_class == "discordant").sum()), border=True)

    if view.empty:
        ui.empty_state("No hub genes match these filters. Try lowering the score threshold or clearing the search.")
        return

    table = view[["symbol", "name", "n_diseases", "diseases", "mean_score", "expression_class", "direction",
                  "drug_only_links"] + [f"log2FC_{d}" for d in kb.DISEASES]].copy()
    table["direction"] = table.direction.fillna("")
    table = table.rename(columns={f"log2FC_{d}": f"log2FC {d}" for d in kb.DISEASES})
    st.dataframe(
        ui.styled(table, {f"log2FC {d}": "{:+.2f}" for d in kb.DISEASES}), hide_index=True, width="stretch", height=520,
        column_config={
            "symbol": st.column_config.TextColumn("Gene", width="small"),
            "name": st.column_config.TextColumn("Name", width="medium"),
            "n_diseases": st.column_config.NumberColumn("# diseases", width="small"),
            "diseases": "Associated with",
            "mean_score": st.column_config.ProgressColumn("Mean OT score", min_value=0, max_value=1, format="%.2f"),
            "expression_class": "Expression",
            "direction": "Direction",
            "drug_only_links": st.column_config.NumberColumn("Drug-only links", width="small",
                                                             help="Links supported only by drug evidence"),

        })
    st.download_button("Download as CSV", view.to_csv(index=False).encode(), file_name=(
        f"hub_genes_score{min_score:.2f}_{'nodrug' if exclude else 'all'}.csv"), mime="text/csv",
        icon=":material/download:")
    st.caption("Expression = consistent: significant (p<0.05) in 2+ associated diseases, same direction. "
               "Discordant: significant in 2+ but opposite directions. Open a gene on the Gene Profile page "
               "for details.")
