import streamlit as st

import kb
import ui

SCHEMA_DOT = """
digraph KB {
  rankdir=LR; bgcolor="transparent"; pad=0.2; nodesep=0.6; ranksep=1.1;
  node [shape=box, style="rounded,filled", fontname="Helvetica", fontsize=13, color="#52514e", penwidth=1.2];
  edge [fontname="Helvetica", fontsize=10, color="#52514e", fontcolor="#52514e", penwidth=1.3];
  Disease [label=<<b>Disease</b><br/><font point-size="10">id · name · source_id (MONDO/EFO)</font>>, fillcolor="#cde2fb"];
  Gene    [label=<<b>Gene</b><br/><font point-size="10">id (Ensembl) · symbol · entrez_id</font>>, fillcolor="#fde3d3"];
  Pathway [label=<<b>Pathway</b><br/><font point-size="10">kegg_id · name · category · size</font>>, fillcolor="#e4e3df"];
  Disease -> Gene [label=<<b>ASSOCIATED_WITH</b><br/>score · evidence_types · drug_only>];
  Gene -> Pathway [label=<<b>PARTICIPATES_IN</b>>];
  Gene -> Disease [label=<<b>EXPRESSED_IN</b><br/>log2FC · direction · p_value · significant>, constraint=false];
}
"""


def render():
    ui.page_header("Comorbidity-Cluster Knowledge Base",
                   "Shared genes and pathways across PCOS, type 2 diabetes, NAFLD and obesity")
    t = kb.tables()
    src = kb.backend()

    st.markdown(
        "PCOS, type 2 diabetes (T2D), non-alcoholic fatty liver disease (NAFLD) and obesity frequently occur "
        "together in the same patients. This knowledge base treats them as one **comorbidity cluster**. It "
        "links each disease to its genes (**Open Targets**), each gene to its biological pathways (**KEGG**), "
        "and records whether each gene is switched up or down in patient tissue (**GEO** expression studies). "
        "**Hub** genes and pathways are those shared by two or more diseases. **Expression-consistent** hubs "
        "also move in the same direction in the patient data. They are candidate explanations for why "
        "these conditions co-occur.")
    ui.disease_legend()
    st.write("")

    with st.spinner("Counting hubs..."):
        hubs_all = kb.hub_genes(0.1, False, src)
        hubs_robust = kb.hub_genes(0.1, True, src)
        pw = kb.pathway_summary(0.1, False, src)

    c = st.columns(5)
    c[0].metric("Genes", f"{len(t['genes']):,}", border=True,
                help="Genes associated with at least one disease (Open Targets score >= 0.1)")
    c[1].metric("Pathways", f"{len(t['pathways']):,}", border=True, help="KEGG pathways containing a KB gene")
    c[2].metric("Hub genes", f"{len(hubs_all):,}", f"{len(hubs_robust):,} excl. drug-only",
                delta_color="off", border=True, help="Associated with 2+ diseases at score >= 0.1")
    n_cons = int((hubs_all.expression_class == "consistent").sum())
    n_cons_r = int((hubs_robust.expression_class == "consistent").sum())
    c[3].metric("Consistent hubs", n_cons, f"{n_cons_r} excl. drug-only", delta_color="off",
                border=True, help="Significant (p<0.05) in 2+ associated diseases, same direction")
    c[4].metric("Hub pathways", int((pw.n_enriched >= 2).sum()), "FDR < 0.05", delta_color="off",
                border=True, help="Hypergeometric over-representation, BH FDR < 0.05")

    left, right = st.columns([1.35, 1])
    with left:
        st.subheader("Schema")
        st.graphviz_chart(SCHEMA_DOT, width="stretch")
        st.caption("Overlap between diseases is not stored. It is computed at query time from shared Gene "
                   "and Pathway nodes. Counts: 4 Disease, "
                   f"{len(t['genes']):,} Gene, {len(t['pathways']):,} Pathway nodes; "
                   f"{len(t['assoc']):,} ASSOCIATED_WITH, {len(t['gene_pathway']):,} PARTICIPATES_IN, "
                   f"{len(t['expr']):,} EXPRESSED_IN edges.")
    with right:
        st.subheader("Data sources")
        d = t["diseases"].set_index("id").reindex(kb.DISEASES)
        rows = []
        for name in kb.DISEASES:
            gse, tissue = kb.GEO[name]
            rows.append({"Disease": name, "Open Targets ID": d.loc[name, "source_id"],
                         "Genes": int((t["assoc"].disease == name).sum()), "GEO": gse, "Tissue": tissue})
        st.dataframe(rows, hide_index=True, width="stretch")
        st.markdown("Use the pages in the sidebar to explore hubs, single genes, disease overlaps, pathways, "
                    "the network, expression data, validation, and a read-only Cypher console.")
