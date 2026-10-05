"""Comorbidity-Cluster KB - Streamlit frontend.

Launch from the repository root:  streamlit run app/app.py
"""
import streamlit as st

st.set_page_config(page_title="Comorbidity KB", page_icon=":material/hub:", layout="wide")

import kb  # noqa: E402  (after set_page_config)
import ui  # noqa: E402
from views import (console, expression, gene, home, hubs, network, overlap, pathways,  # noqa: E402
                   validation)

pages = {
    "Overview": [
        st.Page(home.render, title="Home", icon=":material/home:", url_path="home", default=True),
    ],
    "Explore": [
        st.Page(hubs.render, title="Hub Explorer", icon=":material/table_view:", url_path="hubs"),
        st.Page(gene.render, title="Gene Profile", icon=":material/genetics:", url_path="gene"),
        st.Page(overlap.render, title="Disease Overlap", icon=":material/join_inner:", url_path="overlap"),
        st.Page(pathways.render, title="Pathways", icon=":material/route:", url_path="pathways"),
        st.Page(network.render, title="Network", icon=":material/hub:", url_path="network"),
        st.Page(expression.render, title="Expression", icon=":material/grid_on:", url_path="expression"),
    ],
    "Trust & tools": [
        st.Page(validation.render, title="Validation & Limitations", icon=":material/verified:",
                url_path="validation"),
        st.Page(console.render, title="Query Console", icon=":material/terminal:", url_path="console"),
    ],
}

nav = st.navigation(pages)

with st.sidebar:
    status = kb.neo4j_status()
    if status == "connected":
        st.caption(f":material/database: Connected to Neo4j at `{kb.NEO4J_URI}`")
    else:
        st.info("**Running in demo mode on the bundled database.**", icon=":material/inventory_2:")
        reason = ("Neo4j is configured but unreachable." if status == "unreachable"
                  else "No Neo4j connection is configured.")
        st.caption(f"{reason} All pages read a read-only SQLite snapshot of the knowledge base; "
                   "live Cypher is available when the app runs locally with Neo4j.")
    ui.disease_legend()
    st.caption("Data: Open Targets, KEGG, GEO (snapshot 2026-10-05)")

if not kb.SQLITE_PATH.exists() and kb.backend() != "neo4j":
    st.error("No knowledge base found. Run `./run_all.sh` first to build it.")
    st.stop()

nav.run()
