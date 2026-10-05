"""Shared look-and-feel: disease colours, filter controls, graph rendering, empty states."""
import pandas as pd
import streamlit as st

import kb

# Categorical slots 1-4 of a colour-blind-validated palette; the same mapping is used in results/figures.
DCOLOR = {"PCOS": "#2a78d6", "T2D": "#eb6834", "NAFLD": "#1baf7a", "Obesity": "#eda100"}
UP, DOWN, NEUTRAL = "#c0392b", "#184f95", "#8a8984"
PATHWAY_COLOR = "#b9b8b2"
DIVERGING = [[0.0, "#104281"], [0.2, "#3987e5"], [0.4, "#9ec5f4"], [0.5, "#f0efec"],
             [0.6, "#f4b0a8"], [0.8, "#e34948"], [1.0, "#9b1c1c"]]
CLASS_COLOR = {"consistent up": UP, "consistent down": DOWN, "discordant": NEUTRAL,
               "insufficient data": "#ffffff"}


def disease_chip(d):
    return (f"<span style='display:inline-block;padding:1px 8px;border-radius:10px;margin-right:4px;"
            f"background:{DCOLOR[d]};color:#0b0b0b;font-size:0.8rem;font-weight:600'>{d}</span>")


def disease_legend():
    st.markdown(" ".join(disease_chip(d) for d in kb.DISEASES), unsafe_allow_html=True)


def page_header(title, subtitle=None):
    st.title(title)
    if subtitle:
        st.caption(subtitle)


def evidence_controls(key, default_score=0.1, default_exclude=True, show_min_diseases=False):
    """The filter row shared by several pages: min OT score + drug-only toggle (+ optional min diseases)."""
    cols = st.columns([1.2, 1.2, 1] if show_min_diseases else [1.2, 1.2])
    min_score = cols[0].slider("Minimum Open Targets score", 0.1, 0.8, default_score, 0.05, key=f"{key}_score",
                               help="Associations below 0.1 are not stored in the KB.")
    exclude = cols[1].toggle("Exclude drug-only links", value=default_exclude, key=f"{key}_drug",
                             help="Drop disease-gene links supported only by drug (ChEMBL) evidence, e.g. "
                                  "metformin's mitochondrial complex I targets.")
    min_d = None
    if show_min_diseases:
        min_d = cols[2].select_slider("Minimum diseases", options=[2, 3, 4], value=2, key=f"{key}_mind")
    return min_score, exclude, min_d


def styled(df, formats=None, keep_numeric=()):
    """Display copy of df where missing numbers show as an em dash.

    st.dataframe always prints NaN as 'None', so float columns that contain gaps are rendered as formatted,
    right-aligned text (downloads keep the numeric data). Columns in keep_numeric (e.g. progress bars) stay numeric.
    """
    out = df.copy()
    fmt = {c: "{:.2f}" for c in out.select_dtypes("float").columns}
    fmt.update(formats or {})
    for c, f in fmt.items():
        if c in out.columns and c not in keep_numeric and out[c].isna().any():
            out[c] = out[c].map(lambda v: "\u2014" if pd.isna(v) else f.format(v))
    numeric_fmt = {c: f for c, f in fmt.items() if c in out.columns and out[c].dtype.kind == "f"}
    return out.style.format(numeric_fmt, na_rep="\u2014").set_properties(
        subset=[c for c in fmt if c in out.columns], **{"text-align": "right"})


def empty_state(msg):
    st.info(f"{msg}", icon=":material/search_off:")


def render_network(nodes, edges, height=640, physics=True):
    """nodes: list of dicts(id, label, color, shape, size, title); edges: list of dicts(src, dst, color, width, title)."""
    from pyvis.network import Network
    net = Network(height=f"{height}px", width="100%", bgcolor="#fcfcfb", font_color="#0b0b0b",
                  cdn_resources="in_line", directed=False)
    for n in nodes:
        net.add_node(n["id"], label=n.get("label", n["id"]), color=n.get("color"), shape=n.get("shape", "dot"),
                     size=n.get("size", 12), title=n.get("title", ""), borderWidth=n.get("border", 1),
                     font={"size": n.get("font", 14), "face": "Inter, sans-serif"})
    for e in edges:
        net.add_edge(e["src"], e["dst"], color=e.get("color", "#cccccc"), width=e.get("width", 1),
                     title=e.get("title", ""))
    net.set_options("""{
      "physics": {"enabled": %s, "solver": "forceAtlas2Based",
                  "forceAtlas2Based": {"gravitationalConstant": -60, "springLength": 110, "avoidOverlap": 0.6},
                  "stabilization": {"iterations": 250}},
      "interaction": {"hover": true, "tooltipDelay": 80, "navigationButtons": false},
      "edges": {"smooth": false}
    }""" % ("true" if physics else "false"))
    # The HTML is generated here from KB data (not user input), so embedding it with scripts is safe.
    st.iframe(net.generate_html(notebook=False), height=height + 20)


def gene_node(symbol, cls, direction, size=14, title=""):
    if cls == "consistent":
        color = UP if direction == "up" else DOWN
    elif cls == "discordant":
        color = NEUTRAL
    else:
        color = {"background": "#ffffff", "border": "#52514e"}
    return {"id": f"G:{symbol}", "label": symbol, "color": color, "shape": "dot", "size": size, "title": title}


def disease_node(d, size=34):
    return {"id": f"D:{d}", "label": d, "color": DCOLOR[d], "shape": "dot", "size": size, "font": 20,
            "title": f"{d}", "border": 2}


def pathway_node(kegg_id, name, size=10):
    return {"id": f"P:{kegg_id}", "label": name if len(name) <= 32 else name[:30] + "...",
            "color": {"background": "#ffffff", "border": "#52514e"}, "shape": "square", "size": size,
            "title": f"{name} ({kegg_id})", "font": 11}


def network_legend():
    st.markdown(
        f"<div style='font-size:0.85rem'>"
        f"{' '.join(disease_chip(d) for d in kb.DISEASES)} &nbsp; "
        f"<span style='color:{UP}'>&#9679;</span> consistently up &nbsp; "
        f"<span style='color:{DOWN}'>&#9679;</span> consistently down &nbsp; "
        f"<span style='color:{NEUTRAL}'>&#9679;</span> discordant &nbsp; "
        f"&#9675; insufficient expression data &nbsp; &#9633; KEGG pathway</div>",
        unsafe_allow_html=True)
