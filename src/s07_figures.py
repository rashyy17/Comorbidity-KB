"""Step 7: figures.

fig1_overlap_graph.png   disease-gene-pathway overlap graph (top robust hubs + their hub pathways)
fig2_hub_gene_table.png  hub gene table: association dot-matrix + expression direction + evidence
fig3_expression_heatmap.png  log2FC heatmap of top hub genes across the four GEO datasets
"""
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.lines import Line2D

from config import FIGURES, PROC, RESULTS

DISEASES = ["PCOS", "T2D", "NAFLD", "Obesity"]
# Categorical slots 1-4 of the validated reference palette (fixed order, identity = disease).
DCOLOR = {"PCOS": "#2a78d6", "T2D": "#eb6834", "NAFLD": "#1baf7a", "Obesity": "#eda100"}
INK, INK2, MUTED, GRID, SURFACE = "#0b0b0b", "#52514e", "#8a8984", "#e4e3df", "#fcfcfb"
UP, DOWN, MID = "#c0392b", "#184f95", "#f0efec"
DIVERGING = LinearSegmentedColormap.from_list(
    "bluered", ["#104281", "#3987e5", "#9ec5f4", MID, "#f4b0a8", "#e34948", "#9b1c1c"])

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 9, "text.color": INK, "axes.labelcolor": INK2,
    "xtick.color": INK2, "ytick.color": INK2, "axes.edgecolor": GRID, "figure.facecolor": SURFACE,
    "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE, "savefig.dpi": 200,
})


def robust_hubs():
    """Hubs that stay hubs when drug-only (ChEMBL) edges are removed, ranked for display."""
    h = pd.read_csv(RESULTS / "hub_genes.csv")
    h = h[h.hub_without_drug_only].copy()
    order = {"consistent": 0, "discordant": 1, "insufficient": 2}
    h["_o"] = h.expression_class.map(order)
    return h.sort_values(["n_diseases", "_o", "mean_score"], ascending=[False, True, False])


def fig1_graph(n_genes=24, n_pathways=10):
    hubs = robust_hubs().head(n_genes)
    assoc = pd.read_csv(PROC / "associations.csv")
    assoc = assoc[(assoc.score >= 0.1) & ~assoc.drug_only & assoc.symbol.isin(hubs.symbol)]
    gp = pd.read_csv(PROC / "gene_pathway.csv")
    hp = pd.read_csv(RESULTS / "hub_pathways_excluding_drug_only.csv")
    hp = hp[(hp.category == "pathway") & (hp.n_enriched >= 3)]
    links = gp[gp.symbol.isin(hubs.symbol) & gp.kegg_id.isin(hp.kegg_id)]
    top_p = links.kegg_id.value_counts().head(n_pathways).index
    links = links[links.kegg_id.isin(top_p)]
    pname = hp.set_index("kegg_id").pathway

    G = nx.Graph()
    for d in DISEASES:
        G.add_node(d, kind="disease")
    for g in hubs.symbol:
        G.add_node(g, kind="gene")
    for p in top_p:
        G.add_node(p, kind="pathway")
    for r in assoc.itertuples():
        G.add_edge(r.disease, r.symbol, kind="assoc", disease=r.disease)
    for r in links.itertuples():
        G.add_edge(r.symbol, r.kegg_id, kind="pathway")

    # Diseases pinned on an inner square, pathways on an outer ring, genes placed by spring layout.
    fixed = {"PCOS": (-1, 1), "T2D": (1, 1), "NAFLD": (-1, -1), "Obesity": (1, -1)}
    for i, p in enumerate(top_p):
        a = 2 * np.pi * i / len(top_p) + np.pi / len(top_p)
        fixed[p] = (2.55 * np.cos(a), 2.55 * np.sin(a))
    pos = nx.spring_layout(G, pos=fixed, fixed=list(fixed), seed=7, k=0.9, iterations=500)

    fig, ax = plt.subplots(figsize=(11, 10))
    pe = [e for e in G.edges(data=True) if e[2]["kind"] == "pathway"]
    nx.draw_networkx_edges(G, pos, edgelist=pe, edge_color="#c9c8c3", width=0.6, ax=ax)
    for d in DISEASES:
        de = [e for e in G.edges(data=True) if e[2].get("disease") == d]
        nx.draw_networkx_edges(G, pos, edgelist=de, edge_color=DCOLOR[d], width=1.4, alpha=0.75, ax=ax)

    cls = hubs.set_index("symbol").expression_class
    direction = pd.read_csv(RESULTS / "expression_consistent_hubs.csv").set_index("symbol").direction
    gfill = {g: (UP if direction.get(g) == "up" else DOWN) if cls[g] == "consistent"
             else ("#8a8984" if cls[g] == "discordant" else "#ffffff") for g in hubs.symbol}
    nx.draw_networkx_nodes(G, pos, nodelist=list(hubs.symbol), node_color=[gfill[g] for g in hubs.symbol],
                           edgecolors=INK2, linewidths=1.0, node_size=[90 + 70 * G.degree(g) for g in hubs.symbol],
                           ax=ax)
    nx.draw_networkx_nodes(G, pos, nodelist=list(top_p), node_shape="s", node_color="#ffffff",
                           edgecolors=INK2, linewidths=1.0, node_size=160, ax=ax)
    nx.draw_networkx_nodes(G, pos, nodelist=DISEASES, node_color=[DCOLOR[d] for d in DISEASES],
                           edgecolors=SURFACE, linewidths=2, node_size=2300, ax=ax)
    for d in DISEASES:
        ax.text(*pos[d], d, ha="center", va="center", fontsize=11, fontweight="bold", color=INK)
    for g in hubs.symbol:
        x, y = pos[g]
        ax.text(x, y + 0.085, g, ha="center", va="bottom", fontsize=7.5, color=INK,
                bbox=dict(boxstyle="round,pad=0.12", fc=SURFACE, ec="none", alpha=0.85))
    for p in top_p:
        x, y = pos[p]
        label = pname.get(p, p)
        ha = "left" if x > 0.3 else ("right" if x < -0.3 else "center")
        dx = 0.12 if ha == "left" else (-0.12 if ha == "right" else 0)
        ax.text(x + dx, y + (0 if ha != "center" else (0.14 if y > 0 else -0.14)), f"{label}\n({p})",
                ha=ha, va="center", fontsize=7.5, color=INK2)

    legend = [
        Line2D([], [], marker="o", ls="", mfc=UP, mec=INK2, ms=8, label="Hub gene - consistently UP (GEO)"),
        Line2D([], [], marker="o", ls="", mfc=DOWN, mec=INK2, ms=8, label="Hub gene - consistently DOWN (GEO)"),
        Line2D([], [], marker="o", ls="", mfc="#8a8984", mec=INK2, ms=8, label="Hub gene - discordant direction"),
        Line2D([], [], marker="o", ls="", mfc="#ffffff", mec=INK2, ms=8, label="Hub gene - too little GEO signal"),
        Line2D([], [], marker="s", ls="", mfc="#ffffff", mec=INK2, ms=8, label="KEGG hub pathway (enriched in 3+ diseases)"),
    ] + [Line2D([], [], color=DCOLOR[d], lw=2, label=f"ASSOCIATED_WITH {d}") for d in DISEASES]
    ax.legend(handles=legend, loc="upper left", bbox_to_anchor=(0, -0.01), ncol=3, frameon=False, fontsize=8)
    ax.set_title(f"Comorbidity overlap graph: top {len(hubs)} hub genes shared by 3-4 diseases\n"
                 "and the KEGG pathways they share (drug-only Open Targets edges excluded)",
                 loc="left", fontsize=12, color=INK)
    ax.axis("off")
    ax.set_xlim(-3.6, 3.6)
    ax.set_ylim(-3.0, 3.0)
    fig.tight_layout()
    fig.savefig(FIGURES / "fig1_overlap_graph.png", bbox_inches="tight")
    plt.close(fig)


def fig2_table(n=30):
    hubs = robust_hubs().head(n).reset_index(drop=True)
    assoc = pd.read_csv(PROC / "associations.csv")
    score = assoc.pivot(index="symbol", columns="disease", values="score")
    evid = assoc.pivot(index="symbol", columns="disease", values="evidence_types")

    fig, ax = plt.subplots(figsize=(10.5, 0.32 * n + 1.6))
    cols = {"gene": 0.0, **{d: 1.55 + i * 0.95 for i, d in enumerate(DISEASES)}, "n": 5.5, "mean": 6.3,
            "expr": 7.25}
    y0 = n
    hdr = dict(fontsize=8.5, color=INK2, fontweight="bold", va="center")
    ax.text(cols["gene"], y0, "Gene", ha="left", **hdr)
    for d in DISEASES:
        ax.text(cols[d], y0, d, ha="center", **hdr)
    ax.text(cols["n"], y0, "# dis.", ha="center", **hdr)
    ax.text(cols["mean"], y0, "Mean OT\nscore", ha="center", **hdr)
    ax.text(cols["expr"], y0, "GEO direction (significant datasets)", ha="left", **hdr)
    ax.axhline(y0 - 0.55, color=INK2, lw=0.8)

    for i, r in hubs.iterrows():
        y = n - 1 - i
        if i % 2 == 0:
            ax.axhspan(y - 0.5, y + 0.5, color="#f3f2ee", zorder=0)
        ax.text(cols["gene"], y, r.symbol, ha="left", va="center", fontsize=9, fontweight="bold")
        for d in DISEASES:
            s = score.loc[r.symbol].get(d)
            if pd.notna(s) and s >= 0.1:
                drug = str(evid.loc[r.symbol].get(d)) == "clinical"
                size = 40 + 260 * min(s, 1)
                ax.scatter(cols[d], y, s=size, color=DCOLOR[d] if not drug else SURFACE, edgecolors=DCOLOR[d],
                           linewidths=1.2, zorder=3)
                sig = r.get(f"sig_{d}")
                if sig is True or sig == "True":
                    arrow = "▲" if r[f"log2FC_{d}"] > 0 else "▼"
                    ax.text(cols[d] + 0.27, y, arrow, ha="left", va="center", fontsize=7.5,
                            color=UP if r[f"log2FC_{d}"] > 0 else DOWN)
            else:
                ax.text(cols[d], y, "·", ha="center", va="center", color=MUTED, fontsize=12)
        ax.text(cols["n"], y, str(r.n_diseases), ha="center", va="center")
        ax.text(cols["mean"], y, f"{r.mean_score:.2f}", ha="center", va="center", color=INK2)
        label = {"consistent": "consistent", "discordant": "discordant", "insufficient": "-"}[r.expression_class]
        if r.expression_class == "consistent":
            ups = [d for d in DISEASES if (r.get(f"sig_{d}") in (True, "True")) and r[f"log2FC_{d}"] > 0
                   and pd.notna(score.loc[r.symbol].get(d))]
            label = f"consistent {'UP' if ups else 'DOWN'}"
        ax.text(cols["expr"], y, label, ha="left", va="center",
                color=(UP if "UP" in label else DOWN) if "consistent" in label else INK2,
                fontweight="bold" if "consistent" in label else "normal")

    ax.set_xlim(-0.2, 10.2)
    ax.set_ylim(-1.6, n + 0.7)
    ax.axis("off")
    ax.text(-0.2, -1.05, "Dot = Open Targets association (area ~ score; hollow = drug-only evidence).  "
            "▲/▼ = significant up/down-regulation (moderated t, p<0.05) in that disease's GEO dataset.\n"
            "Ranked: number of diseases, then expression-consistent first, then mean association score. "
            "Full list: results/hub_genes.csv", fontsize=7.5, color=INK2, va="top")
    ax.set_title(f"Top {n} hub genes (associated with 2+ diseases; robust to excluding drug-only edges)",
                 loc="left", fontsize=12)
    fig.tight_layout()
    fig.savefig(FIGURES / "fig2_hub_gene_table.png", bbox_inches="tight")
    plt.close(fig)


def fig3_heatmap(n=40):
    hubs = robust_hubs()
    # Keep genes measured on at least 3 of the 4 arrays (mitochondrially encoded MT-* genes are not).
    full_de = pd.concat({d: pd.read_csv(PROC / f"de_{d}.csv", index_col=0).log2FC for d in DISEASES}, axis=1)
    measured = full_de.reindex(hubs.symbol).notna().sum(axis=1).to_numpy()
    hubs = hubs[measured >= 3].head(n)
    lfc = hubs.set_index("symbol")[[f"log2FC_{d}" for d in DISEASES]]
    lfc.columns = DISEASES
    sig = hubs.set_index("symbol")[[f"sig_{d}" for d in DISEASES]].isin([True, "True"])
    sig.columns = DISEASES
    assoc = pd.read_csv(PROC / "associations.csv")
    linked = assoc[assoc.score >= 0.1].pivot(index="symbol", columns="disease", values="score").reindex(
        index=lfc.index, columns=DISEASES).notna()

    # Show the measured log2FC from the full DE table even where the gene is not associated (hatched).
    full = pd.concat({d: pd.read_csv(PROC / f"de_{d}.csv", index_col=0).log2FC for d in DISEASES}, axis=1)
    lfc = full.reindex(index=lfc.index, columns=DISEASES)

    lim = 2.0
    fig, ax = plt.subplots(figsize=(6.2, 0.26 * n + 1.8))
    im = ax.imshow(lfc.clip(-lim, lim).to_numpy(), cmap=DIVERGING, norm=TwoSlopeNorm(0, -lim, lim),
                   aspect="auto")
    for i in range(lfc.shape[0]):
        for j in range(lfc.shape[1]):
            v = lfc.iat[i, j]
            if pd.isna(v):
                ax.add_patch(plt.Rectangle((j - 0.5, i - 0.5), 1, 1, fc="#d8d7d2", ec=SURFACE, lw=1))
                ax.text(j, i, "n/m", ha="center", va="center", fontsize=6, color=INK2)
                continue
            if not linked.iat[i, j]:
                ax.add_patch(plt.Rectangle((j - 0.5, i - 0.5), 1, 1, fill=False, hatch="////",
                                           ec="#9a9994", lw=0))
            if sig.iat[i, j] and linked.iat[i, j]:
                ax.text(j, i, "*", ha="center", va="center", fontsize=11, fontweight="bold",
                        color="#ffffff" if abs(v) > 1.2 else INK)
            # cell gaps
            ax.add_patch(plt.Rectangle((j - 0.5, i - 0.5), 1, 1, fill=False, ec=SURFACE, lw=1.5))
    ax.set_xticks(range(len(DISEASES)))
    ax.set_xticklabels([f"{d}\n{t}" for d, t in zip(DISEASES, ["GSE34526\ngranulosa", "GSE16415\nvisc. adipose",
                                                                  "GSE48452\nliver", "GSE55200\nsubc. adipose"])],
                       fontsize=8)
    ax.xaxis.tick_top()
    ax.set_yticks(range(lfc.shape[0]))
    ax.set_yticklabels(lfc.index, fontsize=8)
    ax.tick_params(length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    cb = fig.colorbar(im, ax=ax, orientation="horizontal", fraction=0.025, pad=0.02, aspect=30)
    cb.set_label("log2 fold change, disease vs control (clipped at ±2)", fontsize=8)
    cb.outline.set_visible(False)
    ax.set_title(f"Expression change of the top {n} hub genes\n"
                 "* significant (p<0.05) in an associated disease;  hatched = gene not associated with "
                 "that disease;  n/m = not measured", loc="left", fontsize=9.5, pad=46)
    fig.tight_layout()
    fig.savefig(FIGURES / "fig3_expression_heatmap.png", bbox_inches="tight")
    plt.close(fig)


def main():
    fig1_graph()
    fig2_table()
    fig3_heatmap()
    print("figures written to", FIGURES)


if __name__ == "__main__":
    main()
