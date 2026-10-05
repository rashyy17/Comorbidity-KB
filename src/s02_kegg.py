"""Step 2: pull the KEGG human gene->pathway link table and map Open Targets genes onto it."""
import time

import pandas as pd
import requests

from config import KEGG_URL, PROC, RAW


def fetch(endpoint):
    path = RAW / "kegg" / (endpoint.replace("/", "_") + ".tsv")
    if not path.exists():
        for attempt in range(4):
            try:
                r = requests.get(f"{KEGG_URL}/{endpoint}", timeout=120)
                r.raise_for_status()
                path.write_text(r.text)
                break
            except Exception as e:
                if attempt == 3:
                    raise
                print(f"  retry {attempt + 1}: {e}")
                time.sleep(5 * (attempt + 1))
    return path


def main():
    links = pd.read_csv(fetch("link/pathway/hsa"), sep="\t", header=None, names=["kegg_gene", "kegg_path"])
    links["entrez_id"] = links.kegg_gene.str.replace("hsa:", "", regex=False)
    links["kegg_id"] = links.kegg_path.str.replace("path:", "", regex=False)

    pathways = pd.read_csv(fetch("list/pathway/hsa"), sep="\t", header=None, names=["kegg_id", "name"])
    pathways["name"] = pathways.name.str.replace(r" - Homo sapiens \(human\)$", "", regex=True)
    # KEGG pathway classes: 011xx/012xx are global/overview maps, 05xxx are human diseases.
    num = pathways.kegg_id.str[3:].astype(int)
    pathways["category"] = "pathway"
    pathways.loc[(num >= 1100) & (num < 1300), "category"] = "overview"
    pathways.loc[num >= 5000, "category"] = "disease"

    genes = pd.read_csv(fetch("list/hsa"), sep="\t", header=None, names=["kegg_gene", "type", "loc", "desc"])
    genes["entrez_id"] = genes.kegg_gene.str.replace("hsa:", "", regex=False)
    names = genes.desc.str.split(";").str[0].str.split(",")
    genes["primary"] = names.str[0].str.strip()
    genes["aliases"] = names.apply(lambda xs: [x.strip() for x in xs])

    assoc = pd.read_csv(PROC / "associations.csv")
    symbols = set(assoc.symbol)

    # Primary-symbol match first, then alias match for anything left (only if unambiguous).
    sym2entrez = dict(zip(genes.primary, genes.entrez_id))
    alias_map = {}
    for eid, aliases in zip(genes.entrez_id, genes.aliases):
        for a in aliases[1:]:
            alias_map.setdefault(a, set()).add(eid)
    mapping, how = {}, {}
    for s in symbols:
        if s in sym2entrez:
            mapping[s], how[s] = sym2entrez[s], "primary"
        elif s in alias_map and len(alias_map[s]) == 1:
            mapping[s], how[s] = next(iter(alias_map[s])), "alias"
    gmap = pd.DataFrame({"symbol": list(mapping), "entrez_id": list(mapping.values()),
                         "match": [how[s] for s in mapping]})
    gmap.to_csv(PROC / "gene_entrez.csv", index=False)

    gp = links.merge(gmap, on="entrez_id")[["symbol", "entrez_id", "kegg_id"]]
    gp.to_csv(PROC / "gene_pathway.csv", index=False)
    # Universe for enrichment: every human gene with >=1 KEGG pathway.
    pd.Series(links.entrez_id.unique(), name="entrez_id").to_csv(PROC / "kegg_universe.csv", index=False)
    pathways.to_csv(PROC / "pathways.csv", index=False)

    print(f"KEGG links={len(links)} pathways={len(pathways)}; OT symbols={len(symbols)} "
          f"mapped to Entrez={len(gmap)} ({(gmap.match == 'alias').sum()} via alias); "
          f"with >=1 pathway={gp.symbol.nunique()}; gene-pathway edges={len(gp)}")


if __name__ == "__main__":
    main()
