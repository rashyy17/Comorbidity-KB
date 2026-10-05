"""Step 1: resolve disease ontology IDs and pull disease-gene associations from Open Targets."""
import json
import time

import pandas as pd
import requests

from config import DISEASES, EVIDENCE_MIN, EXPECTED_IDS, MIN_ASSOC_SCORE, OPENTARGETS_URL, PROC, RAW

SEARCH_Q = """
query($q: String!) {
  search(queryString: $q, entityNames: ["disease"], page: {index: 0, size: 5}) {
    hits { id name }
  }
}"""

ASSOC_Q = """
query($id: String!, $idx: Int!, $size: Int!) {
  disease(efoId: $id) {
    id name dbXRefs
    associatedTargets(page: {index: $idx, size: $size}) {
      count
      rows { score target { id approvedSymbol approvedName } datatypeScores { id score } }
    }
  }
}"""


def gql(query, variables, retries=4):
    for attempt in range(retries):
        try:
            r = requests.post(OPENTARGETS_URL, json={"query": query, "variables": variables}, timeout=120)
            r.raise_for_status()
            body = r.json()
            if body.get("errors"):
                raise RuntimeError(body["errors"])
            return body["data"]
        except Exception as e:  # network hiccups: back off and retry
            if attempt == retries - 1:
                raise
            print(f"  retry {attempt + 1}: {e}")
            time.sleep(5 * (attempt + 1))


def resolve(short, term):
    hits = gql(SEARCH_Q, {"q": term})["search"]["hits"]
    exact = [h for h in hits if h["name"].lower().replace("disorder", "").strip() in term.lower()
             or term.lower() in h["name"].lower()]
    hit = exact[0] if exact else hits[0]
    if hit["id"] != EXPECTED_IDS[short]:
        print(f"  WARNING: {short} resolved to {hit['id']} ({hit['name']}), expected {EXPECTED_IDS[short]}; "
              f"using expected ID")
        hit = {"id": EXPECTED_IDS[short], "name": hit["name"]}
    return hit


def fetch_associations(disease_id, page_size=1000):
    rows, idx = [], 0
    meta = None
    while True:
        d = gql(ASSOC_Q, {"id": disease_id, "idx": idx, "size": page_size})["disease"]
        meta = meta or d
        page = d["associatedTargets"]["rows"]
        rows.extend(page)
        # Rows come sorted by score descending: stop once we fall below the threshold.
        if not page or page[-1]["score"] < MIN_ASSOC_SCORE or len(rows) >= d["associatedTargets"]["count"]:
            break
        idx += 1
    return meta, rows


def main():
    diseases, assoc = [], []
    for short, term in DISEASES.items():
        hit = resolve(short, term)
        meta, rows = fetch_associations(hit["id"])
        (RAW / "opentargets" / f"{short}.json").write_text(json.dumps(rows))
        efo_xrefs = [x for x in meta["dbXRefs"] if x.upper().startswith("EFO")]
        diseases.append({
            "id": short, "name": short, "full_name": meta["name"], "source_id": meta["id"],
            "efo_xref": ";".join(efo_xrefs), "ot_total_associations": meta["associatedTargets"]["count"],
        })
        kept = [r for r in rows if r["score"] >= MIN_ASSOC_SCORE]
        for r in kept:
            dts = {x["id"]: x["score"] for x in r["datatypeScores"]}
            strong = sorted(k for k, v in dts.items() if v >= EVIDENCE_MIN)
            assoc.append({"disease": short, "ensembl_id": r["target"]["id"],
                          "symbol": r["target"]["approvedSymbol"], "gene_name": r["target"]["approvedName"],
                          "score": r["score"], "evidence_types": ";".join(strong),
                          # Linked only through drugs (ChEMBL "clinical" evidence), e.g. metformin's targets.
                          "drug_only": strong == ["clinical"],
                          "literature_only": strong == ["literature"],
                          **{f"dt_{k}": round(v, 4) for k, v in dts.items()}})
        print(f"{short}: {meta['id']} ({meta['name']}) total={meta['associatedTargets']['count']} "
              f"kept(score>={MIN_ASSOC_SCORE})={len(kept)}")
    pd.DataFrame(diseases).to_csv(PROC / "diseases.csv", index=False)
    pd.DataFrame(assoc).to_csv(PROC / "associations.csv", index=False)


if __name__ == "__main__":
    main()
