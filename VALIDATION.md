# Validation

Source files: `results/validation.json`, `results/validation_genes.csv`, `results/sensitivity.csv`,
`results/dataset_log2fc_correlation.csv`. Run date 2026-10-05.

## Bottom line

| Check | Result |
|---|---|
| INSR is a hub | **Yes**: all 4 diseases at score ≥ 0.1 and ≥ 0.2; 2 diseases at ≥ 0.3 |
| IRS1 is a hub | **Yes**: 4 diseases at ≥ 0.1; 3 at ≥ 0.2 and ≥ 0.3 |
| IGF1 is a hub | **Only at the primary threshold (≥ 0.1)**. All four IGF1 edges score 0.11–0.13 and rest on literature (± animal model) evidence only. At ≥ 0.2 IGF1 has degree 0. **Fragile.** |
| INSR expression-consistent | **Yes, but not in PCOS/T2D**: significantly *down* in NAFLD liver (log2FC −0.17, p = 0.024) and obese adipose (−0.58, p = 1e-8). In PCOS (−1.25, p = 0.07) and T2D (−0.46, p = 0.53) it trends down but is not significant. |
| IRS1 expression-consistent | **No**: significant only in Obesity (down). Down but not significant in T2D (p = 0.08) and PCOS; flat in NAFLD. |
| IGF1 expression-consistent | **No**: significant only in NAFLD (down). Down but not significant in the other three. |
| Hubs enriched for KEGG Insulin resistance genes | **Yes**: 28/535 KEGG-mapped hubs vs 26/1,224 non-hubs; OR 2.54, Fisher p = 6.8e-4 |
| Expression consistency beats a random-gene null | **Yes**: 41 consistent vs 12.1 expected (95% null range 6–20), p < 0.001 |
| Expression consistency beats the genome-wide background | **No**: 82% concordance among hubs vs 75.6% among *all* genes significant in 2+ datasets; binomial p = 0.19 |

**Honest summary.** The pipeline recovers the three literature genes as hubs, and the hub set is
enriched for insulin-resistance biology, so the association layer behaves sensibly. The **expression**
half of the writeup's validation is not met. No validation gene is significantly dysregulated in
*both* the PCOS and T2D datasets, the pair the writeup emphasizes. Those two datasets are small
(3 vs 7 and 5 vs 5). The validation genes all trend *down* in PCOS and T2D, but none of these trends
is significant. IGF1's hub status depends on the choice of threshold.

## 1. Validation genes in detail

| Gene | Disease | OT score | Evidence types (score ≥ 0.1) | GEO log2FC | p | Significant |
|---|---|---|---|---|---|---|
| INSR | NAFLD | 0.354 | clinical, genetic_association, literature | −0.174 | 0.024 | yes |
| INSR | Obesity | 0.293 | animal_model, clinical, literature | −0.580 | 1.3e-8 | yes |
| INSR | PCOS | 0.265 | genetic_association, literature | −1.249 | 0.071 | no |
| INSR | T2D | 0.770 | animal_model, clinical, genetic_association, genetic_literature, literature, somatic_mutation | −0.462 | 0.527 | no |
| IRS1 | NAFLD | 0.360 | genetic_association, literature | +0.047 | 0.750 | no |
| IRS1 | Obesity | 0.500 | animal_model, genetic_association, literature | −0.305 | 0.030 | yes |
| IRS1 | PCOS | **0.105** | **literature only** | −0.343 | 0.657 | no |
| IRS1 | T2D | 0.659 | animal_model, genetic_association, genetic_literature, literature | −1.004 | 0.077 | no |
| IGF1 | NAFLD | **0.114** | **literature only** | −0.704 | 0.004 | yes |
| IGF1 | Obesity | 0.129 | animal_model, literature | −0.312 | 0.075 | no |
| IGF1 | PCOS | **0.110** | **literature only** | −0.358 | 0.451 | no |
| IGF1 | T2D | 0.126 | animal_model, literature | −0.236 | 0.538 | no |

Observations:
- **Direction agrees with the biology.** Reduced insulin-receptor and IRS1 expression in insulin-resistant
  tissue is well described. 11 of 12 fold changes are negative; the exception is IRS1 in NAFLD (+0.05, flat).
- **The PCOS–INSR/IRS1/IGF1 links in Open Targets are weak.** IRS1 and IGF1 reach PCOS only through text
  mining, just above the threshold. INSR–PCOS has genetic support (score 0.265).
- **The T2D dataset (GSE16415) cannot confirm insulin-signalling genes.** INSR p = 0.53 and IRS1 p = 0.08
  with n = 5 vs 5. The same dataset does show the expected GLUT4 (SLC2A4, p = 0.002), ADIPOQ (p = 0.025)
  and PPARG (p = 0.02) decreases, so it is underpowered rather than broken.

## 2. Threshold sensitivity (`results/sensitivity.csv`)

| min score | exclude drug-only | hub genes | in all 4 | in 3+ | expr-consistent | discordant | hub pathways (enriched ≥ 2) | INSR deg. | IRS1 deg. | IGF1 deg.* |
|---|---|---|---|---|---|---|---|---|---|---|
| 0.1 | no | 855 | 99 | 255 | 41 | 9 | 143 | 4 ✓consistent | 4 | 4 |
| 0.2 | no | 483 | 61 | 151 | 21 | 4 | 57 | 4 ✓consistent | 3 | <2 |
| 0.3 | no | 289 | 55 | 76 | 16 | 2 | 36 | 2 | 3 | <2 |
| 0.1 | **yes** | 753 | 48 | 166 | 36 | 8 | 148 | 4 ✓consistent | 4 | 4 |

\* "<2" means not a hub; the hub query does not report the degree of non-hubs.

Excluding drug-only edges halves the 4-disease hub set (99 → 48). That confirms that about 50 of the
"most shared" genes are drug targets, mostly metformin's mitochondrial complex I subunits, propagated to
every disease the drug is trialled in. The validation genes are unaffected.

## 3. Is the hub set biologically meaningful? (unbiased reference set)

Reference set: members of KEGG **hsa04931 Insulin resistance**, restricted to genes in the KB.

|  | in hsa04931 | not in hsa04931 |
|---|---|---|
| hub gene | 28 | 507 |
| non-hub gene | 26 | 1,198 |

Odds ratio 2.54, one-sided Fisher p = 6.8 × 10⁻⁴. Hubs in the set include INSR, IRS1, INS, AKT1, FOXO1,
GSK3B, PTEN, MTOR, PIK3R1, SLC2A2, SLC2A4, PPARA, PPARGC1A, SREBF1, TNF, IL6, NFKB1, STAT3, CD36,
PRKAA1/2 and others. Insulin resistance (hsa04931) is also itself a hub pathway, enriched in all 4 diseases
(FDR 2.5e-11), together with AMPK signalling, adipocytokine signalling, FoxO signalling,
Type II diabetes mellitus and NAFLD (hsa04932).

## 4. Is "expression-consistent" a real signal?

Two nulls answer two different questions:

1. **Random-gene null (gene labels permuted within each dataset, 2,000 permutations).** Hubs produce
   41 consistent and 9 discordant genes. Random genes produce 12.1 consistent on average (95% range 6–20)
   and 48.8% concordance. Both are p < 0.001. Hub genes are more often dysregulated in two diseases at
   once, and more often in the same direction, than random genes.
2. **Genome-wide background (no permutation).** The four datasets are positively correlated genome-wide
   (Spearman ρ of log2FC: Obesity–T2D 0.41, PCOS–T2D 0.27, NAFLD–PCOS 0.20, NAFLD–Obesity 0.16,
   Obesity–PCOS 0.16, NAFLD–T2D 0.09). Among *all* 2,119 genes significant in 2+ datasets, **75.6%** agree
   in direction. Hubs reach 82.0%, which is **not significantly higher** (binomial p = 0.19).

Interpretation: "expression-consistent" is real in the sense that the direction agreement is not a coin flip.
But most of it reflects a shared global disease signature, probably inflammation and adipose remodelling
plus technical factors. It is **not specific to hub genes**. Treat each consistent hub as a lead to follow
up, not as proof of a shared mechanism.

## 5. What would strengthen validation

- Larger PCOS and T2D cohorts on the same tissue (e.g. adipose for all four) to remove tissue confounding.
- A pre-registered positive-control list (e.g. DisGeNET curated PCOS∩T2D genes) and a negative-control
  disease (e.g. an unrelated cancer) run through the same pipeline.
- Fold-change floors and FDR-based significance once sample sizes allow.
