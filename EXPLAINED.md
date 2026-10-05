# The Comorbidity-Cluster Knowledge Base, explained

*A plain-language guide for readers without a background in genomics or databases.*

## The question

Polycystic ovary syndrome (**PCOS**), type 2 diabetes (**T2D**), non-alcoholic fatty liver disease
(**NAFLD**) and **obesity** often turn up in the same patients. A woman with PCOS is several times more
likely to develop T2D, and most people with NAFLD are overweight. Doctors call this clustering
*comorbidity*.

One explanation is that these conditions share some of the same underlying machinery: the same genes and
biological processes go wrong in each. The project tests that idea by putting what is known about all
four diseases into one connected database and asking: **which genes and processes are shared, and do
they misbehave the same way in each disease?**

## Three kinds of evidence

Think of a gene as a single instruction in the body's manual, and a **pathway** as a chapter: a group
of genes that work together on one job, such as responding to insulin.

1. **Which genes are linked to each disease?** From **Open Targets**, a public database that combines
   many kinds of evidence: DNA studies of large patient groups, animal experiments, drug trials and the
   scientific literature. Each disease-gene link gets a confidence score from 0 to 1. Links scoring at
   least 0.1 are kept.
2. **Which pathway does each gene belong to?** From **KEGG**, a long-standing reference that sorts human
   genes into about 370 pathways.
3. **Is the gene turned up or down in patients?** From **GEO**, a public archive of experiments that
   measured the activity of about 20,000 genes at once in tissue from patients and from healthy people.
   One experiment per disease is used:

   | Disease | Tissue measured | Patients vs healthy |
   |---|---|---|
   | PCOS | cells surrounding the egg in the ovary | 7 vs 3 |
   | T2D | belly fat (around the organs) | 5 vs 5 |
   | NAFLD | liver biopsy | 32 vs 14 |
   | Obesity | fat under the skin | 16 vs 7 |

   For every gene the pipeline computes the **fold change**: how much more (up) or less (down) active
   it is in patients. It also computes a **p-value**: how likely a difference that large would be by
   chance alone. Below 0.05 counts as "significant" here.

## How the pieces fit together (the schema)

The data are stored as a **graph**: dots (nodes) connected by labelled arrows (relationships). This
shape suits questions like "what do these diseases have in common?", which are really questions about
paths between dots.

```
           ASSOCIATED_WITH (score)              PARTICIPATES_IN
 [Disease] ───────────────────────► [Gene] ─────────────────────► [Pathway]
     ▲                                 │
     └──────── EXPRESSED_IN ───────────┘
          (fold change, up/down, p-value)
```

- **Disease** (4): PCOS, T2D, NAFLD, Obesity.
- **Gene** (3,303): every gene linked to at least one of the four.
- **Pathway** (366): KEGG pathways containing at least one of those genes.
- **ASSOCIATED_WITH**: "Open Targets links this gene to this disease, with this confidence".
- **PARTICIPATES_IN**: "this gene is part of this pathway".
- **EXPRESSED_IN**: "in this disease's patient tissue, this gene went up or down by this much".

"Shared" is never stored. The database works it out each time you ask, by following the arrows. The
database software is **Neo4j**. An identical copy in plain SQLite tables means it also works without
Neo4j, and every question was asked of both copies to check the answers match. They do.

## What counts as a "hub"

- A **hub gene** is linked to **two or more** of the four diseases.
- A **hub pathway** is one where a disease's genes are **over-represented**: more of them land in the
  pathway than chance predicts, in two or more diseases. Counting any overlap is not enough. T2D alone
  has over 2,000 genes, so it touches almost every pathway by accident (351 of 366 would count).
- An **expression-consistent hub** is a hub gene that is significantly changed in the patient data for
  at least two of its diseases, **in the same direction** (up in both, or down in both). That is
  stronger evidence than "appears on two lists".

## What was found

**Shared genes.** 855 genes are linked to two or more diseases, and 99 to all four. T2D and obesity share
the most (677 genes). NAFLD and PCOS share the fewest (112), partly because PCOS has far fewer
known genes overall.

**A trap uncovered along the way.** At first, the "most shared" genes were dominated by a group of
genes for the cell's energy factories (mitochondrial *complex I*). Checking the evidence showed why.
Metformin, a common diabetes drug that acts on complex I, is used or tested in all four diseases. Open
Targets therefore links its target genes to every one of them, but that reflects prescribing patterns,
not disease biology. The database now flags such "drug-only" links, and results are shown with and
without them. About half of the all-four-disease genes disappear when they are removed.

**Shared processes.** 143 pathways are over-represented in two or more diseases. Those enriched in all
four include *insulin resistance*, *AMPK signalling* (the cell's fuel gauge), *adipocytokine signalling*
(hormones made by fat tissue), *FoxO signalling*, *type II diabetes mellitus* and the *NAFLD* pathway
itself. This fits the long-standing view that insulin resistance and fat-tissue dysfunction are the
common thread.

**Genes that move together.** 41 hub genes change in the same direction in two or more diseases. Examples:
- **Up:** CCL2 and CD4 (immune and inflammatory signals), MMP9 and TGFB1 (tissue remodelling and scarring),
  and GDF15 (a stress hormone).
- **Down:** PPARG (master regulator of healthy fat cells), ADIPOQ (adiponectin, a protective hormone),
  INSR (the insulin receptor itself), RBP4 and TXNIP.

The overall picture: inflammation goes up, while healthy fat-cell function and insulin sensing go down.

**Did it find the expected genes?** The writeup named three genes known to link PCOS and diabetes:
INSR, IRS1 and IGF1.
- **INSR:** found as a hub of all four diseases, and down in liver and fat.
- **IRS1:** found as a hub, but only changed significantly in one dataset.
- **IGF1:** only qualifies at the lowest confidence cutoff, on literature evidence alone.

None of the three changed significantly in *both* the PCOS and diabetes datasets. Those two experiments
are very small (3–5 people per group), which makes real changes hard to detect. All three genes trended
down in both, matching expectations, but the trends were not statistically convincing. Details are in
[VALIDATION.md](VALIDATION.md).

## Limitations

1. **Small patient groups.** The PCOS (3 vs 7) and T2D (5 vs 5) experiments are tiny, so many real
   changes are missed and some flagged changes will be false alarms. A lenient significance cutoff
   (p < 0.05, without correcting for testing 20,000 genes) was used so these datasets contribute at all.
2. **Different tissues.** Ovary cells, two kinds of fat and liver are compared directly. A gene can
   legitimately go up in liver and down in fat, so "inconsistent" does not always mean "unrelated".
3. **The "same direction" signal is not unique to hub genes.** The four experiments resemble each other
   across the whole genome: 76% of *all* genes significant in two datasets move the same way. Hub
   genes reach 82%, which is not a significant difference. Much of the consistency reflects a general
   "sick tissue" signature, such as inflammation, rather than something special about shared genes.
4. **Uneven research attention.** T2D and obesity have been studied far more than PCOS, so they have
   many more known genes and naturally share more. A gene missing from PCOS may simply be unstudied.
5. **Association is not causation.** Open Targets links combine many evidence types, including drug use
   (see the metformin trap) and text mining. A link means "connected in the evidence", not "causes the
   disease".
6. **Thresholds are choices.** Raising the confidence cutoff from 0.1 to 0.2 cuts hub genes from 855 to
   483 and removes IGF1. All key numbers are reported at 0.1, 0.2 and 0.3 so readers can judge.
7. **The NAFLD comparison is partly an obesity comparison.** The liver patients are heavier than the
   controls (median BMI about 48 vs 26). A second comparison against obese people without liver disease
   is included for reference.
8. **Snapshot in time.** Open Targets and KEGG update regularly, so re-running later may change counts
   slightly. Open Targets has also switched most disease IDs from the "EFO" naming scheme in the writeup
   to "MONDO".

## Where to look next

- `results/hub_genes.csv`: every hub gene with diseases, scores, fold changes and consistency class.
- `results/hub_pathways.csv`: shared pathways with enrichment statistics.
- `results/figures/`: the overlap graph, hub table and expression heatmap.
- `queries/hubs.cypher`: ready-made questions to paste into the Neo4j Browser (http://localhost:7474).
