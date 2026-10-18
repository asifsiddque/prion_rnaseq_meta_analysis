# Sample Exclusion Rationale
## Prion Disease RNA-Seq Meta-Analysis Pipeline

**Document version:** 1.0  
**Prepared by:** `04_build_coldata.R`  
**Last updated:** 2026-05-16

---

## 1. Summary

| Metric | N |
|--------|---|
| SRR accessions quantified by Salmon | 867 |
| Samples excluded (pre-analysis) | 125 |
| Samples retained for analysis | **742** |

Exclusions were applied at the coldata-building stage (`04_build_coldata.R`) before
any differential expression analysis. All exclusion decisions are based on
biological incompatibility with the study's primary question — identifying
transcriptomic signatures of prion disease in total RNA from relevant cell
populations — and are independent of any quality metric or statistical outcome.

---

## 2. Per-Study Exclusion Table

| GEO Accession | Stratum | N quantified | N excluded | N retained | Exclusion reason |
|---------------|---------|-------------|-----------|-----------|-----------------|
| GSE149805 | A (CNS Discovery) | 186 | 10 | 176 | noCre negative controls |
| GSE198063 | A (CNS Discovery) | 192 | 101 | 91 | RiboTag IP fraction |
| GSE124932 | D (CNS Replication) | 24 | 8 | 16 | SARM1−/− knockout arms |
| GSE202275 | D (CNS Replication) | 30 | 6 | 24 | VU846-treated arms |
| GSE144738 | A (CNS Discovery) | 60 | 0 | 60 | — |
| GSE189527 | A (CNS Discovery) | 72 | 0 | 72 | — |
| GSE201249 | A (CNS Discovery) | 48 | 0 | 48 | — |
| GSE210128 | B (Peripheral) | 173 | 0 | 173 | — |
| GSE214373 | C (Human sCJD) | 14 | 0 | 14 | — |
| GSE90977 | D (CNS Replication) | 48 | 0 | 48 | — |
| GSE207251 | D (CNS Replication) | 20 | 0 | 20 | — |
| **Total** | | **867** | **125** | **742** | |

---

## 3. Per-Study Exclusion Narratives

### 3.1 GSE149805 — noCre negative controls (N = 10 excluded)

**Study design:** GSE149805 used a Cre-dependent RiboTag system to profile
ribosome-associated mRNA from specific CNS cell types (neurons, astrocytes,
microglia). Mice carrying the RiboTag allele but lacking a Cre driver
(*noCre* animals) were included as negative controls to assess background
immunoprecipitation.

**Exclusion rationale:** The *noCre* samples do not represent any defined
cell-type-specific transcriptome. They capture non-specific background
immunoprecipitation signal and are not biologically comparable to the
Cre-expressing experimental samples. Including them in a prion vs. control
comparison would introduce a spurious cell-type confound and dilute
cell-type-specific prion signatures.

**Decision:** All 10 *noCre* samples excluded. The 176 retained samples
(Cre-expressing, cell-type-stratified) are analysed with a design that
accounts for cell type (`~ time_point + condition`, applied per cell-type
stratum in Script 05).

---

### 3.2 GSE198063 — RiboTag IP fraction (N = 101 excluded)

**Study design:** GSE198063 profiled both the total RNA input fraction and
the RiboTag immunoprecipitation (IP) fraction from the same animals. The IP
fraction captures ribosome-associated (actively translated) mRNA; the input
fraction captures total cellular RNA.

**Exclusion rationale:** This meta-analysis integrates studies that measure
total cellular RNA (or bulk RNA-seq without IP enrichment). The RiboTag IP
fraction represents a fundamentally different RNA population — actively
translated transcripts — and is not directly comparable to total RNA
measurements from other studies. Mixing IP and total RNA fractions in a
meta-analysis would introduce systematic between-study heterogeneity that
cannot be corrected by the random-effects model.

**Decision:** All 101 IP-fraction samples excluded. The 91 retained input-fraction
samples are analysed with `~ condition` (Script 05). This is the largest single
exclusion in the dataset and is the primary reason the retained count (742) is
substantially lower than the quantified count (867).

**Note:** The IP fraction data from GSE198063 is scientifically valuable for
studying translational regulation in prion disease and could be the subject of
a separate analysis. It is excluded here solely because it is incompatible with
the total-RNA meta-analysis framework.

---

### 3.3 GSE124932 — SARM1−/− knockout arms (N = 8 excluded)

**Study design:** GSE124932 compared prion-infected wild-type (WT) mice to
prion-infected SARM1−/− knockout mice, alongside uninfected controls of both
genotypes. SARM1 (Sterile Alpha and TIR Motif Containing 1) is a key regulator
of axon degeneration (Wallerian degeneration) and innate immune signalling.

**Exclusion rationale:** SARM1−/− mice have constitutively altered axon
degeneration pathways, neuroinflammatory responses, and baseline transcriptomes
that differ from WT mice independent of prion infection. Including SARM1−/−
samples in a meta-analysis of prion disease transcriptomics would confound the
prion signal with the SARM1 knockout effect. The scientific question of this
pipeline is the conserved prion transcriptomic signature in a normal genetic
background.

**Decision:** All 8 SARM1−/− samples (4 prion-infected, 4 uninfected) excluded.
The 16 retained WT samples are analysed with `~ condition` (Script 05, Stratum D
replication).

---

### 3.4 GSE202275 — VU846-treated arms (N = 6 excluded)

**Study design:** GSE202275 examined the effect of VU846, a small-molecule
activator of the integrated stress response (ISR), in prion-infected mice.
The study included vehicle-treated and VU846-treated arms, both in prion-infected
and uninfected animals.

**Exclusion rationale:** VU846 is a pharmacological intervention that directly
modulates the ISR — a pathway known to be activated in prion disease. Including
drug-treated samples in a meta-analysis of the natural prion transcriptomic
response would confound the disease signature with the drug effect. The pipeline
aims to characterise the unperturbed prion transcriptome.

**Decision:** All 6 VU846-treated samples excluded. The 24 retained vehicle-treated
samples are analysed with `~ condition` (Script 05, Stratum D replication).

---

## 4. Ready-to-Use Methods Paragraph

The following paragraph is suitable for direct inclusion in the Methods section
of a manuscript. Adjust journal-specific formatting (italics, citation style) as needed.

---

> **Sample selection and exclusion criteria.** Of the 867 RNA-seq samples
> quantified by Salmon across 11 GEO datasets, 125 samples were excluded prior
> to differential expression analysis on the basis of biological incompatibility
> with the study's primary question. Exclusions were determined before any
> statistical analysis and are independent of expression data or quality metrics.
>
> From GSE149805 (N = 186 quantified), 10 *noCre* negative-control samples were
> excluded because they lack a Cre driver and therefore do not represent a
> defined cell-type-specific transcriptome; they capture non-specific background
> immunoprecipitation signal and are not comparable to the Cre-expressing
> experimental samples.
>
> From GSE198063 (N = 192 quantified), 101 RiboTag immunoprecipitation (IP)
> fraction samples were excluded. The IP fraction captures ribosome-associated
> (actively translated) mRNA, which represents a fundamentally different RNA
> population from the total cellular RNA measured in all other studies in this
> meta-analysis. Including IP-fraction samples would introduce systematic
> between-study heterogeneity incompatible with the random-effects meta-analysis
> framework. The 91 retained input-fraction samples from GSE198063 represent
> total cellular RNA and are directly comparable to the other datasets.
>
> From GSE124932 (N = 24 quantified), 8 SARM1−/− knockout samples were excluded.
> SARM1 is a key regulator of Wallerian degeneration and innate immune signalling;
> SARM1−/− mice have constitutively altered axon degeneration pathways and
> baseline transcriptomes that differ from wild-type animals independent of prion
> infection. Only wild-type samples (N = 16) were retained for replication analysis.
>
> From GSE202275 (N = 30 quantified), 6 VU846-treated samples were excluded.
> VU846 is a pharmacological activator of the integrated stress response (ISR),
> a pathway known to be dysregulated in prion disease; drug-treated samples
> would confound the natural prion transcriptomic response. Only vehicle-treated
> samples (N = 24) were retained.
>
> After exclusions, 742 samples across 11 datasets were retained for analysis
> (Stratum A: 447 samples across 5 CNS discovery studies; Stratum B: 173 samples,
> GSE210128 peripheral longitudinal; Stratum C: 14 samples, GSE214373 human sCJD;
> Stratum D: 108 samples across 4 CNS replication studies).

---

## 5. Implementation

Exclusions are implemented in `scripts/04_build_coldata.R` via the following
per-study filters applied during coldata construction:

| Study | Filter applied |
|-------|---------------|
| GSE149805 | `coldata <- coldata[coldata$cre_status != "noCre", ]` |
| GSE198063 | `coldata <- coldata[coldata$fraction == "input", ]` |
| GSE124932 | `coldata <- coldata[coldata$genotype == "WT", ]` |
| GSE202275 | `coldata <- coldata[coldata$treatment == "vehicle", ]` |

The resulting per-study coldata CSVs (written to `config/coldata/`) contain only
the retained samples. Scripts 05–08 read these CSVs and never see the excluded
samples.

---

## 6. Cross-Reference

- **Script implementing exclusions:** `scripts/04_build_coldata.R` (header, Section 2)
- **Sample metadata (all 867 SRRs):** `config/sample_metadata.tsv`
- **Study registry:** `config/study_registry.csv`
- **Pipeline overview:** `README.md`
