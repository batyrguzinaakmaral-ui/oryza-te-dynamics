# Transposable Element Dynamics and Genome Rearrangements in *Oryza* Genomes

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.10+](https://img.shields.io/badge/Python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![R 4.2+](https://img.shields.io/badge/R-4.2+-blue.svg)](https://www.r-project.org/)
[![Platform: CyVerse / OpenOnDemand](https://img.shields.io/badge/Platform-CyVerse%20%7C%20UA%20HPC-red.svg)](https://cyverse.org/)

## 📌 Project Overview
This repository contains the bioinformatic pipeline, data analysis scripts, and visualization tools for investigating Transposable Element (TE) dynamics across the genus *Oryza*, with a specific focus on species harbouring **CC, BBCC, and CCDD genomes**.

The project aims to assess TE abundance, superfamily composition, and evolutionary dynamics (expansions/contractions), test and validate the **RiTE database**, and correlate TE distribution with large-scale structural variations (inversions and genome rearrangements).

---

## 🎯 Deliverables & Objectives
1. **Standardized TE Annotation:** Genome-wide annotation of representative *Oryza* species using established pipelines (EDTA, RepeatMasker).
2. **RiTE Database Validation:** Evaluation of the RiTE database against newly annotated genomes to identify missing, redundant, or misclassified TE families.
3. **Comparative Analysis (CC-genome Focus):** Comparative profiling of TE content across 9 global species containing CC genome types (3 CC, 3 BBCC, 3 CCDD).
4. **TE Dynamics & Genome Rearrangements:** Correlating TE locations and lineage-specific activity with structural inversions and synteny breakpoints.
5. **Visualization & Reproducible Workflow:** Publication-quality figures, summary tables, and documented computational workflows.
6. **Final Deliverable & Poster:** Written research report and poster presentation summarizing methods, results, and biological insights.

---

## 🔬 Target Species & Genome Assemblies (CC-Clade Focus)

### 🧬 Rationale for Focusing on the CC Genome Clade
Species possessing the **CC genome** (and its polyploid combinations **BBCC** and **CCDD**) represent a critical evolutionary lineage in the genus *Oryza* for several key reasons:

1. **Global Geographic & Ecological Diversity:** Unlike the cultivated AA genome species, CC-bearing species are widely distributed across diverse global ecosystems—spanning Asia, Africa, and the Americas. This global distribution provides a unique framework to study how TE dynamics correlate with biogeography, ecology, and environmental adaptation.
2. **Polyploid & Hybridization Dynamics:** Comparing diploid CC species (*O. officinalis*) with allopolyploids (**BBCC**: *O. minuta*, *O. malampuzhaensis*; **CCDD**: *O. alta*, *O. latifolia*, *O. grandiglumis*) allows us to investigate how transposable elements behave during genome merging (polyploidization) and whether TE family expansions/contractions drive structural divergence between subgenomes.
3. **Drivers of Structural Variation & Inversions:** TEs are major catalysts for chromosomal inversions, translocations, and synteny breakpoints. Studying TE activity across these representative genomes helps elucidate the precise mechanisms underlying genome-rearrangement dynamics across the *Oryza* genus.

---

### 🌐 Conceptual Overview of Target Genotypes

| Genome Type | Representative Species | Geographic Distribution |
| :--- | :--- | :--- |
| **CC** | *Oryza officinalis*, *Oryza rhizomatis*, *Oryza eichingeri* | Asia / Africa |
| **BBCC** | *Oryza punctata* (allopolyploid), *Oryza minuta*, *Oryza malampuzhaensis* | Africa / Asia |
| **CCDD** | *Oryza latifolia*, *Oryza alta*, *Oryza grandiglumis* | Americas |

---

### 📊 Assembly Selection & Retained Status

The initial target set comprises **9 species** across CC, BBCC, and CCDD genome types. Currently, **6 species have retained high-quality assemblies** for downstream TE analysis, while 3 species remain pending/not retained.

| Genome Type | Species | Geographic Distribution | Assembly Accession | Status |
| :--- | :--- | :--- | :--- | :--- |
| **CC** | *Oryza officinalis* | Asia / Africa | `GCA_008326285.1` | 🟢 Retained |
| **CC** | *Oryza eichingeri* | Africa | — | 🔴 Not Retained (Pending assembly) |
| **CC** | *Oryza rhizomatis* | Asia | — | 🔴 Not Retained (Pending assembly) |
| **BBCC** | *Oryza minuta* | Asia | `GCA_048166525.1` | 🟢 Retained |
| **BBCC** | *Oryza malampuzhaensis* | Asia | `GCA_048564985.1` | 🟢 Retained |
| **BBCC** | *Oryza punctata* | Africa | — | 🔴 Not Retained (Pending assembly) |
| **CCDD** | *Oryza alta* | Americas | `GCA_047899615.1` | 🟢 Retained |
| **CCDD** | *Oryza latifolia* | Americas | `GCA_048174585.1` | 🟢 Retained |
| **CCDD** | *Oryza grandiglumis* | Americas | `GCA_048188845.1` | 🟢 Retained |

> 🟢 **Retained Assemblies:** 6 / 9 genomes configured for TE dynamics, RiTE DB validation, and structural rearrangement profiling.  
> 🔴 **Not Retained Assemblies:** 3 / 9 genomes awaiting high-resolution assembly integration.
---

## 📂 Repository Structure

| Directory / File | Description |
| :--- | :--- |
| `data/` | Genome assembly metadata, TE annotations (GFF3/BED), and RiTE db benchmarks |
| `notebooks/` | Interactive Jupyter and RStudio notebooks for comparative analysis and plots |
| `workflow/` | Reproducible execution pipelines (Snakemake / Shell scripts) |
| `src/` | Helper Python modules and R plotting functions |
| `docs/poster/` | Project poster layout, figures, and presentation materials |

---

## 🛠 Tech Stack & Tools
* **Environments:** CyVerse Discovery Environment / UArizona HPC (Open OnDemand)
* **Languages:** Python (Pandas, NumPy, Biopython, Seaborn), R (ggplot2, GenomicRanges, tidyverse)
* **TE Tools:** EDTA, RepeatMasker, BLAST+, RiTE Database

---

## 👥 Project Team & Acknowledgments

* **Investigators:**
  * **Akmaral Batyrguzhina** — [akmaralb@arizona.edu](mailto:akmaralb@arizona.edu)
  * **Daniil Gerassimov** — [daniilgerassimov@arizona.edu](mailto:daniilgerassimov@arizona.edu)
* **PIs:**
  * **Andrea Zuccolo** — [azuccolo@arizona.edu](mailto:azuccolo@arizona.edu)
  * **Md Nafis Ul Alam (Michael)** — [mdalam@arizona.edu](mailto:mdalam@arizona.edu)
  * **Rod A. Wing** — [rwing@arizona.edu](mailto:rwing@arizona.edu)
* **Institution & Facility:**
  * **University of Arizona (UArizona)** — BIO5 Institute (Keating Bioresearch Building) / W.M. Keck Center for Genome Sciences
