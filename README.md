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

## 🔬 Target Species & Genomes Analyzed (CC-Clade)

### 🧬 Rationale for Focusing on the CC Genome Clade
Species possessing the **CC genome** (and its polyploid combinations **BBCC** and **CCDD**) represent a critical lineage in the *Oryza* genus for several key evolutionary reasons:

1. **Global Geographic & Ecological Diversity:** Unlike the cultivated AA genome species, CC-bearing species are widely distributed across diverse global ecosystems—spanning Asia, Africa, and the Americas. This global distribution provides a unique evolutionary framework to study how TE dynamics correlate with geographic adaptation and environmental stress.
2. **Polyploid & Hybridization Dynamics:** Comparing diploid CC species (*O. officinalis*) with allopolyploids (**BBCC**: *O. minuta*, *O. malampuzhaensis*; **CCDD**: *O. alta*, *O. latifolia*, *O. grandiglumis*) allows us to investigate how transposable elements behave during genome merging (polyploidization) and whether TE family expansions/contractions drive structural changes between subgenomes.
3. **Drivers of Structural Variation & Inversions:** TEs are major catalysts for chromosomal inversions, translocations, and synteny breakpoints. Studying TE activity across these 6 representative genomes helps elucidate the precise mechanisms underlying genome-rearrangement dynamics across the *Oryza* genus.

---

### 📊 Included Assemblies

The dataset currently includes **6 target assemblies** spanning CC, BBCC, and CCDD genome types:

| Genome Type | Species | Assembly Accession | Status |
| :--- | :--- | :--- | :--- |
| **CC** | *Oryza officinalis* | `GCA_008326285.1` | 🟢 Retained |
| **BBCC** | *Oryza minuta* | `GCA_048166525.1` | 🟢 Retained |
| **BBCC** | *Oryza malampuzhaensis* | `GCA_048564985.1` | 🟢 Retained |
| **CCDD** | *Oryza alta* | `GCA_047899615.1` | 🟢 Retained |
| **CCDD** | *Oryza latifolia* | `GCA_048174585.1` | 🟢 Retained |
| **CCDD** | *Oryza grandiglumis* | `GCA_048188845.1` | 🟢 Retained |

> 🟢 **Total Genomes:** 6 / 6 assemblies configured for downstream TE dynamics, RiTE DB validation, and structural rearrangement analysis.

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
* **Investigators:** Akmaral Batyrguzhina, Daniil Gerassimov
* **Advisors / PIs:** Michael, Dr. Rod A. Wing
* **Institution:** University of Arizona (UArizona) / W.M. Keck Center for Genome Sciences
