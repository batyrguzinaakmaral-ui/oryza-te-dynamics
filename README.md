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

## 🔬 Target Species (CC-Genome Focus)

| Genome Type | Species | Geographic Distribution |
| :--- | :--- | :--- |
| **CC** | *Oryza officinalis*, *Oryza rhizomatis*, *Oryza eichingeri* | Asia / Africa |
| **BBCC** | *Oryza punctata* (allopolyploid), *Oryza minuta*, *Oryza malampuzhaensis* | Africa / Asia |
| **CCDD** | *Oryza latifolia*, *Oryza alta*, *Oryza grandiglumis* | Americas |

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
