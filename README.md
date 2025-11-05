# 🧬 Asgard Protein Annotation Pipeline

A modular workflow for annotation and structural prediction of **Asgard archaeal proteins**, integrating sequence clustering, annotation, and structure inference with **ESMFold**.

---

## 📘 Overview

This repository automates the process of:
- Filtering and quality-checking Asgard archaeal genomes
- Clustering protein sequences using **MMseqs2**
- Annotating proteins via **EggNOG-Mapper** and **Prokka**
- Predicting 3D structures using **ESMFold**
- Summarizing structural confidence (pLDDT) and cluster information

The pipeline supports both **local Conda environments** and **Docker-based ESMFold inference**.

---

## 🧰 Repository Structure

| File / Script | Description |
|----------------|-------------|
| `snakefile` | Main Snakemake pipeline controlling all steps |
| `Genome_QA.sh` | Checks genome completeness and contamination |
| `filter_genomes.py` | Filters genomes based on QC thresholds |
| `get_filtered_genomes_names.py` | Extracts filtered genome names |
| `gtdb_filter.py` | Filters genomes using GTDB taxonomy |
| `merge_faa.sh` | Merges all `.faa` files into one dataset |
| `asgard_cluster_pipeline.sh` | Main script for MMseqs2-based clustering |
| `clustering_diag.sh` | Diagnostics and stats for cluster distribution |
| `cluster_selection.sh` | Selects clusters with ≥5 members |
| `mmseqs_clustering_summary.py` | Summarizes MMseqs2 clustering output |
| `run_eggnog_mapper.sh` | Functional annotation using EggNOG-Mapper |
| `run_prokka_all.sh` | Runs Prokka annotation for all genomes |
| `esm_inference.py` | ESMFold structure prediction script (Python) |
| `esmfold_runner.sh` | Docker wrapper for running ESMFold on FASTA inputs |
| `missing_faas.txt` | List of missing protein files for debugging |
| `compile_cds_info.py` | Compiles all the proteins of the given completeness level; filter_genomes.py |

---

## ⚙️ Installation

### 1️⃣ Create environment
```bash
conda create -n asgard_anno python=3.9 -y
conda activate asgard_anno
