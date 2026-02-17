# Component Analysis Pipeline

This repository contains a comprehensive pipeline for analyzing single-cell RNA-seq data using **Consensus Non-negative Matrix Factorization (cNMF)**. The workflow identifies latent gene programs (components) and characterizes them using downstream functional enrichment (ORA, GSEA) and differential expression (DE) comparisons.

## 📂 Repository Structure

### Core Analysis
*   **`runNMF.ipynb`**: The entry point for the analysis.
    *   **Function**: Performs data preprocessing, runs cNMF factorization across multiple `k` (number of topics), and generates usage/spectra matrices.
    *   **Input**: `data/adata_ge1.h5ad` (AnnData object).
    *   **Key Parameters**: `K` range (e.g., 22, 46), worker threads for parallelization.
    *   **Output**: Factorization results in `cnmf_run/` or `ml_program/`.

### Functional Enrichment
These scripts assign biological meaning to the discovered programs.

*   **`runORA.py`** (Over-Representation Analysis):
    *   **Function**: Tests for enrichment of known pathways (e.g., KEGG, GO) in the **top 300 genes** of each program.
    *   **Method**: Uses `gseapy.enrichr`.
    *   **Output**: `pathway_analysis_results/ORA_summary_all_programs.csv`.
*   **`runGSEA.py`** (Gene Set Enrichment Analysis):
    *   **Function**: Performs enrichment on the **entire ranked list** of genes for each program (based on spectra scores).
    *   **Method**: Uses `gseapy.prerank`.
    *   **Output**: `GO_prerank_results_res2d.csv`, `GO_all_top_programs.csv`.

### Downstream & Comparative Analysis
*   **`de_program_analysis.ipynb`**:
    *   **Function**: Integrates cNMF programs with external Differential Expression (DE) data (e.g., Koenig, Chaffin datasets).
    *   **Method**: Checks if DE gene sets (Up/Down/Sig) are enriched in specific cNMF programs using both ORA and GSEA.
    *   **Input**: `data/Koenig_DE.csv`, `data/chaffin_DE.csv` and cNMF results.
*   **`matchGuides.ipynb`**:
    *   **Function**: Correlates program usage with CRISPR guide assignments to identify phenotype-driving programs.
    *   **Input**: Guided AnnData (`ml_program/adata_ge1_HQstrict_guides.h5ad`) and Usage Matrix (`ml_program/k22_usage.csv`).
*   **`makeUpsetPlot.ipynb`** (formerly `upsetPlot.ipynb`):
    *   **Function**: Visualizes the intersections of gene sets (e.g., shared genes between programs or DE lists) using UpSet plots.

## 🚀 Usage Guide

### 1. Environment Setup
Ensure you have a Python environment with the necessary dependencies:
*   `scanpy`
*   `cnmf` (install via `pip install cnmf` or from source)
*   `gseapy`
*   `pandas`, `numpy`, `matplotlib`, `tqdm`, `upsetplot`

### 2. Running the Pipeline

**Step 1: Factorization**
Open `runNMF.ipynb` and execute the cells to:
1.  Load your `.h5ad` data.
2.  Decompose variation into `K` programs.
3.  Select the optimal `K` based on stability/error plots.

**Step 2: Characterization**
Run the enrichment scripts to label your programs:
```bash
# Run Over-Representation Analysis
python runORA.py

# Run Gene Set Enrichment Analysis
python runGSEA.py
```
*Note: Ensure the script parameters (`name`, `output_dir`, `chosen_k`) match your cNMF run settings.*

**Step 3: Integration**
Use `de_program_analysis.ipynb` to overlap your programs with differential expression lists. This helps answer: *"Does Program X correspond to a known disease signature?"*

**Step 4: Guide Mapping**
Use `matchGuides.ipynb` to associate your labeled programs with specific CRISPR perturbations.

## 📊 Data & Outputs
*   **`data/`**: Raw inputs (AnnData, DE CSVs).
*   **`results/`** & **`cnmf_run/`**: Intermediate factorization files (usage `H`, spectra `W`).
*   **`p22_data/`** & **`pathway_analysis_results/`**: Final CSV tables of enriched pathways and statistics.

