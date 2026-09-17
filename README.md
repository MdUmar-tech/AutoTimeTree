# ⏳ AutoTimeTree: Automated Molecular Clock & Time Tree Analysis Pipeline

[![Language: Python](https://img.shields.io/badge/Language-Python_3.8+-3776AB?style=flat&logo=python&logoColor=white)](https://www.python.org/)
[![Language: R](https://img.shields.io/badge/Language-R_4.0+-276DC3?style=flat&logo=r&logoColor=white)](https://www.r-project.org/)
[![Jupyter Notebook](https://img.shields.io/badge/Jupyter-Notebook-orange?style=flat&logo=jupyter&logoColor=white)](AutoTimeTree_Pipeline.ipynb)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

**AutoTimeTree** is an automated, end-to-end bioinformatics software suite for sequence alignment, Maximum Likelihood (ML) phylogenetics, relaxed molecular clock rate-smoothing (timetree/chronogram computation), pairwise genetic distance matrix calculation, publication-ready multi-panel visualization, and Word (`.docx`) report generation.

---

## 📌 Features

- **⚡ Automated Alignment:** Standardizes header identifiers and executes MAFFT multiple sequence alignment (`mafft --auto`).
- **🌲 Maximum Likelihood Phylogenetics:** Reconstructs ML trees using FastTree (GTR+CAT nucleotide substitution model).
- **⏳ Molecular Clock Rate-Smoothing:** Computes ultrametric Time Trees (chronograms) using `ape::chronos` relaxed molecular clock rate smoothing.
- **Outgroup Rooting:** Automatically identifies outgroup sequences by keyword (e.g., `OUTGROUP`) and roots the tree.
- **📊 Pairwise Distance Metrics:** Computes complete nucleotide p-distance matrices and group-level divergence statistics.
- **🎨 Publication Figures:** Renders a 4-panel publication figure in vector (`.svg`) and raster (`.png`) formats.
- **📄 Summary Report:** Generates an executive Word (`.docx`) report containing figures and divergence tables.
- **📓 Jupyter Notebook Included:** Run interactively step-by-step using [`AutoTimeTree_Pipeline.ipynb`](file:///Volumes/T-omics/usama_time_tree/AutoTimeTree/AutoTimeTree_Pipeline.ipynb).

---

## 📁 Repository Structure

```text
AutoTimeTree/
├── run_pipeline.sh             # 🚀 Main execution runner script
├── AutoTimeTree_Pipeline.ipynb  # 📓 Interactive Jupyter Notebook
├── README.md                   # 📄 Software documentation
├── .gitignore                  # 🙈 Git ignore configuration
├── data/                       # 📂 Input FASTA datasets
│   └── example_input.fasta     # Sample input sequence dataset
└── scripts/                    # 📁 Sub-module processing scripts
    ├── 01_prepare_and_align.py  # Step 1: Cleaning & MAFFT alignment
    ├── 02_build_ml_and_timetree.R   # Step 2: FastTree ML & ape chronogram
    ├── 03_analyze_and_visualize.py# Step 3: Distance calculation & plots
    └── 04_generate_report.py     # Step 4: Word doc report compiler
```

---

## ⚙️ Prerequisites & Dependencies

### System Dependencies
- **MAFFT** (`mafft`) — [MAFFT Installation Guide](https://mafft.cbrc.jp/alignment/software/)
- **FastTree** (`fasttree`) — [FastTree Installation Guide](http://www.microbesonline.org/fasttree/)
- **R (>= 4.0)** with `ape` package:
  ```r
  install.packages("ape")
  ```

### Python Dependencies (Python 3.8+)
Install required Python packages via `pip`:
```bash
pip install biopython pandas numpy matplotlib seaborn python-docx jupyter
```

---

## 🚀 Quick Start

### Option 1: Command Line Interface (CLI)
```bash
chmod +x run_pipeline.sh
./run_pipeline.sh data/example_input.fasta OUTGROUP
```

### Option 2: Jupyter Notebook
Open [`AutoTimeTree_Pipeline.ipynb`](file:///Volumes/T-omics/usama_time_tree/AutoTimeTree/AutoTimeTree_Pipeline.ipynb) in VS Code, Jupyter Lab, or Google Colab and execute the cells sequentially!

---

## 📄 License

Distributed under the MIT License. See `LICENSE` for more information.
