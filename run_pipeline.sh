#!/usr/bin/env bash
# ==============================================================================
# AutoTimeTree: Automated Molecular Clock & Time Tree Analysis Pipeline
# ==============================================================================

set -e

# Export PATH to include standard conda and system binary locations
export PATH="/opt/miniconda3/bin:/usr/local/bin:/opt/homebrew/bin:${PATH}"

INPUT_FASTA="${1:-data/example_input.fasta}"
OUTGROUP_KEY="${2:-OUTGROUP}"
RESULTS_DIR="results"

echo "========================================================================"
echo " AutoTimeTree: Automated Molecular Clock & Time Tree Pipeline"
echo "========================================================================"
echo " Input File:        ${INPUT_FASTA}"
echo " Outgroup Keyword:  ${OUTGROUP_KEY}"
echo " Results Directory: ${RESULTS_DIR}/"
echo "========================================================================"

mkdir -p "${RESULTS_DIR}"

chmod +x scripts/*.py 2>/dev/null || true
chmod +x scripts/*.R 2>/dev/null || true

# Step 1: Sequence Cleaning & MAFFT Alignment
echo ""
echo "[Step 1/4] Cleaning Sequence Headers and Running MAFFT Alignment..."
python3 scripts/01_prepare_and_align.py --input "${INPUT_FASTA}" --output "${RESULTS_DIR}/aligned_sequences.fasta" --cleaned "${RESULTS_DIR}/cleaned_sequences.fasta"

# Step 2: ML Tree & Chronogram Construction
echo ""
echo "[Step 2/4] Building Maximum Likelihood Tree and Time Tree Chronogram..."
Rscript scripts/02_build_ml_and_timetree.R "${RESULTS_DIR}/aligned_sequences.fasta" "${OUTGROUP_KEY}" "${RESULTS_DIR}/ml_tree.nwk" "${RESULTS_DIR}/rooted_ml_tree.nwk" "${RESULTS_DIR}/timetree.nwk"

# Step 3: Distance Matrix & Multi-Panel Visualization
echo ""
echo "[Step 3/4] Computing Distance Matrices and Rendering Multi-Panel Plots..."
python3 scripts/03_analyze_and_visualize.py --align "${RESULTS_DIR}/aligned_sequences.fasta" --tree "${RESULTS_DIR}/rooted_ml_tree.nwk" --timetree "${RESULTS_DIR}/timetree.nwk" --dist_out "${RESULTS_DIR}/pairwise_pdistances.csv" --fig_out "${RESULTS_DIR}/timetree_analysis_multipanel.png" --outgroup "${OUTGROUP_KEY}"

# Step 4: Word Document Report Generation
echo ""
echo "[Step 4/4] Generating Word Summary Report..."
python3 scripts/04_generate_report.py --fig "${RESULTS_DIR}/timetree_analysis_multipanel.png" --dist "${RESULTS_DIR}/pairwise_pdistances.csv" --output "${RESULTS_DIR}/Evolutionary_TimeTree_Report.docx"

echo ""
echo "========================================================================"
echo " Pipeline Finished Successfully!"
echo " Output Artifacts stored in ${RESULTS_DIR}/:"
echo "   - Time Tree Newick:         ${RESULTS_DIR}/timetree.nwk"
echo "   - Rooted ML Tree Newick:    ${RESULTS_DIR}/rooted_ml_tree.nwk"
echo "   - Multi-Panel Figure:       ${RESULTS_DIR}/timetree_analysis_multipanel.png and .svg"
echo "   - Pairwise Distance Matrix:  ${RESULTS_DIR}/pairwise_pdistances.csv"
echo "   - Summary Report:           ${RESULTS_DIR}/Evolutionary_TimeTree_Report.docx"
echo "========================================================================"
