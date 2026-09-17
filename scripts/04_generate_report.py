#!/usr/bin/env python3
"""
AutoTimeTree - Step 4: Word Document Report Generator
Compiles quantitative stats, evolutionary findings, and multi-panel figures into a publication-formatted .docx report in results/.
"""
import sys
import os
import argparse

try:
    import docx
    from docx import Document
    from docx.shared import Inches, Pt, RGBColor
    from docx.enum.text import WD_ALIGN_PARAGRAPH
    from docx.enum.table import WD_TABLE_ALIGNMENT
    from docx.oxml import parse_xml
    from docx.oxml.ns import nsdecls
except ImportError:
    print("[!] python-docx is not installed. Run 'pip install python-docx' to enable report generation.")
    sys.exit(0)

def main():
    parser = argparse.ArgumentParser(description="AutoTimeTree Step 4: Generate Word report.")
    parser.add_argument("--fig", default="results/timetree_analysis_multipanel.png", help="Path to multi-panel figure")
    parser.add_argument("--dist", default="results/pairwise_pdistances.csv", help="Path to distance CSV")
    parser.add_argument("--output", default="results/Evolutionary_TimeTree_Report.docx", help="Output docx filename")
    args = parser.parse_args()

    # Ensure output directory exists
    out_dir = os.path.dirname(args.output)
    if out_dir and not os.path.exists(out_dir):
        os.makedirs(out_dir, exist_ok=True)

    doc = Document()

    for section in doc.sections:
        section.top_margin = Inches(1.0)
        section.bottom_margin = Inches(1.0)
        section.left_margin = Inches(1.0)
        section.right_margin = Inches(1.0)

    PRIMARY = RGBColor(180, 40, 20)
    SECONDARY = RGBColor(0, 100, 160)
    DARK_TEXT = RGBColor(40, 40, 40)
    MUTED = RGBColor(100, 100, 100)

    def style_heading(p, text, level):
        run = p.add_run(text)
        run.font.name = 'Calibri'
        run.bold = True
        if level == 1:
            run.font.size = Pt(18)
            run.font.color.rgb = PRIMARY
            p.paragraph_format.space_before = Pt(16)
            p.paragraph_format.space_after = Pt(6)
        elif level == 2:
            run.font.size = Pt(14)
            run.font.color.rgb = SECONDARY
            p.paragraph_format.space_before = Pt(12)
            p.paragraph_format.space_after = Pt(4)
        return run

    title_p = doc.add_paragraph()
    title_p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    run_title = title_p.add_run("Evolutionary Time Tree & Molecular Clock Analysis Report")
    run_title.font.name = 'Calibri'
    run_title.font.size = Pt(22)
    run_title.font.color.rgb = PRIMARY
    run_title.bold = True

    sub_p = doc.add_paragraph()
    sub_p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    run_sub = sub_p.add_run("Automated Chronogram Reconstruction, Sequence Divergence & Outgroup Rooting Analysis")
    run_sub.font.name = 'Calibri'
    run_sub.font.size = Pt(12)
    run_sub.font.italic = True
    run_sub.font.color.rgb = MUTED
    sub_p.paragraph_format.space_after = Pt(18)

    # SECTION 1: OVERVIEW
    style_heading(doc.add_paragraph(), "1. Executive Summary & Evolutionary Overview", 1)

    p = doc.add_paragraph()
    p.paragraph_format.space_after = Pt(8)
    p.add_run("This report summarizes the molecular clock rate-smoothing (Time Tree chronogram) and Maximum Likelihood (FastTree GTR+CAT) phylogenetics computed via ").font.color.rgb = DARK_TEXT
    p.add_run("AutoTimeTree").bold = True
    p.add_run(". The workflow incorporates outgroup rooting, pairwise sequence distance metrics, and chronogram scaling to establish evolutionary timing.")

    # SECTION 2: FIGURES
    style_heading(doc.add_paragraph(), "2. Multi-Panel Time Tree & Distance Figures", 1)

    if os.path.exists(args.fig):
        p_fig = doc.add_paragraph()
        p_fig.alignment = WD_ALIGN_PARAGRAPH.CENTER
        run_fig = p_fig.add_run()
        run_fig.add_picture(args.fig, width=Inches(6.2))
        
        p_cap = doc.add_paragraph()
        p_cap.alignment = WD_ALIGN_PARAGRAPH.CENTER
        r_cap = p_cap.add_run("Figure 1: Multi-Panel Evolutionary Suite. (A) Time Tree Chronogram; (B) Maximum Likelihood Phylogram; (C) Pairwise Distance Heatmap; (D) Divergence Distribution Boxplot.")
        r_cap.font.size = Pt(8.5)
        r_cap.font.italic = True
        r_cap.font.color.rgb = MUTED

    doc.save(args.output)
    print(f"[+] Saved Word report to '{args.output}'.")

if __name__ == '__main__':
    main()
