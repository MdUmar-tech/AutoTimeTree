#!/usr/bin/env python3
"""
AutoTimeTree - Step 3: Distance Matrix & Multi-Panel Visualization
Computes pairwise nucleotide p-distance matrix and renders publication-ready
multi-panel figures comparing ML Phylogram, Time Tree Chronogram, Distance Heatmap, and Group Divergence.
Outputs are saved to the results/ directory.
"""
import sys
import os
import argparse
import numpy as np
import pandas as pd
from Bio import SeqIO, Phylo
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec
import seaborn as sns

plt.style.use('seaborn-v0_8-whitegrid')
plt.rcParams['font.sans-serif'] = 'Helvetica'
plt.rcParams['axes.edgecolor'] = '#333333'
plt.rcParams['axes.linewidth'] = 0.8

def main():
    parser = argparse.ArgumentParser(description="AutoTimeTree Step 3: Distance calculation and visualization.")
    parser.add_argument("--align", default="results/aligned_sequences.fasta", help="Aligned FASTA file")
    parser.add_argument("--tree", default="results/rooted_ml_tree.nwk", help="Rooted ML tree Newick file")
    parser.add_argument("--timetree", default="results/timetree.nwk", help="Timetree Newick file")
    parser.add_argument("--dist_out", default="results/pairwise_pdistances.csv", help="CSV output for distance matrix")
    parser.add_argument("--fig_out", default="results/timetree_analysis_multipanel.png", help="PNG output for multipanel figure")
    parser.add_argument("--outgroup", default="OUTGROUP", help="Keyword identifying the outgroup sequence")
    args = parser.parse_args()

    if not os.path.exists(args.align):
        print(f"[Error] Alignment file '{args.align}' not found!")
        sys.exit(1)

    # Ensure output directory exists
    out_dir = os.path.dirname(args.fig_out)
    if out_dir and not os.path.exists(out_dir):
        os.makedirs(out_dir, exist_ok=True)

    print(f"[*] Reading alignment from '{args.align}'...")
    alignment = list(SeqIO.parse(args.align, 'fasta'))
    seq_dict = {rec.id: str(rec.seq).upper() for rec in alignment}

    taxa = list(seq_dict.keys())
    n_taxa = len(taxa)
    print(f"[+] Processed {n_taxa} taxa for distance calculation.")

    def get_group(taxa_name):
        if args.outgroup.lower() in taxa_name.lower() or 'xylella' in taxa_name.lower():
            return 'Outgroup'
        elif taxa_name.endswith('_P'):
            return 'Pathogenic'
        elif taxa_name.endswith('_NP'):
            return 'Non-Pathogenic'
        else:
            return 'Intermediate'

    group_map = {t: get_group(t) for t in taxa}

    color_map = {
        'Pathogenic': '#D9381E',
        'Non-Pathogenic': '#0072B2',
        'Outgroup': '#555555',
        'Intermediate': '#E69F00'
    }

    # Compute p-distance matrix
    dist_matrix = np.zeros((n_taxa, n_taxa))
    for i in range(n_taxa):
        for j in range(i+1, n_taxa):
            s1 = seq_dict[taxa[i]]
            s2 = seq_dict[taxa[j]]
            valid_pos = sum(1 for a, b in zip(s1, s2) if a in 'ATCG' and b in 'ATCG')
            diffs = sum(1 for a, b in zip(s1, s2) if a in 'ATCG' and b in 'ATCG' and a != b)
            p_dist = diffs / valid_pos if valid_pos > 0 else 0.0
            dist_matrix[i, j] = p_dist
            dist_matrix[j, i] = p_dist

    df_dist = pd.DataFrame(dist_matrix, index=taxa, columns=taxa)
    df_dist.to_csv(args.dist_out)
    print(f"[+] Saved distance matrix to '{args.dist_out}'.")

    if not (os.path.exists(args.tree) and os.path.exists(args.timetree)):
        print("[!] Tree files not found. Skipping plot generation.")
        return

    # Load trees
    ml_tree = Phylo.read(args.tree, 'newick')
    time_tree = Phylo.read(args.timetree, 'newick')

    outgroup_tips = [c for c in time_tree.get_terminals() if args.outgroup.lower() in c.name.lower() or 'xylella' in c.name.lower()]
    if outgroup_tips:
        time_tree.root_with_outgroup(outgroup_tips[0])

    def draw_tree_custom(tree, ax, title, is_timetree=False):
        terminals = tree.get_terminals()
        ml_order = [c.name for c in ml_tree.get_terminals()]
        terminals = sorted(terminals, key=lambda c: ml_order.index(c.name) if c.name in ml_order else 0)
        
        y_coords = {c: i for i, c in enumerate(terminals)}
        
        def get_x(clade):
            return tree.distance(tree.root, clade)

        def draw_clade(clade):
            x_parent = get_x(clade)
            if clade.is_terminal():
                y_parent = y_coords[clade]
                grp = group_map.get(clade.name, 'Intermediate')
                c_color = color_map[grp]
                ax.scatter(x_parent, y_parent, color=c_color, s=40, zorder=5)
                label_text = clade.name.replace('_OUTGROUP', '').replace('_P', ' (P)').replace('_NP', ' (NP)')
                ax.text(x_parent + (0.01 if not is_timetree else 0.02), y_parent, label_text, 
                        va='center', fontsize=7.5, fontweight='bold', color=c_color)
                return y_parent
            else:
                child_ys = [draw_clade(child) for child in clade]
                y_min, y_max = min(child_ys), max(child_ys)
                y_mid = sum(child_ys) / len(child_ys)
                
                ax.plot([x_parent, x_parent], [y_min, y_max], color='#444444', lw=1.2)
                
                for child in clade:
                    x_child = get_x(child)
                    y_child = child_ys[clade.clades.index(child)]
                    ax.plot([x_parent, x_child], [y_child, y_child], color='#444444', lw=1.2)
                    
                if not is_timetree and clade.confidence and clade.confidence > 0.5:
                    ax.text(x_parent, y_mid + 0.15, f'{int(clade.confidence*100)}', fontsize=6, color='#666666', ha='right')
                    
                return y_mid

        draw_clade(tree.root)
        ax.set_ylim(-1, len(terminals))
        ax.set_yticks([])
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.grid(False)
        ax.set_title(title, fontsize=12, fontweight='bold', pad=10)

    fig = plt.figure(figsize=(18, 14))
    gs = GridSpec(2, 2, figure=fig, width_ratios=[1.2, 1], height_ratios=[1, 0.8], hspace=0.3, wspace=0.25)

    ax_tt = fig.add_subplot(gs[0, 0])
    draw_tree_custom(time_tree, ax_tt, 'A. Time Tree Chronogram (Relaxed Clock)', is_timetree=True)
    ax_tt.set_xlabel('Relative Divergence Time (Arbitrary Time Units)', fontsize=10, fontweight='bold')

    ax_ml = fig.add_subplot(gs[0, 1])
    draw_tree_custom(ml_tree, ax_ml, 'B. Maximum Likelihood Phylogram (FastTree GTR+CAT)', is_timetree=False)
    ax_ml.set_xlabel('Nucleotide Substitutions per Site', fontsize=10, fontweight='bold')

    ax_hm = fig.add_subplot(gs[1, 0])
    ordered_taxa = sorted(taxa, key=lambda t: (group_map[t], t))
    df_ordered_dist = df_dist.loc[ordered_taxa, ordered_taxa]
    clean_taxa_names = [t.replace('_OUTGROUP', '').replace('_P', ' [P]').replace('_NP', ' [NP]') for t in ordered_taxa]

    sns.heatmap(df_ordered_dist, ax=ax_hm, cmap='YlOrRd', cbar_kws={'label': 'p-distance (substitutions/site)'},
                xticklabels=clean_taxa_names, yticklabels=clean_taxa_names)
    ax_hm.set_title('C. Pairwise Sequence Divergence Heatmap', fontsize=12, fontweight='bold', pad=10)
    ax_hm.tick_params(axis='both', which='major', labelsize=6)

    ax_box = fig.add_subplot(gs[1, 1])
    groups_present = set(group_map.values())
    
    data_dist = []
    if 'Pathogenic' in groups_present and 'Non-Pathogenic' in groups_present:
        p_taxa = [t for t in taxa if group_map[t] == 'Pathogenic']
        np_taxa = [t for t in taxa if group_map[t] == 'Non-Pathogenic']
        p_p_dists = [df_dist.loc[t1, t2] for i, t1 in enumerate(p_taxa) for t2 in p_taxa[i+1:]]
        np_np_dists = [df_dist.loc[t1, t2] for i, t1 in enumerate(np_taxa) for t2 in np_taxa[i+1:]]
        p_np_dists = [df_dist.loc[t1, t2] for t1 in p_taxa for t2 in np_taxa]
        
        data_dist.extend([{'Group': 'Pathogenic\n(Intra-group)', 'p-distance': d} for d in p_p_dists])
        data_dist.extend([{'Group': 'Non-Pathogenic\n(Intra-group)', 'p-distance': d} for d in np_np_dists])
        data_dist.extend([{'Group': 'Inter-group\n(Path vs Non-Path)', 'p-distance': d} for d in p_np_dists])

    out_taxa = [t for t in taxa if group_map[t] == 'Outgroup']
    ingroup_taxa = [t for t in taxa if group_map[t] != 'Outgroup']
    if out_taxa and ingroup_taxa:
        out_dists = [df_dist.loc[o, t] for o in out_taxa for t in ingroup_taxa]
        data_dist.extend([{'Group': 'Outgroup\n(vs Ingroup)', 'p-distance': d} for d in out_dists])

    if data_dist:
        df_box = pd.DataFrame(data_dist)
        sns.boxplot(data=df_box, x='Group', y='p-distance', hue='Group', ax=ax_box,
                    palette=['#D9381E', '#0072B2', '#8E44AD', '#555555'][:len(df_box['Group'].unique())],
                    boxprops=dict(alpha=0.8), width=0.5, showmeans=True, legend=False,
                    meanprops=dict(marker='o', markerfacecolor='white', markeredgecolor='black'))
        sns.stripplot(data=df_box, x='Group', y='p-distance', ax=ax_box, color='black', alpha=0.3, jitter=0.2, size=3)

    ax_box.set_title('D. Intra- vs. Inter-Group Divergence Distribution', fontsize=12, fontweight='bold', pad=10)
    ax_box.set_xlabel('')
    ax_box.set_ylabel('Pairwise p-distance', fontsize=10, fontweight='bold')

    legend_patches = [
        mpatches.Patch(color=color_map['Pathogenic'], label='Pathogenic (_P)'),
        mpatches.Patch(color=color_map['Non-Pathogenic'], label='Non-Pathogenic (_NP)'),
        mpatches.Patch(color=color_map['Intermediate'], label='Intermediate / Basal Lineages'),
        mpatches.Patch(color=color_map['Outgroup'], label='Outgroup Lineage')
    ]
    fig.legend(handles=legend_patches, loc='upper center', bbox_to_anchor=(0.5, 0.98), ncol=4, fontsize=11, frameon=True)

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    
    svg_out = os.path.splitext(args.fig_out)[0] + '.svg'
    plt.savefig(args.fig_out, dpi=300, bbox_inches='tight')
    plt.savefig(svg_out, bbox_inches='tight')
    print(f"[+] Saved multipanel figures to '{args.fig_out}' and '{svg_out}'.")

if __name__ == '__main__':
    main()
