#Title: Chromosomal plot with hypergeometric testing for enrichment
#Description: Tests if the overlapping DE genes and TEs are enriched in genic regions of the genome

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.patches import Ellipse
import seaborn as sns
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests

# Set the colorblind-friendly color palette
colorblind_palette = sns.color_palette("colorblind")
diverging_colors = sns.color_palette("RdBu", 10)

# Read the matched CSV file
matched_df = pd.read_csv("FishTEA/matched_chromfilter.csv")

# Read the chrom_end.txt file and remove duplicate markers
end_data = []
prev_chrom = None
with open('chrom_end.txt') as f:
    for line in f:
        chrom, pos = line.split()
        if prev_chrom != chrom:
            end_data.append((chrom, int(pos)))
            prev_chrom = chrom

# Load significant regions
significant_regions = pd.read_csv("significant_te_de_gene_enrichment.csv")

# Create the plot
fig, ax = plt.subplots(figsize=(15, 8))

# Set colors based on TE_family column
colors = []
for family in matched_df['TE_family']:
    if 'SINE' in family or 'LINE' in family or 'LTR' in family:
        colors.append(diverging_colors[0])  # Blue
    elif 'DNA' in family or 'RC' in family:
        colors.append(diverging_colors[6])  # Red
    else:
        colors.append(diverging_colors[3])  # Green (default)

# Add scatter plot for TE locations with different colors
marker_sizes = []
for family in matched_df['TE_family']:
    if 'SINE' in family or 'LINE' in family or 'LTR' in family:
        marker_sizes.append(22)
    elif 'DNA' in family or 'RC' in family:
        marker_sizes.append(40)
    else:
        marker_sizes.append(0)

te_scatter = ax.scatter(matched_df['chromosome_TE'], matched_df['TE_start'], color=colors, alpha=0.7, s=marker_sizes, zorder=2)

# Add scatter plot for gene locations
ax.scatter(matched_df['Gene_chromosome'], matched_df['Gene_start'], label='Gene locations', color=diverging_colors[9], alpha=0.8, s=100)

# Add centromere markers
centromere_data = pd.read_csv('chrcen.txt', sep='\s+', header=None, names=['chromosome', 'cen_position'])
for _, row in centromere_data.iterrows():
    ax.plot(row['chromosome'], row['cen_position'], marker='D', markersize=5, color=colorblind_palette[7], zorder=2, alpha=1)

# Highlight significant regions with a purple marker slightly left of the chromosome
ax.scatter(significant_regions['chrom'] - 0.3, significant_regions['start'], color='orange', alpha=0.5, marker='D', label='Significant Regions', s=50, zorder=3)

# Set labels and title
plt.xlabel('Chromosome', fontsize=18, fontweight='bold')
plt.ylabel('Position', fontsize=18, fontweight='bold')
plt.title('RNA-seq: Significantly DE overlapping Genes and TEs', fontsize=20, fontweight='bold')

# Set x-axis tick positions and labels dynamically based on unique chromosome values
unique_chromosomes = sorted(matched_df['chromosome_TE'].unique())
ax.set_xticks(range(1, len(unique_chromosomes) + 1))
ax.set_xticklabels(unique_chromosomes, fontsize=14, fontweight='bold')

# Add ellipses for each half of the chromosome
for i, (chrom, end) in enumerate(end_data, start=1):
    centromere_position = centromere_data[centromere_data['chromosome'] == int(chrom)]['cen_position'].iloc[0]
    ax.add_patch(Ellipse((i, centromere_position / 2), 0.2, centromere_position, edgecolor='grey', linewidth=0.4, fill=False, zorder=0))
    ax.add_patch(Ellipse((i, (end + centromere_position) / 2), 0.2, end - centromere_position, edgecolor='grey', linewidth=0.3, fill=False, zorder=0))

plt.margins(y=0.02)
ax.set_xlim(0.4, len(unique_chromosomes) + 0.6)

# Create custom legend
legend_elements = [
    te_scatter,
    plt.Line2D([0], [0], marker='o', color='w', alpha=0.7, markerfacecolor=diverging_colors[0], markersize=8, label='DE Class II RNA TE'),
    plt.Line2D([0], [0], marker='o', color='w', alpha=0.7, markerfacecolor=diverging_colors[6], markersize=11, label='DE Class I DNA TE'),
    plt.Line2D([0], [0], marker='o', color='w', alpha=0.7, markerfacecolor=diverging_colors[9], markersize=15, label='DE Gene'),
    plt.Line2D([0], [0], marker='D', color='w', alpha=1, markersize=6, markerfacecolor=colorblind_palette[7], markeredgewidth=0, label='Centromere'),
    plt.Line2D([0], [0], marker='D', color='w', alpha=0.5, markersize=6, markerfacecolor='orange', label='Significant Regions')
]

# Add legend
ax.legend(handles=legend_elements, loc='upper right', fontsize=12)

# Show the plot
plt.tight_layout()
plt.xticks(fontsize=18, fontweight='bold')
plt.yticks(fontsize=18, fontweight='bold')
plt.savefig("te_genes_plot_phenogram_socstress.png", dpi=600)
plt.show()
