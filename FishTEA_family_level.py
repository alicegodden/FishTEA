# Title: FishPi overlapping genes and TEs plot
# Subtitle: Plotting at the family level not class
# Author: Dr. Alice M. Godden

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.patches import Ellipse
import seaborn as sns

# Set the rocket color palette
rocket_palette = sns.color_palette("rocket", 6)

# Read the matched CSV file
matched_df = pd.read_csv("matched_Ov_tetrans.csv")

# Read the chrom_end.txt file and remove duplicate markers
end_data = []
prev_chrom = None
with open('chrom_end.txt') as f:
    for line in f:
        chrom, pos = line.split()
        if prev_chrom != chrom:
            end_data.append((chrom, int(pos)))
            prev_chrom = chrom

# Create the plot
fig, ax = plt.subplots(figsize=(15, 8))

# Set colors based on TE_family column
te_colors = {
    'DNA': rocket_palette[0],       # Color 1 from rocket palette
    'LINE': rocket_palette[1],      # Color 2 from rocket palette
    'LTR': rocket_palette[2],       # Color 3 from rocket palette
    'RC': rocket_palette[3],        # Color 4 from rocket palette
    'SATELLITE': rocket_palette[4], # Color 5 from rocket palette
    'Unknown': rocket_palette[5]    # Color 6 from rocket palette (default for unknown)
}

colors = [te_colors.get(family.split(' ')[0].upper(), te_colors['Unknown']) for family in matched_df['TE_family']]

# Add scatter plot for TE locations with different colors
te_scatter = ax.scatter(matched_df['chromosome_TE'], matched_df['TE_start'], color=colors, alpha=0.7, s=25, zorder=2)

# Add scatter plot for gene locations with a valid color and adjusted transparency
ax.scatter(matched_df['Gene_chromosome'], matched_df['Gene_start'], label='Gene locations', color=rocket_palette[5], alpha=0.8, s=85)

# Add centromere markers
centromere_data = pd.read_csv('chrcen.txt', sep='\s+', header=None, names=['chromosome', 'cen_position'])
for _, row in centromere_data.iterrows():
    ax.plot(row['chromosome'], row['cen_position'], marker='D', markersize=5, color='grey', zorder=2, alpha=1)

# Set labels and title
plt.xlabel('Chromosome', fontsize=18, fontweight='bold')
plt.ylabel('Position', fontsize=18, fontweight='bold')
plt.title('FishTEA: Significantly DE overlapping Genes and TEs Ovaries Temperature', fontsize=20, fontweight='bold')

# Set x-axis tick positions and labels dynamically based on unique chromosome values
unique_chromosomes = sorted(matched_df['chromosome_TE'].unique())
ax.set_xticks(range(1, len(unique_chromosomes) + 1))
ax.set_xticklabels(unique_chromosomes, fontsize=14, fontweight='bold')

# Add ellipses for each half of the chromosome
for i, (chrom, end) in enumerate(end_data, start=1):
    centromere_position = centromere_data[centromere_data['chromosome'] == int(chrom)]['cen_position'].iloc[0]

    # Add ellipse for the first half of the chromosome (from 0 to centromere)
    ax.add_patch(
        Ellipse((i, centromere_position / 2), 0.2, centromere_position, edgecolor='grey', linewidth=0.4, fill=False,
                zorder=0))

    # Add ellipse for the second half of the chromosome (from centromere to end)
    ax.add_patch(
        Ellipse((i, (end + centromere_position) / 2), 0.2, end - centromere_position, edgecolor='grey', linewidth=0.3,
                fill=False, zorder=0))

# Adjust the y-axis limit to add space under each chromosome plot
max_end_position = max([end for _, end in end_data])
ax.set_ylim(-0.01 * max_end_position, max_end_position*1.01)
ax.set_xlim(0.4, len(unique_chromosomes) + 0.6)

# Create custom legend
legend_elements = [te_scatter] + [
    plt.Line2D([0], [0], marker='o', color='w', alpha=0.7, markerfacecolor=color, markersize=10, label=te_type)
    for te_type, color in te_colors.items() if te_type != 'Unknown'
] + [
    plt.Line2D([0], [0], marker='o', color='w', alpha=0.7, markerfacecolor=rocket_palette[5], markersize=15, label='DE Gene'),
    plt.Line2D([0], [0], marker='D', color='w', alpha=1, markersize=6, markerfacecolor='grey', markeredgewidth=0, label='Centromere')
]

# Add legend with custom elements
ax.legend(handles=legend_elements, loc='upper right', fontsize=12)

# Show the plot
plt.tight_layout()
plt.xticks(fontsize=18, fontweight='bold')
plt.yticks(fontsize=18, fontweight='bold')
plt.savefig("te_genes_plot_phenogram_fishtea_family_Ov_temp.png", dpi=600)
plt.show()
