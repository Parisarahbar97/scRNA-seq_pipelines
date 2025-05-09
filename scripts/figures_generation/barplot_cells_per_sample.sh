#First Load and check the metadata CSV file

cd /rds/general/user/pr422/projects/puklandmarkproject/live/Users/Parisa/ephemeral_files/Parisa_scdownstream/output_full/run3_celltypist_new_model/finalized

import pandas as pd

# Load the metadata CSV file
metadata = pd.read_csv("merged_metadata.csv")

# Display the first few rows to understand the structure
print(metadata.head())

# Identify the Sample Pool Column
print(metadata.columns)
# Index(['Unnamed: 0', 'label', 'n_counts', 'n_genes',
#       'celltypist:Human_AdultAged_Hippocampus:conf', 'batch',
#       'total_counts_mt', 'sample_original', 'total_counts',
#       'celltypist:Human_AdultAged_Hippocampus', 'n_genes_by_counts',
#       'pct_counts_mt', 'sample', 'scvi-global-0.5_leiden',
#       'scvi-unknown-0.5_leiden', 'scvi-global-1.0_leiden',
#       'scvi-unknown-1.0_leiden'],
#      dtype='object')

# Count Cells per Sample Pool:
cell_counts = metadata['sample'].value_counts().sort_index()
print(cell_counts)

# sample
# S10A_S18_mapped    23599
# S11A_S19_mapped    25161
# S11B_S20_mapped    21620
# S12B_S21_mapped    23713
# S13A_S22_mapped    20064
# S13B_S23_mapped    29992
# S14A_S24_mapped    19360
# S14B_S25_mapped    23012
# S15A_S26_mapped    10253
# S16A_S27_mapped    31692
# S16B_S28_mapped    33244
# S17B_S29_mapped    32156
# S18B_S30_mapped    21010
# S19A_S31_mapped      150
# S19B_S32_mapped     7432
# S1A_S1_mapped      27095
# S1B_S2_mapped      27191
# S2A_S3_mapped      31390
# S2B_S4_mapped      28544
# S3A_S5_mapped      21518
# S4A_S6_mapped      25840
# S4B_S7_mapped      35475
# S5A_S8_mapped      20959
# S5B_S9_mapped      23634
# ...
# S8B_S15_mapped     24982
# S9A_S16_mapped     22183
# S9B_S17_mapped     27721
# Name: count, dtype: int64

# Generating .png and.pdf plots/ final script :
# Add cell count labels on top of each bar
# Format the figure for publication (NatGen):
# - Width < 7 inches
# - Font size max 7 pt

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Load the metadata CSV
metadata = pd.read_csv("/rds/general/user/pr422/projects/puklandmarkproject/live/Users/Parisa/ephemeral_files/Parisa_scdownstream/output_full/run3_celltypist_new_model/finalized/merged_metadata.csv")  # Update the path as needed

# Count cells per sample
cell_counts = metadata['sample'].value_counts().sort_index()

# Set up the figure size (in inches)
fig, ax = plt.subplots(figsize=(6.8, 4))  # <7 inch width

# Set seaborn style
sns.set(style="whitegrid", font_scale=0.7)  # ~7pt font

# Create barplot
bars = sns.barplot(x=cell_counts.index, y=cell_counts.values, color="steelblue", edgecolor="black", ax=ax)

# Annotate bars with the cell count
for i, val in enumerate(cell_counts.values):
    ax.text(i, val + max(cell_counts.values) * 0.01, f'{val:,}', 
            ha='center', va='bottom', fontsize=6)

# Rotate x-axis labels
ax.set_xticklabels(ax.get_xticklabels(), rotation=90, fontsize=6)
ax.set_ylabel("Number of Cells", fontsize=7)
ax.set_xlabel("Sample", fontsize=7)
ax.tick_params(axis='y', labelsize=6)

# Tight layout for better spacing
plt.tight_layout()
plt.savefig("/rds/general/user/pr422/projects/puklandmarkproject/live/Users/Parisa/ephemeral_files/Parisa_scdownstream/output_full/run3_celltypist_new_model/finalized/cell_counts_per_sample.pdf", format='pdf', dpi=300)
plt.show() 