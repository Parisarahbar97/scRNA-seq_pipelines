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
import pandas as pd
import matplotlib.pyplot as plt

# Load the metadata
metadata_path = "/rds/general/user/pr422/projects/puklandmarkproject/live/Users/Parisa/ephemeral_files/Parisa_scdownstream/output_full/run3_celltypist_new_model/finalized/merged_metadata.csv"
df = pd.read_csv(metadata_path)

# Count number of cells per sample
cell_counts = df['sample'].value_counts().sort_index()

# Plot
plt.figure(figsize=(14, 6))
bars = plt.bar(cell_counts.index, cell_counts.values)
plt.xticks(rotation=90, ha='center')  
plt.ylabel("Number of cells")
plt.xlabel("Sample pool")
plt.title("Number of cells per sample (from finalized metadata)")
plt.tight_layout()
plt.show()

# Save as PNG and PDF
plt.savefig("/rds/general/user/pr422/projects/puklandmarkproject/live/Users/Parisa/ephemeral_files/Parisa_scdownstream/output_full/run3_celltypist_new_model/finalized/cell_counts_per_sample.png", dpi=300)
plt.savefig("/rds/general/user/pr422/projects/puklandmarkproject/live/Users/Parisa/ephemeral_files/Parisa_scdownstream/output_full/run3_celltypist_new_model/finalized/cell_counts_per_sample.pdf")  

plt.show()