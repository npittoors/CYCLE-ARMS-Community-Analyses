#!/usr/bin/env python3
from Bio import SeqIO
import sys
import pandas as pd

fasta_file = sys.argv[1]  # Original ASV FASTA file
results_file = sys.argv[2]  # Final taxonomy results
output_file = sys.argv[3]  # Output file for ASVs with no hits

# Get all ASV IDs from FASTA
all_asvs = set()
seq_records = {}
for record in SeqIO.parse(fasta_file, "fasta"):
    all_asvs.add(record.id)
    seq_records[record.id] = record

print(f"Found {len(all_asvs)} total ASVs in FASTA file")

# Get ASVs with hits from results file
results_df = pd.read_csv(results_file, sep='\t')
hit_asvs = set(results_df['ASV_ID'])

print(f"Found {len(hit_asvs)} ASVs with taxonomy hits")

# Find ASVs with no hits
missing_asvs = all_asvs - hit_asvs
print(f"Found {len(missing_asvs)} ASVs with no taxonomy hits")

# Write the no-hit ASVs to a file
with open(output_file, 'w') as f:
    f.write("ASV_ID\n")
    for asv_id in sorted(missing_asvs):
        f.write(f"{asv_id}\n")

print(f"Saved {len(missing_asvs)} ASVs with no hits to {output_file}")

# Extract sequences for no-hit ASVs
fasta_output = output_file.replace('.tsv', '.fasta')
with open(fasta_output, 'w') as f:
    for asv_id in sorted(missing_asvs):
        if asv_id in seq_records:
            SeqIO.write(seq_records[asv_id], f, "fasta")

print(f"Saved sequences for {len(missing_asvs)} ASVs with no hits to {fasta_output}")
