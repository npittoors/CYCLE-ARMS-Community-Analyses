import pandas as pd
import re
import os
import numpy as np

# Define file paths
input_file = "Path/to/CYCLE_ARMS_2mm_counts.csv"
output_dir = "Your/output/path"

# Load the CSV data
print(f"Loading data from {input_file}...")
motile_df = pd.read_csv(input_file)

# Create mapping for site locality
site_mapping = {
    "DIAback": "Diaphus_background",
    "DIAcoral": "Diaphus_coral",
    "ALDcoral": "Alderdice_coral",
    "ALDback": "Alderdice_background",
    "ALDshal": "Alderdice_shallow",
    "MCG": "McGrail_coral",
    "BRIcoral": "Bright_coral",
    "BRIback": "Bright_background",
    "BRIshal": "Bright_shallow",
    "STE": "Stetson_coral",
    "EFGshal": "EFGB_shallow",
    "EFGdeep": "EFGB_deep"
}

# Define symbols for taxonomy forward filling
SYMBOLS = {
    'phylum': '',
    'class': '!',
    'order': '#',
    'family': '^',
    'genus': '*',
    'species': '**'
}

# Function to extract sample_id_r and site_locality from eventID
def extract_info_from_eventid(eventid):
    # Example: CYCLE_2021_ARMS_05_DIAcoral
    # Look for pattern like 05_DIAcoral in the name
    match = re.search(r'ARMS_(\d+_\w+)', eventid)
    if match:
        full_match = match.group(1)  # e.g., 05_DIAcoral
        
        # Extract site code
        site_match = re.search(r'_([A-Za-z]+)', full_match)
        if site_match:
            site_code = site_match.group(1)
            # Check all known site codes to find the match (case-insensitive)
            for known_code in site_mapping.keys():
                if known_code.lower() in site_code.lower():
                    site_locality = site_mapping[known_code]
                    # Create sample_id_r
                    sample_id_r = full_match + "_2mm"
                    return sample_id_r, site_locality
    
    # If we can't extract the information, return empty strings
    return "", ""

# Apply the extraction function to the eventID column
print("Extracting sample_id_r and sitelocality from eventID...")
motile_df[['sample_id_r', 'sitelocality']] = pd.DataFrame(
    motile_df['eventID'].apply(extract_info_from_eventid).tolist(), 
    index=motile_df.index
)

# Define taxonomic rank columns in order
tax_ranks = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'scientificName']

# Forward fill taxonomy with lowest taxonomic level available
print("Forward filling taxonomy...")
for i, row in motile_df.iterrows():
    last_valid_rank = None
    last_valid_value = None
    
    # Find the lowest available taxonomic rank
    for rank in tax_ranks:
        if pd.notna(row[rank]) and row[rank] != '':
            last_valid_rank = rank
            last_valid_value = row[rank]
    
    # Skip if no valid taxonomy found
    if last_valid_rank is None:
        continue
    
    # Forward fill empty higher ranks with the last valid one + symbol
    for rank in tax_ranks:
        if rank == last_valid_rank:
            continue  # Skip the rank that already has a value
        
        if pd.isna(row[rank]) or row[rank] == '':
            # Get the symbol for this rank
            symbol = SYMBOLS.get(rank, '')
            if rank == 'scientificName':  # Special case for scientific name
                symbol = SYMBOLS.get('species', '')
            
            # Fill with the last valid value + symbol
            motile_df.at[i, rank] = last_valid_value + symbol

# Create a taxonomy column that will be used for aggregation
# Use the scientificName column as it will have the most specific identification with appropriate symbol
motile_df['taxonomy'] = motile_df['scientificName']

# Check for NA values in the key columns and report
print(f"Rows with missing sample_id_r: {motile_df['sample_id_r'].isna().sum()}")
print(f"Rows with missing sitelocality: {motile_df['sitelocality'].isna().sum()}")
print(f"Rows with missing taxonomy: {motile_df['taxonomy'].isna().sum()}")

# Create count matrix for taxonomic units across samples
print("Creating count matrix...")
count_matrix = motile_df.groupby(['sample_id_r', 'taxonomy'])['individualCount'].sum().unstack(fill_value=0)

# Calculate proportions
print("Calculating proportions...")
prop_matrix = count_matrix.div(count_matrix.sum(axis=1), axis=0)

# Add site information to the matrices
sample_to_site = motile_df.drop_duplicates('sample_id_r').set_index('sample_id_r')['sitelocality'].to_dict()
site_info = pd.DataFrame({'sitelocality': [sample_to_site.get(idx, "") for idx in count_matrix.index]},
                        index=count_matrix.index)

# Combine with count and proportion matrices
count_with_site = pd.concat([site_info, count_matrix], axis=1)
prop_with_site = pd.concat([site_info, prop_matrix], axis=1)

# Create R-ready files
print("Preparing R-ready files...")
r_count_file = os.path.join(output_dir, "Motile_count_matrix.csv")
r_meta_file = os.path.join(output_dir, "Motile_metadata.csv")

# Save the count matrix for R
count_matrix.to_csv(r_count_file)

# Create and save metadata for R
metadata = pd.DataFrame(index=count_matrix.index)
metadata['SampleID'] = metadata.index
metadata['Site'] = [sample_to_site.get(idx, "") for idx in metadata.index]

# Save metadata
metadata.to_csv(r_meta_file)

# Save the processed dataframe with taxonomic information
processed_file = os.path.join(output_dir, "Motile_processed.csv")
motile_df.to_csv(processed_file, index=False)

# Save count and proportion matrices with site info
count_file = os.path.join(output_dir, "Motile_counts_with_site.csv")
prop_file = os.path.join(output_dir, "Motile_proportions_with_site.csv")
count_with_site.to_csv(count_file)
prop_with_site.to_csv(prop_file)

print(f"Processed data saved to {processed_file}")
print(f"Count matrix saved to {r_count_file}")
print(f"Metadata saved to {r_meta_file}")
print(f"Counts with site info saved to {count_file}")
print(f"Proportions with site info saved to {prop_file}")

# Aggregate counts by order
print("Creating order-level count matrix...")
# Extract order information for each taxonomic unit
order_info = motile_df.drop_duplicates('taxonomy')[['taxonomy', 'order']].set_index('taxonomy')['order'].to_dict()

# Create order-level count matrix
order_counts = pd.DataFrame(0, index=count_matrix.index, columns=sorted(motile_df['order'].unique()))
for taxa_col in count_matrix.columns:
    if taxa_col in order_info:
        order = order_info[taxa_col]
        if pd.notna(order) and order != '':
            order_counts[order] += count_matrix[taxa_col]

# Calculate order-level proportions
order_props = order_counts.div(order_counts.sum(axis=1), axis=0)

# Add site information
order_counts_with_site = pd.concat([site_info, order_counts], axis=1)
order_props_with_site = pd.concat([site_info, order_props], axis=1)

# Save order-level matrices
order_counts_file = os.path.join(output_dir, "Motile_order_counts_with_site.csv")
order_props_file = os.path.join(output_dir, "Motile_order_proportions_with_site.csv")
order_counts_with_site.to_csv(order_counts_file)
order_props_with_site.to_csv(order_props_file)

print(f"Order counts saved to: {order_counts_file}")
print(f"Order proportions saved to: {order_props_file}")
