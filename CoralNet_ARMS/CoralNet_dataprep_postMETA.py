# Python script to process CoralNet Output for downstream analysis
# Import and File Loading:

# 1.Import the necessary libraries (pandas, os, and re for regular expressions)
# - load both CSV files using pandas

# 2.Mapping Dictionaries:
# - label_mapping: Maps Label_code to Code_deff and Phylum values
# - site_mapping: Maps site codes to full site locality names

# 3. Extraction Function:
# - extract_info_from_name() parses the Name column to extract sample_id_r and sitelocality
# - It uses regex to find patterns like "01_DIAback" in names like "CYCLE_2021_ARMS_01_DIAback_P9T_stitched_CNedit.jpg"
# - It then identifies the site code and maps it to the full site locality

# 4. Column Population:
# - apply the extraction function to fill sample_id_r and sitelocality
# - map Label_code values to get Code_deff and Phylum
# - identify and copy over the environmental columns from the metadata file

# Output: The updated dataframe is saved to a new CSV file

import pandas as pd
import os
import re

# Define file paths
coral_net_file = "./ALL_ARMS_CoralNet_Annotations.csv"
metadata_file = "./CoralNet_R_sample_metadata.csv"

# Load CSV files
coral_df = pd.read_csv(coral_net_file)
metadata_df = pd.read_csv(metadata_file)

# Create mappings for Label_code to Code_deff and Phylum
label_mapping = {
    "_CCA": ("Rhodophyta*CCA", "Rhodophyta"),
    "_BSED": ("Bound*Sediment", "Sediment"),
    "_UNAV": ("Unavailble", "Unavailable"),
    "Bact": ("Bacterial_Biofilm", "Biofilm"),
    "_GRAL": ("Chlorophyta", "Chlorophyta"),
    "_RDAL": ("Rhodophyta", "Rhodophyta"),
    "_BRYN": ("Bryozoa*encrusting", "Bryozoa"),
    "_BRY": ("Bryozoa*erect", "Bryozoa"),
    "_RDEN": ("Rhodophyta*encrusting", "Rhodophyta"),
    "_CAWT": ("Calcareous*worm_tube", "Annelida"),
    "_FORM": ("Forminifera", "Forminifera"),
    "_RDUP": ("Rhodophyta*upright", "Rhodophyta"),
    "_HYD": ("Hydrozoa", "Hydrozoa"),
    "_BRAL": ("Ochrophyta", "Ochrophyta"),
    "_BRUP": ("Ochrophyta*upright", "Ochrophyta"),
    "_NR": ("No*recruitment", "No_recruitment"),
    "_SOWT": ("Soft*wormtube", "Annelida"),
    "_SP": ("Porifera", "Porifera"),
    "_SED": ("Sediment", "Sediment"),
    "_OMOL": ("Other*mollusc", "Mollusca"),
    "_BI": ("Bivalia*empty", "Mollusca"),
    "Film": ("Clear_biofilm", "Biofilm"),
    "_GREN": ("Chlorophyta*encrusting", "Chlorophyta"),
    "_UNK": ("Unkown", "Unkown"),
    "_TUNC": ("Ascidiacea*colonial", "Chordata"),
    "_GRUP": ("Chlorophyta*upright", "Chlorophyta"),
    "_TUNS": ("Ascidiacea*solitary", "Chordata"),
    "_BIEM": ("Bivalia*empty", "Mollusca"),
    "_BREN": ("Ochrophyta*encrusting", "Ochrophyta"),
    "_PTRO": ("Pterobranch", "Hemichordata"),
    "_GAS": ("Gastropoda", "Mollusca"),
    "_MOBF": ("Other*mobile_fauna", "Other_mobile_fauna")
}

# Create mapping for site locality (keeping the same as before)
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

# Function to extract sample_id_r and site_locality from Name
def extract_info_from_name(name):
    # Example: CYCLE_2021_ARMS_01_DIAback_P9T_stitched_CNedit.jpg
    # Looking for pattern like 01_DIAback in the name
    match = re.search(r'ARMS_(\d+_\w+)_', name)
    if match:
        full_match = match.group(1)  # e.g., 01_DIAback
        
        # Extract site code
        site_match = re.search(r'_([A-Za-z]+)', full_match)
        if site_match:
            site_code = site_match.group(1)
            # Check all known site codes to find the match (case-insensitive)
            for known_code in site_mapping.keys():
                if known_code.lower() in site_code.lower():
                    site_locality = site_mapping[known_code]
                    # Create sample_id_r - matching the format in the new metadata file
                    sample_id_r = full_match + "_CN"
                    return sample_id_r, site_locality
    
    # If we can't extract the information, return empty strings
    return "", ""

# Apply the extraction function to the Name column
coral_df[['sample_id_r', 'sitelocality']] = pd.DataFrame(
    coral_df['Name'].apply(extract_info_from_name).tolist(), 
    index=coral_df.index
)

# Map Label_code to Code_deff and Phylum
coral_df['Code_deff'] = coral_df['Label_code'].map(lambda x: label_mapping.get(x, ("", ""))[0])
coral_df['Phylum'] = coral_df['Label_code'].map(lambda x: label_mapping.get(x, ("", ""))[1])

# Environmental columns to copy from the new metadata file
# Based on the structure shown in the second document
env_columns = [
    'turbidity_std_rank', 'latitude', 'longitude', 'temp_mean', 'depth', 
    'temp', 'salinity', 'turbidity_m'
]

# Print information about the metadata file
print("Columns in metadata file:", metadata_df.columns.tolist())
print("Unique sites in metadata file:", metadata_df['Site'].unique())
print("Sample of full_sample_id values:", metadata_df['full_sample_id'].head().tolist())

# Create a mapping from site locality to environmental data
# Using the 'Site' column from the metadata file which matches our site_mapping values
site_to_env_data = {}

# Group by Site to get unique environmental values for each site
for site in metadata_df['Site'].unique():
    site_data = metadata_df[metadata_df['Site'] == site].iloc[0]  # Take first occurrence
    site_to_env_data[site] = {}
    for env_col in env_columns:
        if env_col in metadata_df.columns:
            site_to_env_data[site][env_col] = site_data[env_col]

# Apply the environmental data to the coral dataframe
for env_col in env_columns:
    if env_col in metadata_df.columns:
        coral_df[env_col] = coral_df['sitelocality'].apply(
            lambda x: site_to_env_data.get(x, {}).get(env_col, None)
        )

# Print some summary statistics to verify
print("\nSummary of filled data:")
print(f"Total rows in coral_df: {len(coral_df)}")
print(f"Rows with sample_id_r: {coral_df['sample_id_r'].notna().sum()}")
print(f"Rows with sitelocality: {coral_df['sitelocality'].notna().sum()}")
print(f"Rows with Code_deff: {coral_df['Code_deff'].notna().sum()}")
print(f"Rows with Phylum: {coral_df['Phylum'].notna().sum()}")

# Check unique sitelocality values found
print(f"Unique sitelocality values found: {coral_df['sitelocality'].unique()}")

# For each environmental column we added, print how many rows have data
for env_col in env_columns:
    if env_col in coral_df.columns:
        non_null_count = coral_df[env_col].notna().sum()
        print(f"Rows with {env_col}: {non_null_count}")

# Print a sample of the merged data to verify
print("\nSample of merged data:")
sample_cols = ['Name', 'sample_id_r', 'sitelocality', 'Label_code', 'Code_deff', 'Phylum', 'temp', 'depth', 'latitude']
available_cols = [col for col in sample_cols if col in coral_df.columns]
print(coral_df[available_cols].head())

# Save the updated dataframe
output_file = "/Users/nicolepittoors/Documents/Bioinformatics/CoralNet/CoralNet_AllPlates_Counts.csv"
coral_df.to_csv(output_file, index=False)

print(f"\nUpdated CSV saved to {output_file}")
print(f"Added columns: sample_id_r, sitelocality, Code_deff, Phylum and environmental variables")
