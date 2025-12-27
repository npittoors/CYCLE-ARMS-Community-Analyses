import pandas as pd
import os
import re
import numpy as np

# Define file paths
coral_net_file = "./ALL_ARMS_CoralNet_Annotations.csv"
metadata_file = "./CYCLE_ARMS_metadata_FINAL.csv"
output_dir = "./CoralNet/"

# Load CSV files
coral_df = pd.read_csv(coral_net_file)
metadata_df = pd.read_csv(metadata_file)

# Create mapping for Label_code to Code_deff and Phylum
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
                    # Create sample_id_r
                    sample_id_r = full_match + "_CN"
                    return sample_id_r, site_locality
    
    # If we can't extract the information, return empty strings
    return "", ""

# Function to extract plate ID from name
def extract_plate_id(name):
    # Example: CYCLE_2021_ARMS_01_DIAback_P9T_stitched_CNedit.jpg
    match = re.search(r'_P(\d+[TB])_', name)
    if match:
        return match.group(1)  # Returns plate number like "9T"
    return ""

# Apply the extraction functions to the Name column
coral_df[['sample_id_r', 'sitelocality']] = pd.DataFrame(
    coral_df['Name'].apply(extract_info_from_name).tolist(), 
    index=coral_df.index
)
coral_df['plate_id'] = coral_df['Name'].apply(extract_plate_id)

# Map Label_code to Phylum
coral_df['Phylum'] = coral_df['Label_code'].map(lambda x: label_mapping.get(x, ("", ""))[1])

# Check for unmapped Label_code values
unique_labels = coral_df['Label_code'].unique()
missing_labels = [label for label in unique_labels if label not in label_mapping]
if missing_labels:
    print(f"WARNING: The following Label_code values have no mapping: {missing_labels}")

# Create a full sample identifier combining ARMS ID and plate
coral_df['full_sample_id'] = coral_df['sample_id_r'].str.replace('_CN', '') + '_P' + coral_df['plate_id']

# Count occurrences of each Phylum per plate
print("Calculating phylum counts per plate...")
phylum_counts = coral_df.groupby(['full_sample_id', 'Phylum']).size().unstack(fill_value=0)

# Calculate proportions (relative abundance)
print("Calculating phylum proportions per plate...")
phylum_props = phylum_counts.div(phylum_counts.sum(axis=1), axis=0)

# Create R-ready data format for PERMANOVA and NMDS
print("Preparing data for R analysis...")

# 1. Create a community data matrix (samples as rows, phyla as columns)
community_matrix = phylum_counts.copy()

# 2. Create a sample metadata dataframe
sample_metadata = pd.DataFrame(index=community_matrix.index)
sample_metadata['SampleID'] = sample_metadata.index
sample_metadata['Site'] = coral_df.drop_duplicates('full_sample_id').set_index('full_sample_id')['sitelocality'].to_dict()

# Print available columns in metadata_df for debugging
print("Available columns in metadata file:", metadata_df.columns.tolist())

# Identify site column in metadata
site_col = None
potential_site_columns = ['locality', 'site', 'sitelocality', 'location']
for col in potential_site_columns:
    if col in metadata_df.columns:
        if any(site in metadata_df[col].values for site in site_mapping.values()):
            site_col = col
            break

if site_col:
    print(f"Using column '{site_col}' from metadata file to match sites")
    
    # Create a dictionary to map from site to environmental values
    env_cols = ['trubidity_std_rank', 'latitude', 'longitude', 'temp_mean', 
                'depth', 'temp', 'salinity', 'turbidity_m']
    
    site_env_data = {}
    for col in env_cols:
        if col in metadata_df.columns:
            # Create a mapping from site to environmental value
            site_env_data[col] = {}
            for _, row in metadata_df.iterrows():
                site = row[site_col]
                if pd.notna(site) and pd.notna(row.get(col, np.nan)):
                    site_env_data[col][site] = row[col]
            
            # Add environmental data to sample metadata
            sample_metadata[col] = sample_metadata['Site'].map(site_env_data[col])
            print(f"Added column '{col}' to sample metadata")
else:
    print("Could not find a suitable column in metadata file for matching sites")
    print("Available columns:", metadata_df.columns.tolist())

# Save files in R-ready format
r_data_file = os.path.join(output_dir, "CoralNet_R_community_matrix_ALL.csv")
r_metadata_file = os.path.join(output_dir, "CoralNet_R_sample_metadata_ALL.csv")

community_matrix.to_csv(r_data_file)
sample_metadata.to_csv(r_metadata_file)

print(f"R-ready community matrix saved to: {r_data_file}")
print(f"R-ready sample metadata saved to: {r_metadata_file}")
print("Sample metadata columns:", sample_metadata.columns.tolist())
