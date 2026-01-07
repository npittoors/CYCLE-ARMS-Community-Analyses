#!/usr/bin/env python3
"""
Configurable ASV Taxonomy Processing Script

This script processes ASV taxonomy assignments from BLAST results to:
1. Count and extract ASVs with pident = 0 or < threshold
2. Remove ASVs with pident < threshold
3. Remove ASVs where genus is 'Vertebrata' and species is 'environmental'
4. Backfill missing taxonomy values using WoRMS or GBIF databases
5. Forward fill missing taxonomy ranks using symbols

Key features:
- Supports both COI and 18S markers with appropriate thresholds
- Configurable backfill source (WoRMS or GBIF)
- Compatible with both full and reduced database modes
- Produces separate file for ASVs with no hits meeting criteria

Usage:
python final_process_taxonomy.py input.tsv output.tsv [options]
"""

import argparse
import pandas as pd
import requests
import time
import sys
import logging
import re
import os
import json

# Optional import of pygbif - will be checked at runtime
try:
    from pygbif import species as gbif_species
    GBIF_AVAILABLE = True
except ImportError:
    GBIF_AVAILABLE = False

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Define metazoan phyla plus Rhodophyta that we want to keep
METAZOAN_PLUS_RHODOPHYTA = {
    'Annelida', 'Arthropoda', 'Brachiopoda', 'Bryozoa', 'Chaetognatha',
    'Chordata', 'Cnidaria', 'Ctenophora', 'Echinodermata', 'Entoprocta',
    'Metazoa', 'Mollusca', 'Nemertea', 'Platyhelminthes', 'Porifera',
    'Rotifera', 'Sipuncula', 'Rhodophyta'
}

# Add a flag for tracking changes
BACKFILL_CHANGES = []  # Will store information about taxonomy changes

# Global cache for taxonomy lookup results
TAXONOMY_CACHE = {}  # Maps query string to taxonomy result

# Marker-specific settings
MARKER_SETTINGS = {
    'COI': {
        'min_identity': 85,
        'min_alignment': 200
    },
    '18S': {
        'min_identity': 85,
        'min_alignment': 300
    }
}

# Taxonomy rank symbols for forward filling
SYMBOLS = {
    'phylum': '!',
    'class': '#',
    'order': '^',
    'family': '*',
    'genus': '**',
    'species': ''}

def clean_data(df):
    """
    Clean the input data by fixing formatting issues and inconsistencies.
    
    Args:
        df (pandas.DataFrame): Input dataframe
        
    Returns:
        pandas.DataFrame: Cleaned dataframe
    """
    # Make a copy to avoid modifying the original
    cleaned_df = df.copy()
    
    # Fix numeric columns that might have non-numeric characters
    numeric_cols = ['pident', 'length', 'qcovs', 'evalue', 'bitscore', 
                   'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send']
    
    for col in numeric_cols:
        if col in cleaned_df.columns:
            # Convert to string first to handle any potential non-numeric characters
            cleaned_df[col] = cleaned_df[col].astype(str).str.replace(',', '.')
            # Remove any non-numeric characters (except decimal points and scientific notation)
            cleaned_df[col] = cleaned_df[col].apply(lambda x: re.sub(r'[^\d.e-]', '', str(x)) if pd.notna(x) else x)
            # Convert to float
            cleaned_df[col] = pd.to_numeric(cleaned_df[col], errors='coerce')
    
    return cleaned_df


def query_worms(taxon_name):
    """
    Query the WoRMS database for complete taxonomy information.
    
    Args:
        taxon_name (str): Taxonomic name to query
        
    Returns:
        dict: Complete taxonomy information or None if not found
    """
    base_url = "https://www.marinespecies.org/rest/AphiaRecordsByName"
    params = {
        "scientificname": taxon_name,
        "like": "false",
        "marine_only": "true",
        "offset": 1
    }
    
    try:
        response = requests.get(base_url, params=params)
        if response.status_code == 200 and isinstance(response.json(), list) and response.json():
            aphia_id = response.json()[0].get("AphiaID")
        if aphia_id:
                # Get full classification
                classification_url = f"https://www.marinespecies.org/rest/AphiaClassificationByAphiaID/{aphia_id}"
                class_response = requests.get(classification_url)
                if class_response.status_code == 200:
                    return class_response.json()
        return None
    except Exception as e:
        logger.error(f"Error querying WoRMS: {e}")
        return None


def query_gbif(taxon_name):
    """
    Query the GBIF database for complete taxonomy information.
    
    Args:
        taxon_name (str): Taxonomic name to query
        
    Returns:
        dict: Complete taxonomy information or None if not found
    """
    if not GBIF_AVAILABLE:
        logger.error("pygbif is not available. Install with: pip install pygbif")
        return None
    
    try:
        # Query GBIF backbone taxonomy
        result = gbif_species.name_backbone(name=taxon_name, strict=False)
        
        # Check if we got a valid result
        if result and result.get('matchType') in ['EXACT', 'FUZZY']:
            return result
        return None
    except Exception as e:
        logger.error(f"Error querying GBIF: {e}")
        return None

def extract_taxonomy_from_worms(worms_data):
    """
    Extract standardized taxonomy from WoRMS API response.
    
    Args:
        worms_data (dict): WoRMS API response
        
    Returns:
        dict: Extracted taxonomy information
    """
    taxonomy = {
        'phylum': None,
        'class': None,
        'order': None,
        'family': None,
        'genus': None,
        'species': None
    }
    
    if not worms_data:
        return taxonomy
    
    # Process classification data
    try:
        # Handle current level classification
        if 'rank' in worms_data and 'scientificname' in worms_data:
            rank_lower = worms_data.get('rank', '').lower()
            
            # Special handling for species
            if rank_lower == 'species':
                taxonomy['species'] = worms_data.get('scientificname')
                # Extract genus from species name
                parts = worms_data.get('scientificname', '').split()
                if parts:
                    taxonomy['genus'] = parts[0]
            elif rank_lower in taxonomy:
                taxonomy[rank_lower] = worms_data.get('scientificname')
        
        # Process the classification hierarchy
        # Check children first (lower ranks)
        current = worms_data
        while current and 'child' in current and current['child']:
            current = current['child']
            rank_lower = current.get('rank', '').lower()
            if rank_lower in taxonomy:
                taxonomy[rank_lower] = current.get('scientificname')
        
        # Then check parents (higher ranks)
        current = worms_data
        while current and 'parent' in current and current['parent']:
            current = current['parent']
            rank_lower = current.get('rank', '').lower()
            if rank_lower in taxonomy and not taxonomy[rank_lower]:
                taxonomy[rank_lower] = current.get('scientificname')
    
    except Exception as e:
        logger.error(f"Error extracting taxonomy from WoRMS data: {e}")
    
    return taxonomy


def extract_taxonomy_from_gbif(gbif_data):
    """
    Extract standardized taxonomy from GBIF API response.
    
    Args:
        gbif_data (dict): GBIF API response
        
    Returns:
        dict: Extracted taxonomy information
    """
    taxonomy = {
        'phylum': None,
        'class': None,
        'order': None,
        'family': None,
        'genus': None,
        'species': None
    }
    
    if not gbif_data:
        return taxonomy
    
    # Process GBIF data
    try:
        # Map GBIF fields directly to our taxonomy dictionary
        for rank in taxonomy.keys():
            if rank in gbif_data and gbif_data[rank]:
                taxonomy[rank] = gbif_data[rank]
        
        # Special handling for species rank
        if 'species' in gbif_data and gbif_data['species']:
            taxonomy['species'] = gbif_data['species']
        elif 'scientificName' in gbif_data and gbif_data.get('rank', '').lower() == 'species':
            taxonomy['species'] = gbif_data['scientificName']
            
        # If we have genus in the response, make sure it's captured
        if 'genus' in gbif_data and gbif_data['genus']:
            taxonomy['genus'] = gbif_data['genus']
        # If we have a species but no genus, extract genus from species name
        elif taxonomy['species'] and not taxonomy['genus']:
            parts = taxonomy['species'].split()
            if parts:
                taxonomy['genus'] = parts[0]
                
    except Exception as e:
        logger.error(f"Error extracting taxonomy from GBIF data: {e}")
    
    return taxonomy

def backfill_taxonomy(asv_data, source='worms'):
    """
    Backfill missing taxonomy values using taxonomy database with caching.
    
    Args:
        asv_data (pandas.DataFrame): ASV data with taxonomy columns
        source (str): Source database to use ('worms' or 'gbif')
        
    Returns:
        pandas.DataFrame: Updated ASV data with backfilled taxonomy
    """
    if source not in ['worms', 'gbif']:
        logger.error(f"Invalid source '{source}'. Using 'worms' as fallback.")
        source = 'worms'
        
    if source == 'gbif' and not GBIF_AVAILABLE:
        logger.error("pygbif is not installed. Install with: pip install pygbif. Falling back to WoRMS.")
        source = 'worms'
    
    total_asvs = len(asv_data)
    processed = 0
    
    # Make a copy to avoid modifying the original during iteration
    result_df = asv_data.copy()
    
    # Track unique taxonomy queries to avoid redundant API calls
    queries_made = 0
    cache_hits = 0
    
    for idx, row in asv_data.iterrows():
        processed += 1
        if processed % 10 == 0:
            logger.info(f"Processed {processed}/{total_asvs} ASVs for {source.upper()} backfilling")
        
        # Check if we need backfilling (any taxonomy ranks are missing)
        needs_backfill = False
        for rank in ['phylum', 'class', 'order', 'family', 'genus', 'species']:
            if pd.isna(row[rank]) or not row[rank]:
                needs_backfill = True
                break
                
        if not needs_backfill:
            logger.debug(f"Skipping ASV {row['ASV_ID']} - all taxonomy ranks already present")
            continue
        
        # Determine the best query to use, in priority order
        query_value = None
        query_description = None
        
        # Option 1: scientific_name
        if pd.notna(row.get('scientific_name')) and row.get('scientific_name'):
            # Skip if scientific_name contains "environmental" or "sp."
            sci_name = str(row['scientific_name'])
            if 'environmental' not in sci_name.lower() and ' sp.' not in sci_name.lower():
                query_value = sci_name
                query_description = "scientific_name"
        
        # Option 2: genus + species
        if not query_value and pd.notna(row.get('genus')) and row.get('genus') and pd.notna(row.get('species')) and row.get('species'):
            # Skip if species is "environmental" or "sp."
            species = str(row['species'])
            if 'environmental' not in species.lower() and 'sp.' not in species.lower():
                query_value = f"{row['genus']} {species}"
                query_description = "genus+species"
        
        # Option 3: Lowest taxonomic rank (excluding species by itself, which is often not useful)
        if not query_value:
            for rank in ['genus', 'family', 'order', 'class', 'phylum']:
                if pd.notna(row[rank]) and row[rank]:
                    query_value = row[rank]
                    query_description = rank
                    break
        
        # Skip if we couldn't find a good query value
        if not query_value:
            continue
        
        # Check if we have this query in the cache
        cache_key = f"{source}:{query_value}"
        if cache_key in TAXONOMY_CACHE:
            taxonomy_data = TAXONOMY_CACHE[cache_key]
            cache_hits += 1
            logger.debug(f"Cache hit for {query_value} (using: {query_description})")
        else:
            # Query the appropriate database for taxonomy
            taxonomy_data = None
            if source == 'worms':
                logger.info(f"Querying WoRMS for {query_value} (using: {query_description})")
                db_response = query_worms(query_value)
                if db_response:
                    taxonomy_data = extract_taxonomy_from_worms(db_response)
            else:  # GBIF
                logger.info(f"Querying GBIF for {query_value} (using: {query_description})")
                db_response = query_gbif(query_value)
                if db_response:
                    taxonomy_data = extract_taxonomy_from_gbif(db_response)
            
            # Store in cache
            TAXONOMY_CACHE[cache_key] = taxonomy_data
            queries_made += 1
            
            # Be nice to the API (sleep only for WoRMS - GBIF is more robust)
            if source == 'worms':
                time.sleep(1)
        
        if taxonomy_data:
            # Only update values that are missing or empty
            for rank, value in taxonomy_data.items():
                # Only update if the current value is missing or empty and we have a value to use
                if value and (pd.isna(result_df.at[idx, rank]) or result_df.at[idx, rank] == ''):
                    old_value = result_df.at[idx, rank] if pd.notna(result_df.at[idx, rank]) else '(empty)'
                    logger.debug(f"  Filling {rank}: {value} for ASV {row['ASV_ID']}")
                    result_df.at[idx, rank] = value
                    
                    # Track this change
                    BACKFILL_CHANGES.append({
                        'ASV_ID': row['ASV_ID'],
                        'database': row.get('database', 'unknown'),
                        'rank': rank,
                        'old_value': old_value,
                        'new_value': value,
                        'query_value': query_value,
                        'source': source.upper()
                    })
    
    logger.info(f"Made {queries_made} unique taxonomy queries with {cache_hits} cache hits")
    return result_df

def forward_fill_taxonomy(asv_data):
    """
    Forward fill missing taxonomy ranks using the value from the lowest filled rank
    and append the symbol from one rank higher.
    
    Args:
        asv_data (pandas.DataFrame): ASV data with taxonomy columns
        
    Returns:
        pandas.DataFrame: Updated ASV data with forward-filled taxonomy
    """
    # Create a copy to avoid modifying the original
    result_df = asv_data.copy()
    
    # Ranks in order from highest to lowest
    ranks = ['phylum', 'class', 'order', 'family', 'genus', 'species']
    
    # Process each row
    for idx, row in result_df.iterrows():
        # Step 1: Find the lowest rank with actual data (not forward-filled)
        # and clear any existing forward-filled values
        lowest_rank = None
        lowest_idx = -1
        
        # First pass: clear any forward-filled values and find lowest genuine rank
        for i, rank in enumerate(ranks):
            value = str(row[rank]) if pd.notna(row[rank]) else ""
            # Check if this value contains any symbol (indicating it was forward-filled)
            if any(sym in value for sym in SYMBOLS.values() if sym):
                # Clear this value as it was forward-filled
                result_df.at[idx, rank] = None
            elif value:  # If it has a non-empty value that doesn't contain symbols
                lowest_rank = rank
                lowest_idx = i
        
        # Skip if no taxonomy data
        if lowest_rank is None:
            continue
        
        # Step 2: Get the value from the lowest filled rank
        value_to_use = row[lowest_rank]
        
        # Step 3: Get the symbol from one rank higher
        symbol_idx = max(0, lowest_idx - 1)  # Ensure we don't go below 0
        symbol_rank = ranks[symbol_idx]
        symbol_to_use = SYMBOLS[symbol_rank]
        
        # Step 4: Fill all empty ranks below the lowest filled rank
        for i in range(lowest_idx + 1, len(ranks)):
            rank_to_fill = ranks[i]
            # Only fill if empty
            if pd.isna(result_df.at[idx, rank_to_fill]) or not result_df.at[idx, rank_to_fill]:
                result_df.at[idx, rank_to_fill] = f"{value_to_use}{symbol_to_use}"
    
    return result_df

def separate_taxa_by_phylum(data, keep_phyla):
    """
    Separate ASVs into two groups based on phylum
    
    Args:
        data (pandas.DataFrame): ASV data with taxonomy columns
        keep_phyla (set): Set of phyla to keep in the main output
        
    Returns:
        tuple: (kept_data, separated_data) - DataFrames with kept and separated ASVs
    """
    # Handle case sensitivity by converting to lowercase for comparison
    phyla_lower = {p.lower() for p in keep_phyla}
    
    # Create mask for ASVs to keep
    mask = data['phylum'].str.lower().isin(phyla_lower)
    
    # Split the data
    kept_data = data[mask].copy()
    separated_data = data[~mask].copy()
    
    return kept_data, separated_data

def main():
    parser = argparse.ArgumentParser(description='Process ASV taxonomy assignments.')
    parser.add_argument('--cache-file', help='File to load/save taxonomy cache')
    parser.add_argument('input_file', help='Input TSV file with ASV data')
    parser.add_argument("output_file", help="Base output file name (without marker suffix)")
    parser.add_argument('--marker', choices=['COI', '18S'],
                        help='Marker type: COI or 18S (required)')
    parser.add_argument('--min-pident', type=float, default=None,
                        help='Minimum percent identity threshold (default: depends on marker)')
    parser.add_argument('--min-alignment', type=int, default=None,
                        help='Minimum alignment length (default: depends on marker)')
    parser.add_argument('--no-backfill', action='store_true', 
                        help='Skip taxonomy database queries entirely')
    parser.add_argument('--source', choices=['worms', 'gbif'], default='gbif', 
                        help='Source database for taxonomy backfilling (default: gbif)')
    parser.add_argument('--only-fill-missing', action='store_true', 
                        help='Only fill in missing values, do not replace existing ones', default=True)
    parser.add_argument('--database-filter', 
                        help='Only process ASVs from a specific database (e.g., NCBI, MZG)', default=None)
    parser.add_argument('--changes-log', 
                        help='Output file to log taxonomy changes', default='taxonomy_changes.log')
    parser.add_argument('--no-hits-file', 
                        help='Output file for ASVs with no hits meeting criteria (optional)')
    parser.add_argument('--separate-by-phylum', action='store_true',
                    help='Separate metazoan/Rhodophyta from other phyla into separate files')
    args = parser.parse_args()
    


    # Set marker-specific settings if not provided
    marker_settings = MARKER_SETTINGS.get(args.marker, MARKER_SETTINGS['COI'])
    min_identity = args.min_pident if args.min_pident is not None else marker_settings['min_identity']
    min_alignment = args.min_alignment if args.min_alignment is not None else marker_settings['min_alignment']
    
    logger.info(f"Processing taxonomy for {args.marker} marker")
    logger.info(f"Using thresholds: min_identity={min_identity}%, min_alignment={min_alignment}bp")
    logger.info(f"Backfill source: {'None' if args.no_backfill else args.source.upper()}")
    
    try:
        # Load the data
        logger.info(f"Loading data from {args.input_file}")
        asv_data = pd.read_csv(args.input_file, sep='\t')
        
        # Clean the data
        logger.info("Cleaning input data")
        asv_data = clean_data(asv_data)
        
        # 1. Count ASVs with pident = 0 or < min_identity
        zero_pident = asv_data[asv_data['pident'] == 0]
        low_pident = asv_data[(asv_data['pident'] > 0) & (asv_data['pident'] < min_identity)]
        
        logger.info(f"ASVs with pident = 0: {len(zero_pident)}")
        logger.info(f"ASVs with 0 < pident < {min_identity}: {len(low_pident)}")
        
        # Save these to a file if no_hits_file is specified
        below_threshold = pd.concat([zero_pident, low_pident])
        if args.no_hits_file:
            below_threshold[['ASV_ID', 'pident']].to_csv(args.no_hits_file, sep='\t', index=False)
            logger.info(f"Saved {len(below_threshold)} ASVs with pident < {min_identity} to {args.no_hits_file}")
        
        # 2. Remove ASVs with pident < min_identity
        filtered_data = asv_data[asv_data['pident'] >= min_identity].copy()
        logger.info(f"Removed {len(asv_data) - len(filtered_data)} ASVs with pident < {min_identity}")
        logger.info(f"Remaining ASVs: {len(filtered_data)}")
        
        # 3. Remove ASVs with genus 'Vertebrata'/'Metazoa' and species 'environmental'
        env_seqs = filtered_data[
            ((filtered_data['genus'] == 'Vertebrata') & (filtered_data['species'] == 'environmental')) |
            ((filtered_data['genus'] == 'Metazoa') & (filtered_data['species'] == 'environmental'))
        ]
        logger.info(f"ASVs with genus 'Vertebrata'/'Metazoa' and species 'environmental': {len(env_seqs)}")
        
        filtered_data = filtered_data.drop(env_seqs.index)
        logger.info(f"Removed {len(env_seqs)} environmental sequences")
        logger.info(f"Remaining ASVs: {len(filtered_data)}")
        
        # Filter by database if specified
        if args.database_filter:
            pre_filter_count = len(filtered_data)
            filtered_data = filtered_data[filtered_data['database'] == args.database_filter]
            logger.info(f"Filtered to only include {args.database_filter} database: {len(filtered_data)}/{pre_filter_count} ASVs remaining")
        
        # 4. Backfill missing taxonomy using selected database
        def save_taxonomy_cache(filename):
            """Save taxonomy cache to file"""
            with open(filename, 'w') as f:
                json.dump(TAXONOMY_CACHE, f)
            logger.info(f"Saved {len(TAXONOMY_CACHE)} taxonomy entries to cache file {filename}")

        def load_taxonomy_cache(filename):
            """Load taxonomy cache from file"""
            global TAXONOMY_CACHE
            if os.path.exists(filename):
                with open(filename, 'r') as f:
                    TAXONOMY_CACHE = json.load(f)
                logger.info(f"Loaded {len(TAXONOMY_CACHE)} taxonomy entries from cache file {filename}")
            else:
                logger.info(f"No cache file found at {filename}, starting with empty cache")
        
        if not args.no_backfill:
            logger.info(f"Backfilling missing taxonomy using {args.source.upper()} database...")
            
            # Check if GBIF is available when requested
            if args.source == 'gbif' and not GBIF_AVAILABLE:
                logger.warning("pygbif is not installed. Install with: pip install pygbif")
                logger.warning("Falling back to WoRMS for taxonomy backfilling")
                source = 'worms'
            else:
                source = args.source
                
            # Perform the backfilling
            filtered_data = backfill_taxonomy(filtered_data, source=source)
            logger.info(f"{source.upper()} backfilling complete")
            
            # Save a log of all taxonomy changes made
            if BACKFILL_CHANGES:
                changes_df = pd.DataFrame(BACKFILL_CHANGES)
                changes_log = args.changes_log
                # Ensure the changes log path exists
                changes_dir = os.path.dirname(changes_log)
                if changes_dir and not os.path.exists(changes_dir):
                    os.makedirs(changes_dir)
                changes_df.to_csv(changes_log, sep='\t', index=False)
                logger.info(f"Saved {len(BACKFILL_CHANGES)} taxonomy changes to {changes_log}")
            else:
                logger.info("No taxonomy changes were made during backfilling")
        else:
            logger.info("Skipping taxonomy backfilling (--no-backfill flag used)")
        
        # 5. Forward fill missing taxonomy ranks
        logger.info("Forward filling missing taxonomy ranks...")
        final_data = forward_fill_taxonomy(filtered_data)
        logger.info("Forward filling complete")
        
        # 6. Separate metazoan/Rhodophyta from other taxa
        if args.separate_by_phylum:
            logger.info("Separating metazoan and non-metazoan taxa...")
            metazoan_file = args.output_file.replace('.tsv', f'_metazoan_{args.marker}.tsv')
            non_metazoan_file = args.output_file.replace('.tsv', f'_non_metazoan_{args.marker}.tsv')
            metazoan_data, non_metazoan_data = separate_taxa_by_phylum(final_data, METAZOAN_PLUS_RHODOPHYTA)

            # Save the separated data
            metazoan_data.to_csv(metazoan_file, sep='\t', index=False)
            non_metazoan_data.to_csv(non_metazoan_file, sep='\t', index=False)
            logger.info(f"Saving metazoan data to: {metazoan_file}")
            logger.info(f"Saving non-metazoan data to: {non_metazoan_file}")

        # Print phylum summary
        phylum_counts = pd.Series(dtype=int)
        if not final_data.empty and 'phylum' in final_data.columns:
            phylum_counts = final_data.groupby('phylum').size().sort_values(ascending=False)
        logger.info("Phylum distribution in full dataset:")
        for phylum, count in phylum_counts.items():
            logger.info(f"  {phylum}: {count} ASVs")

        # Save the processed data (original combined output)
        # Ensure the output path exists
        output_dir = os.path.dirname(args.output_file)
        if output_dir and not os.path.exists(output_dir):
            os.makedirs(output_dir)
        final_data.to_csv(args.output_file, sep='\t', index=False)
        logger.info(f"Saved processed data to {args.output_file}")

        # Save taxonomy cache
        if args.cache_file:
            save_taxonomy_cache(args.cache_file)

            return 0
        
    except Exception as e:
        logger.error(f"Error processing data: {e}")
        raise
        return 1

if __name__ == "__main__":
    sys.exit(main())