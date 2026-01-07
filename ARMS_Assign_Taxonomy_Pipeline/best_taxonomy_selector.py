#!/usr/bin/env python3
# Author: Modified from Nicole Pittoors' script
# Date: April 22, 2025
"""
Configurable Best Taxonomy Hit Selector

This script combines taxonomy information from different sources
(NMNH voucher, NCBI, and MetaZooGene) and selects the best hit for each ASV
based on database priority and percent identity.

Key features:
- Configurable for different marker types (COI or 18S)
- Support for full (NMNH+MZG+NCBI) or reduced (NMNH+MZG) database modes
- Configurable identity and alignment length thresholds
- Taxonomic filtering to exclude non-target organisms
- Option to prioritize hits by database source (MZG > NMNH > NCBI)

Input files:
- nmnh_taxonomy_results.tsv: NMNH voucher BLAST results with taxonomy
- ncbi_filtered_results.tsv: NCBI BLAST results (optional)
- MZG_taxonomy_results.tsv: MetaZooGene BLAST results with taxonomy

Output:
- combined_taxonomy_results.tsv: Best hit for each ASV with complete taxonomy
"""

import pandas as pd
import re
import argparse
from Bio import Entrez
import time
import sys
import logging
import os

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Define taxonomic groups to filter out from NCBI results
# Common for both COI and 18S
EXCLUDED_ORDERS = {
    # Birds
    'Accipitriformes', 'Anseriformes', 'Charadriiformes', 'Cuculiformes',
    'Dinornithiformes', 'Galliformes', 'Passeriformes', 'Pelecaniformes',
    'Piciformes', 'Procellariiformes', 'Psittaciformes', 'Sphenisciformes',
    'Trogoniformes', 'Apodiformes',
    
    # Mammals, reptiles, amphibians
    'Carnivora', 'Primates', 'Rodentia', 'Squamata', 'Anura',
    
    # Insects and mites
    'Coleoptera', 'Diptera', 'Hymenoptera', 'Lepidoptera',
    'Acariformes', 'Parasitiformes'
}

# Marker-specific settings
MARKER_SETTINGS = {
    'COI': {
        'min_identity': 85,
        'min_alignment': 200,
        'excluded_orders': EXCLUDED_ORDERS
    },
    '18S': {
        'min_identity': 85,
        'min_alignment': 300,
        'excluded_orders': EXCLUDED_ORDERS
    }
}

# Default database priority order (highest to lowest)
DATABASE_PRIORITY = {
    'MZG': 1,
    'NMNH': 2, 
    'NCBI': 3
}

def extract_accession(sseqid):
    """Extract the accession number from the NCBI BLAST sseqid field"""
    if not isinstance(sseqid, str):
        sseqid = str(sseqid)
        
    # Try to match the pattern in gi|number|gb|XXXXX| format
    match = re.search(r'gb\|([^|]+)', sseqid)
    if match:
        return match.group(1)
    
    # Try to match the pattern in gi|number|emb|XXXXX| format
    match = re.search(r'emb\|([^|]+)', sseqid)
    if match:
        return match.group(1)
    
    # Check for direct accession format (e.g., KX866738.1)
    match = re.search(r'[A-Z]{1,2}\d{5,}\.\d+', sseqid)
    if match:
        return match.group(0)
    
    # If no pattern matches, return the whole sseqid as fallback
    return sseqid

def fetch_taxonomy_from_ncbi(accession, excluded_orders):
    """Fetch taxonomy information from NCBI for a given accession number"""
    try:
        # Fetch data from NCBI
        handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="xml")
        record = Entrez.read(handle)
        handle.close()
        
        if not record or len(record) == 0:
            logger.warning(f"No record found for {accession}")
            return {}
        
        rec = record[0]
        organism = rec.get('GBSeq_organism', '')
        taxonomy_string = rec.get('GBSeq_taxonomy', '')
        
        # Parse taxonomy string
        taxonomy = {}
        if taxonomy_string:
            taxa = [t.strip() for t in taxonomy_string.split(';')]
            
            # Initialize taxonomy dict with empty values
            taxonomy = {
                'kingdom': '',
                'phylum': '',
                'class': '',
                'order': '',
                'family': '',
                'genus': '',
                'species': '',
                'scientific_name': ''
            }
            
            # Basic taxonomy processing - assign levels based on position
            if len(taxa) >= 1:
                taxonomy['kingdom'] = taxa[0]
            if len(taxa) >= 2:
                taxonomy['phylum'] = taxa[1]
            if len(taxa) >= 3:
                taxonomy['class'] = taxa[2]
            if len(taxa) >= 4:
                taxonomy['order'] = taxa[3]
            if len(taxa) >= 5:
                taxonomy['family'] = taxa[4]
        
        # Process organism name to extract genus and species
        if organism:
            parts = organism.split()
            if len(parts) >= 1:
                taxonomy['genus'] = parts[0]
            if len(parts) >= 2:
                # Check if the second part is actually a species name
                if not parts[1].startswith(('sp.', 'cf.', 'aff.')):
                    taxonomy['species'] = parts[1]
            
            # Set scientific name
            if taxonomy['genus'] and taxonomy['species']:
                taxonomy['scientific_name'] = f"{taxonomy['genus']} {taxonomy['species']}"
            elif taxonomy['genus']:
                taxonomy['scientific_name'] = taxonomy['genus']
        
        # Check if order is in excluded list
        order = taxonomy.get('order', '')
        if order in excluded_orders:
            taxonomy['excluded'] = True
        else:
            taxonomy['excluded'] = False
        
        logger.info(f"Fetched taxonomy for {accession}: {organism}")
        return taxonomy
    
    except Exception as e:
        logger.error(f"Error fetching taxonomy for {accession}: {str(e)}")
        return {}

def read_nmnh_data(file_path, min_identity, min_alignment):
    """Read and process NMNH voucher data with filtering"""
    try:
        if not os.path.exists(file_path):
            logger.warning(f"NMNH data file not found: {file_path}")
            return pd.DataFrame()
            
        nmnh_df = pd.read_csv(file_path, sep='\t')
        
        # Standardize column names
        if 'scientific_nammismatch' in nmnh_df.columns:  # Handle truncated header
            # Fix truncated column name and split into two
            nmnh_df['scientific_name'] = ''
            nmnh_df['mismatch'] = 0
            nmnh_df = nmnh_df.drop(columns=['scientific_nammismatch'])
        
        # Ensure all expected columns exist
        required_columns = ['qseqid', 'sseqid', 'pident', 'length', 'database', 
                           'phylum', 'class', 'order', 'family', 'genus', 'species']
        
        for col in required_columns:
            if col not in nmnh_df.columns:
                nmnh_df[col] = ''
        
        # Add scientific_name if not present
        if 'scientific_name' not in nmnh_df.columns:
            nmnh_df['scientific_name'] = nmnh_df.apply(
                lambda x: f"{x['genus']} {x['species']}" if x['genus'] and x['species'] else 
                         (x['genus'] if x['genus'] else ''),
                axis=1
            )
        
        # Standardize column names
        nmnh_df = nmnh_df.rename(columns={'qseqid': 'ASV_ID'})
        
        # Apply filters for minimum identity and alignment length
        if 'pident' in nmnh_df.columns:
            nmnh_df['pident'] = pd.to_numeric(nmnh_df['pident'], errors='coerce')
            nmnh_df = nmnh_df[nmnh_df['pident'] >= min_identity]
        
        if 'length' in nmnh_df.columns:
            nmnh_df['length'] = pd.to_numeric(nmnh_df['length'], errors='coerce')
            nmnh_df = nmnh_df[nmnh_df['length'] >= min_alignment]
        
        logger.info(f"Read {len(nmnh_df)} NMNH voucher records passing filters")
        return nmnh_df
    
    except Exception as e:
        logger.error(f"Error reading NMNH data: {str(e)}")
        return pd.DataFrame()

def read_ncbi_data(file_path, min_identity, min_alignment):
    """Read and process NCBI data with filtering"""
    try:
        if not os.path.exists(file_path):
            logger.warning(f"NCBI data file not found: {file_path}")
            return pd.DataFrame()
            
        # Define column names for NCBI BLAST results
        columns = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
                   'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'qcovs']
        
        # Read the BLAST results
        ncbi_df = pd.read_csv(file_path, sep='\t', header=None)
        
        # Assign appropriate column names based on the number of columns
        num_cols = ncbi_df.shape[1]
        if num_cols <= len(columns):
            ncbi_df.columns = columns[:num_cols]
        else:
            ncbi_df.columns = columns + [f'col_{i}' for i in range(len(columns), num_cols)]
        
        # Extract accession numbers
        ncbi_df['accession'] = ncbi_df['sseqid'].apply(extract_accession)
        
        # Apply filters for minimum identity and alignment length
        if 'pident' in ncbi_df.columns:
            ncbi_df['pident'] = pd.to_numeric(ncbi_df['pident'], errors='coerce')
            ncbi_df = ncbi_df[ncbi_df['pident'] >= min_identity]
        
        if 'length' in ncbi_df.columns:
            ncbi_df['length'] = pd.to_numeric(ncbi_df['length'], errors='coerce')
            ncbi_df = ncbi_df[ncbi_df['length'] >= min_alignment]
        
        # Select best hit per ASV based on percent identity
        if not ncbi_df.empty:
            best_hits = ncbi_df.loc[ncbi_df.groupby('qseqid')['pident'].idxmax()]
        else:
            best_hits = pd.DataFrame(columns=ncbi_df.columns)
        
        # Add database column
        best_hits['database'] = 'NCBI_nt'
        
        # Rename columns for consistency
        best_hits = best_hits.rename(columns={'qseqid': 'ASV_ID'})
        
        logger.info(f"Read {len(ncbi_df)} NCBI records, selected {len(best_hits)} best hits passing filters")
        return best_hits
    
    except Exception as e:
        logger.error(f"Error reading NCBI data: {str(e)}")
        return pd.DataFrame()

def read_mzg_data(file_path, min_identity, min_alignment):
    """Read and process MetaZooGene data with filtering"""
    try:
        if not os.path.exists(file_path):
            logger.warning(f"MZG data file not found: {file_path}")
            return pd.DataFrame()
            
        mzg_df = pd.read_csv(file_path, sep='\t')
        
        # Add database column if not present
        if 'database' not in mzg_df.columns:
            mzg_df['database'] = 'MetaZooGene'
        
        # Rename percent identity column if needed
        if 'Percent_Identity' in mzg_df.columns:
            mzg_df = mzg_df.rename(columns={'Percent_Identity': 'pident'})
        
        # Rename Length column if needed
        if 'Length' in mzg_df.columns:
            mzg_df = mzg_df.rename(columns={'Length': 'length'})
        
        # Rename taxonomy fields if needed (capitalize to lowercase)
        column_map = {
            'Kingdom': 'kingdom',
            'Phylum': 'phylum',
            'Class': 'class',
            'Order': 'order',
            'Family': 'family',
            'Genus': 'genus',
            'Species': 'species',
            'Scientific_Name': 'scientific_name'
        }
        mzg_df = mzg_df.rename(columns={k: v for k, v in column_map.items() if k in mzg_df.columns})
        
        # Apply filters for minimum identity and alignment length
        if 'pident' in mzg_df.columns:
            mzg_df['pident'] = pd.to_numeric(mzg_df['pident'], errors='coerce')
            mzg_df = mzg_df[mzg_df['pident'] >= min_identity]
        
        if 'length' in mzg_df.columns:
            mzg_df['length'] = pd.to_numeric(mzg_df['length'], errors='coerce')
            mzg_df = mzg_df[mzg_df['length'] >= min_alignment]
        
        logger.info(f"Read {len(mzg_df)} MetaZooGene records passing filters")
        return mzg_df
    
    except Exception as e:
        logger.error(f"Error reading MetaZooGene data: {str(e)}")
        return pd.DataFrame()

def select_best_hit(nmnh_df, ncbi_df, mzg_df, marker_type, db_mode, min_identity, excluded_orders, prioritize_db=True):
    """
    Select the best hit for each ASV based on database priority and percent identity
    
    Args:
        nmnh_df (DataFrame): NMNH data
        ncbi_df (DataFrame): NCBI data
        mzg_df (DataFrame): MZG data
        marker_type (str): Marker type (COI or 18S)
        db_mode (str): Database mode (full or reduced)
        min_identity (float): Minimum identity threshold
        excluded_orders (set): Set of taxonomic orders to exclude
        prioritize_db (bool): Whether to prioritize by database source over percent identity
        
    Returns:
        list: Selected best hits
    """
    # Get the unique list of ASVs across all datasets
    all_asvs = set()
    
    if not nmnh_df.empty:
        all_asvs.update(nmnh_df['ASV_ID'].unique())
    
    if not mzg_df.empty:
        all_asvs.update(mzg_df['ASV_ID'].unique())
    
    if db_mode == 'full' and not ncbi_df.empty:
        all_asvs.update(ncbi_df['ASV_ID'].unique())
    
    logger.info(f"Found {len(all_asvs)} unique ASVs across all datasets")
    
    # Process NCBI accessions to get taxonomy (only in full mode)
    ncbi_taxonomy = {}
    if db_mode == 'full' and not ncbi_df.empty:
        logger.info("Fetching taxonomy for NCBI accessions...")
        for i, row in ncbi_df.iterrows():
            accession = row['accession']
            if accession not in ncbi_taxonomy:
                logger.info(f"Processing {i+1}/{len(ncbi_df)}: {accession}")
                taxonomy = fetch_taxonomy_from_ncbi(accession, excluded_orders)
                ncbi_taxonomy[accession] = taxonomy
                time.sleep(0.5)  # Avoid overwhelming NCBI servers
    
    # Combine all hits into a single list for comparison
    all_hits = []
    
    # Track ASVs with no hits meeting criteria
    no_hits_asvs = set(all_asvs)
    
    # Process each ASV
    for asv_id in sorted(all_asvs):
        hits = []
        
        # Get NMNH hits
        if not nmnh_df.empty:
            nmnh_hits = nmnh_df[nmnh_df['ASV_ID'] == asv_id]
            for _, hit in nmnh_hits.iterrows():
                try:
                    pident = float(hit['pident'])
                    if pident >= min_identity:
                        hits.append({
                            'ASV_ID': asv_id,
                            'source': 'NMNH',
                            'pident': pident,
                            'data': hit,
                            'db_priority': DATABASE_PRIORITY.get('NMNH', 2)
                        })
                except (ValueError, TypeError):
                    logger.warning(f"Invalid percent identity value for NMNH hit on ASV {asv_id}")
        
        # Get NCBI hits (only in full mode)
        if db_mode == 'full' and not ncbi_df.empty:
            ncbi_hits = ncbi_df[ncbi_df['ASV_ID'] == asv_id]
            for _, hit in ncbi_hits.iterrows():
                try:
                    pident = float(hit['pident'])
                    if pident >= min_identity:
                        accession = hit['accession']
                        taxonomy = ncbi_taxonomy.get(accession, {})
                        
                        # Check if this hit belongs to an excluded order
                        excluded = taxonomy.get('excluded', False)
                        
                        hits.append({
                            'ASV_ID': asv_id,
                            'source': 'NCBI',
                            'pident': pident,
                            'excluded': excluded,
                            'data': hit,
                            'taxonomy': taxonomy,
                            'db_priority': DATABASE_PRIORITY.get('NCBI', 3)
                        })
                except (ValueError, TypeError):
                    logger.warning(f"Invalid percent identity value for NCBI hit on ASV {asv_id}")
        
        # Get MetaZooGene hits
        if not mzg_df.empty:
            mzg_hits = mzg_df[mzg_df['ASV_ID'] == asv_id]
            for _, hit in mzg_hits.iterrows():
                try:
                    pident = float(hit['pident'])
                    if pident >= min_identity:
                        hits.append({
                            'ASV_ID': asv_id,
                            'source': 'MZG',
                            'pident': pident,
                            'data': hit,
                            'db_priority': DATABASE_PRIORITY.get('MZG', 1)
                        })
                except (ValueError, TypeError):
                    logger.warning(f"Invalid percent identity value for MZG hit on ASV {asv_id}")
        
        # Sort hits by database priority (if enabled) and percent identity
        if prioritize_db:
            # First by database priority, then by percent identity
            hits.sort(key=lambda x: (x['db_priority'], -x['pident']))
            if hits:
                logger.debug(f"ASV {asv_id}: Prioritizing by database. Selected {hits[0]['source']} with {hits[0]['pident']}% identity")
        else:
            # Sort only by percent identity (descending)
            hits.sort(key=lambda x: x['pident'], reverse=True)
            if hits:
                logger.debug(f"ASV {asv_id}: Prioritizing by identity. Selected {hits[0]['source']} with {hits[0]['pident']}% identity")
        
        if hits:
            # Remove this ASV from the no_hits list
            no_hits_asvs.discard(asv_id)
            
            # If the best hit is from NCBI and is excluded, select the next best non-NCBI hit
            if db_mode == 'full' and hits[0]['source'] == 'NCBI' and hits[0].get('excluded', False):
                logger.info(f"ASV {asv_id}: Best hit is from excluded NCBI order, looking for alternative")
                
                # Look for the next best hit from NMNH or MZG
                alternative_hit = None
                for hit in hits[1:]:
                    if hit['source'] != 'NCBI':
                        alternative_hit = hit
                        break
                
                if alternative_hit:
                    logger.info(f"  Selected alternative hit from {alternative_hit['source']} with {alternative_hit['pident']}% identity")
                    all_hits.append(alternative_hit)
                else:
                    # If no alternative, use the NCBI hit despite being excluded
                    logger.info(f"  No alternative found, using excluded NCBI hit")
                    all_hits.append(hits[0])
            else:
                # Use the best hit
                all_hits.append(hits[0])
    
    # Log number of ASVs with no hits
    logger.info(f"{len(no_hits_asvs)} ASVs had no hits meeting criteria")
    
    # Log statistics about database sources in the selected hits
    source_counts = {}
    for hit in all_hits:
        source = hit['source']
        source_counts[source] = source_counts.get(source, 0) + 1
    
    logger.info("Selected hits by source:")
    for source, count in source_counts.items():
        logger.info(f"  {source}: {count} ASVs ({count/len(all_hits)*100:.1f}%)")
    
    return all_hits, no_hits_asvs

def create_combined_taxonomy(best_hits):
    """
    Create a standardized combined taxonomy dataframe from the best hits
    
    Args:
        best_hits (list): List of best hit dictionaries
        
    Returns:
        pandas.DataFrame: Combined taxonomy dataframe
    """
    # Create a list to hold standardized records
    records = []
    
    for hit in best_hits:
        source = hit['source']
        data = hit['data']
        
        # Start with a base record with empty values
        record = {
            'ASV_ID': hit['ASV_ID'],
            'database': source,
            'pident': hit['pident'],
            'sseqid': '',
            'length': '',
            'qcovs': '',
            'evalue': '',
            'bitscore': '',
            'specimen_id': '',
            'phylum': '',
            'class': '',
            'order': '',
            'family': '',
            'genus': '',
            'species': '',
            'scientific_name': '',
            'mismatch': '',
            'gapopen': '',
            'qstart': '',
            'qend': '',
            'sstart': '',
            'send': '',
            'accession': ''
        }
        
        # Fill in available data based on source
        if source == 'NCBI':
            taxonomy = hit.get('taxonomy', {})
            
            record.update({
                'sseqid': data['sseqid'],
                'length': data['length'],
                'qcovs': data.get('qcovs', ''),
                'evalue': data['evalue'],
                'bitscore': data['bitscore'],
                'mismatch': data.get('mismatch', ''),
                'gapopen': data.get('gapopen', ''),
                'qstart': data.get('qstart', ''),
                'qend': data.get('qend', ''),
                'sstart': data.get('sstart', ''),
                'send': data.get('send', ''),
                'accession': data['accession'],
                'phylum': taxonomy.get('phylum', ''),
                'class': taxonomy.get('class', ''),
                'order': taxonomy.get('order', ''),
                'family': taxonomy.get('family', ''),
                'genus': taxonomy.get('genus', ''),
                'species': taxonomy.get('species', ''),
                'scientific_name': taxonomy.get('scientific_name', '')
            })
            
        elif source == 'NMNH':
            record.update({
                'sseqid': data['sseqid'],
                'length': data['length'],
                'qcovs': data.get('qcovs', ''),
                'evalue': data.get('evalue', ''),
                'bitscore': data.get('bitscore', ''),
                'specimen_id': data.get('specimen_id', ''),
                'phylum': data.get('phylum', ''),
                'class': data.get('class', ''),
                'order': data.get('order', ''),
                'family': data.get('family', ''),
                'genus': data.get('genus', ''),
                'species': data.get('species', ''),
                'scientific_name': data.get('scientific_name', ''),
                'accession': data.get('specimen_id', '')
            })
            
        elif source == 'MZG':
            record.update({
                'sseqid': data.get('Reference_ID', ''),
                'length': data.get('length', ''),
                'phylum': data.get('phylum', ''),
                'class': data.get('class', ''),
                'order': data.get('order', ''),
                'family': data.get('family', ''),
                'genus': data.get('genus', ''),
                'species': data.get('species', ''),
                'scientific_name': data.get('scientific_name', '')
            })
        
        records.append(record)
    
    # Create dataframe
    df = pd.DataFrame(records)
    
    # Ensure columns have appropriate types
    numeric_cols = ['pident', 'length', 'qcovs', 'evalue', 'bitscore', 
                   'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send']
    
    for col in numeric_cols:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors='coerce')
    
    return df

def create_no_hits_file(no_hits_asvs, output_file):
    """Create a file with ASVs that had no hits meeting criteria"""
    no_hits_df = pd.DataFrame({'ASV_ID': list(no_hits_asvs)})
    no_hits_df.to_csv(output_file, sep='\t', index=False)
    logger.info(f"Saved {len(no_hits_asvs)} ASVs with no hits to {output_file}")

def main():
    parser = argparse.ArgumentParser(description="Select best taxonomic hit from multiple BLAST sources")
    
    # Basic options
    parser.add_argument("--output", default="combined_taxonomy_results.tsv",
                        help="Output file for combined taxonomy results (default: combined_taxonomy_results.tsv)")
    parser.add_argument("--email", default=None,
                        help="Your email for NCBI Entrez (required for NCBI API access)")
    parser.add_argument("--marker", choices=['COI', '18S'], default='COI',
                        help="Marker type: COI or 18S (default: COI)")
    parser.add_argument("--db-mode", choices=['full', 'reduced'], default='full',
                        help="Database mode: full (NMNH+MZG+NCBI) or reduced (NMNH+MZG) (default: full)")
    
    # Input files
    parser.add_argument("--nmnh", default="nmnh_taxonomy_results.tsv", 
                        help="NMNH voucher taxonomy results file (default: nmnh_taxonomy_results.tsv)")
    parser.add_argument("--ncbi", default="ncbi_filtered_results.tsv",
                        help="NCBI filtered BLAST results file (default: ncbi_filtered_results.tsv)")
    parser.add_argument("--mzg", default="MZG_taxonomy_results.tsv",
                        help="MetaZooGene taxonomy results file (default: MZG_taxonomy_results.tsv)")
    
    # Filtering parameters
    parser.add_argument("--min-identity", type=float, default=None,
                        help="Minimum percent identity threshold (default: depends on marker)")
    parser.add_argument("--min-alignment", type=int, default=None,
                        help="Minimum alignment length (default: depends on marker)")
    
    # New option for database prioritization
    parser.add_argument("--prioritize-db", action="store_true", default=True,
                        help="Prioritize by database source (MZG > NMNH > NCBI) over percent identity (default: enabled)")
    parser.add_argument("--no-prioritize-db", action="store_false", dest="prioritize_db",
                        help="Sort hits only by percent identity, not by database source")
    
    # Database priority order customization
    parser.add_argument("--db-order", default="MZG,NMNH,NCBI",
                        help="Custom database priority order (comma-separated, highest priority first, default: MZG,NMNH,NCBI)")
    
    # Output for no-hits ASVs
    parser.add_argument("--no-hits-file", default="no_hits.tsv",
                        help="Output file for ASVs with no hits meeting criteria (default: no_hits.tsv)")
    
    args = parser.parse_args()
    
    # Set marker-specific settings
    marker_settings = MARKER_SETTINGS.get(args.marker, MARKER_SETTINGS['COI'])
    min_identity = args.min_identity if args.min_identity is not None else marker_settings['min_identity']
    min_alignment = args.min_alignment if args.min_alignment is not None else marker_settings['min_alignment']
    excluded_orders = marker_settings['excluded_orders']
    
    # Set up custom database priority order if provided
    if args.db_order:
        try:
            db_priorities = args.db_order.split(',')
            # Reset the global database priority dictionary
            global DATABASE_PRIORITY
            DATABASE_PRIORITY = {db: i+1 for i, db in enumerate(db_priorities)}
            logger.info(f"Using custom database priority order: {args.db_order}")
        except Exception as e:
            logger.error(f"Error parsing database order: {e}")
            logger.info("Using default database priority: MZG > NMNH > NCBI")
    
    # Set email for NCBI
    if args.email:
        Entrez.email = args.email
    elif args.db_mode == 'full':
        logger.warning("Email not provided for NCBI Entrez. This may limit NCBI API access.")
    
    # Log configuration
    logger.info(f"Running with marker: {args.marker}, database mode: {args.db_mode}")
    logger.info(f"Minimum identity threshold: {min_identity}%, minimum alignment length: {min_alignment}bp")
    logger.info(f"Database prioritization: {'Enabled' if args.prioritize_db else 'Disabled'}")
    
    # Read input files
    logger.info("Reading input files...")
    
    nmnh_df = read_nmnh_data(args.nmnh, min_identity, min_alignment)
    mzg_df = read_mzg_data(args.mzg, min_identity, min_alignment)
    
    # Only read NCBI data in full mode
    if args.db_mode == 'full':
        ncbi_df = read_ncbi_data(args.ncbi, min_identity, min_alignment)
    else:
        ncbi_df = pd.DataFrame()  # Empty dataframe in reduced mode
    
    # Check if essential files were read successfully
    if nmnh_df.empty and mzg_df.empty and (args.db_mode == 'reduced' or ncbi_df.empty):
        logger.error("Failed to read any input files. Exiting.")
        sys.exit(1)
    
    # Select best hits
    logger.info("Selecting best hits...")
    best_hits, no_hits_asvs = select_best_hit(
        nmnh_df, ncbi_df, mzg_df, 
        args.marker, args.db_mode, min_identity, excluded_orders,
        prioritize_db=args.prioritize_db
    )
    
    # Create combined taxonomy
    logger.info("Creating combined taxonomy...")
    combined_df = create_combined_taxonomy(best_hits)
    
    # Save results
    combined_df.to_csv(args.output, sep='\t', index=False)
    logger.info(f"Saved combined taxonomy results for {len(combined_df)} ASVs to {args.output}")
    
    # Save ASVs with no hits
    create_no_hits_file(no_hits_asvs, args.no_hits_file)
    
    # Print summary
    sources = combined_df['database'].value_counts()
    logger.info("Summary of data sources in final results:")
    for source, count in sources.items():
        logger.info(f"  {source}: {count} ASVs")
    
    logger.info(f"Total ASVs with hits: {len(combined_df)}")
    logger.info(f"Total ASVs without hits: {len(no_hits_asvs)}")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())