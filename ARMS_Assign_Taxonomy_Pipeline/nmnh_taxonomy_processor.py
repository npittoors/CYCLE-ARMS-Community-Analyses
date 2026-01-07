#!/usr/bin/env python3

"""
Configurable NMNH Voucher Taxonomy Processing Script

This script processes BLAST results against the NMNH voucher database and
properly formats taxonomic information for downstream analysis.

Key features:
- Support for multiple marker types (COI, 18S)
- Configurable thresholds for percent identity and alignment length
- Proper taxonomy parsing for different reference formats
- Output of no-hit sequences to a separate file

Input:
- BLAST results against NMNH voucher database

Output:
- Formatted taxonomy information for each ASV's best hit
- Optional file with ASVs that had no hits meeting criteria

The script:
1. Processes BLAST results from NMNH voucher database
2. Parses taxonomic information from the sseqid field
3. Correctly assigns taxonomic ranks based on a classification lookup
4. Handles Genus_species format appropriately
5. Formats output for compatibility with other datasets
6. Filters results based on marker-specific thresholds
"""

import pandas as pd
import re
import argparse
from collections import defaultdict
import os
import logging
import sys

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

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

# Classification lookup tables
# These help determine the appropriate taxonomic rank for each classification

PHYLUM_LEVEL = {
    'Arthropoda', 'Mollusca', 'Chordata', 'Bryozoa', 'Cnidaria', 
    'Porifera', 'Echinodermata', 'Annelida', 'Rhodophyta', 'Chlorophyta',
    'Sipuncula', 'Entoprocta'
}

CLASS_LEVEL = {
    'Demospongiae', 'Hydrozoa', 'Ophiuroidea', 'Gastropoda', 'Polyplacophora',
    'Ascidiacea', 'Polychaeta', 'Phlebobranchia', 'Aplousobranchia',
    'Heteroscleromorpha', 'Coronatae', 'Reptantia', 'Brachyura'
}

ORDER_LEVEL = {
    'Suberitida', 'Poecilosclerida', 'Chondrillida', 'Haplosclerida', 
    'Axinellida', 'Chondrillida', 'Dictyoceratida', 'Dendroceratida',
    'Leptothecata', 'Amphipoda', 'Tanaidacea', 'Isopoda', 'Caridea',
    'Cheilostomatida', 'Tubuliporina', 'Comatulida', 'Corallinales',
    'Polycladida', 'Myzostomida'
}

FAMILY_LEVEL = {
    'Scopalinidae', 'Polynoidae', 'Nereididae', 'Glyceridae', 'Facelinidae',
    'Lafoeidae', 'Amphinomidae', 'Paguridae', 'Calloporidae', 'Smittinidae',
    'Syllidae', 'Hesionidae', 'Eunicidae', 'Cerithiopsidae', 'Terebellidae',
    'Haleciidae', 'Pleurobranchidae', 'Bugulidae', 'Muricidae', 'Fissurellidae',
    'Lumbrineridae', 'Spionidae', 'Ischyroceridae', 'Serpulidae', 'Caprellidae',
    'Styelidae', 'Didemnidae', 'Pyuridae', 'Sabellidae', 'Schizoporellidae',
    'Arcidae', 'Limidae', 'Tubuliporidae', 'Clavelinidae', 'Isodictyidae',
    'Polymastiidae', 'Beaniidae', 'Capitellidae', 'Polyclinidae', 'Polycitoridae',
    'Palicidae', 'Palaemonidae', 'Comatulidae', 'Crisiidae', 'Podospongiidae',
    'Xanthidae', 'Zygophylacidae', 'Dictyodendrillidae', 'Ampeliscidae',
    'Linnaeoxanthidae', 'Discodorididae', 'Bugulididae', 'Gobiidae',
    'Peyssonneliaceae'
}

# Regular expression patterns
GENUS_SPECIES_PATTERN = re.compile(r'^([A-Z][a-z]+)_([a-z]+)$')
CF_PATTERN = re.compile(r'^cf\._(.+)$')

def determine_taxonomic_rank(taxon, marker_type='COI'):
    """
    Determine the appropriate taxonomic rank for a given taxon
    
    Args:
        taxon (str): The taxonomic name to evaluate
        marker_type (str): The marker type (COI or 18S)
        
    Returns:
        tuple: (rank, genus, species) where rank is the taxonomic level
    """
    # Check for Genus_species pattern
    genus_species_match = GENUS_SPECIES_PATTERN.match(taxon)
    if genus_species_match:
        return 'genus_species', genus_species_match.group(1), genus_species_match.group(2)
    
    # Check for cf. pattern (e.g., cf._Costoanachis_translirata)
    cf_match = CF_PATTERN.match(taxon)
    if cf_match:
        inner_taxon = cf_match.group(1)
        inner_match = GENUS_SPECIES_PATTERN.match(inner_taxon)
        if inner_match:
            return 'genus_species', inner_match.group(1), inner_match.group(2)
    
    # Check against known taxonomic levels
    if taxon in PHYLUM_LEVEL:
        return 'phylum', taxon, None
    elif taxon in CLASS_LEVEL:
        return 'class', taxon, None
    elif taxon in ORDER_LEVEL:
        return 'order', taxon, None
    elif taxon in FAMILY_LEVEL:
        return 'family', taxon, None
    else:
        # Default to genus if capitalized
        if taxon and taxon[0].isupper():
            return 'genus', taxon, None
        else:
            return 'unknown', taxon, None

def parse_nmnh_sseqid(sseqid, marker_type='COI'):
    """
    Parse the NMNH sseqid format based on marker type
    
    Args:
        sseqid (str): The subject sequence ID from BLAST results
        marker_type (str): The marker type (COI or 18S)
        
    Returns:
        dict: Dictionary with extracted taxonomic information
    """
    parts = sseqid.split('|')
    
    taxonomy = {
        'specimen_id': parts[0] if len(parts) > 0 else '',
        'phylum': '',
        'class': '',
        'order': '',
        'family': '',
        'genus': '',
        'species': '',
        'scientific_name': ''
    }
    
    # Extract phylum if available (second part)
    if len(parts) > 1:
        phylum = parts[1]
        taxonomy['phylum'] = phylum
    
    # Process the third part which could be various taxonomic ranks
    if len(parts) > 2:
        taxon = parts[2]
        rank, value, species = determine_taxonomic_rank(taxon, marker_type)
        
        if rank == 'genus_species':
            taxonomy['genus'] = value
            taxonomy['species'] = species
            taxonomy['scientific_name'] = f"{value} {species}"
        elif rank == 'phylum':
            taxonomy['phylum'] = value
        elif rank == 'class':
            taxonomy['class'] = value
        elif rank == 'order':
            taxonomy['order'] = value
        elif rank == 'family':
            taxonomy['family'] = value
        elif rank == 'genus':
            taxonomy['genus'] = value
            taxonomy['scientific_name'] = value
    
    # Handle the case where scientific_name is missing but genus is present
    if not taxonomy['scientific_name'] and taxonomy['genus']:
        if taxonomy['species']:
            taxonomy['scientific_name'] = f"{taxonomy['genus']} {taxonomy['species']}"
        else:
            taxonomy['scientific_name'] = taxonomy['genus']
    
    return taxonomy

def process_nmnh_blast(blast_file, output_file, marker_type='COI', min_identity=85, min_alignment=200, no_hits_file=None):
    """
    Process NMNH BLAST results and create a formatted taxonomy file
    
    Args:
        blast_file (str): Path to BLAST results file
        output_file (str): Path to output file
        marker_type (str): Marker type (COI or 18S)
        min_identity (float): Minimum percent identity threshold
        min_alignment (int): Minimum alignment length threshold
        no_hits_file (str): Optional path to file for ASVs with no hits
    """
    logger.info(f"Reading NMNH BLAST results from {blast_file}...")
    
    try:
        # Check if the file exists
        if not os.path.exists(blast_file):
            logger.error(f"BLAST file not found: {blast_file}")
            return
        
        # Define column names based on the format of NMNH BLAST results
        column_names = ['qseqid', 'sseqid', 'pident', 'length', 'qcovs', 'evalue', 'bitscore']
        
        # Read the file
        blast_df = pd.read_csv(blast_file, sep='\t', header=None)
        
        # Assign column names based on number of columns
        if blast_df.shape[1] <= len(column_names):
            blast_df.columns = column_names[:blast_df.shape[1]]
        else:
            blast_df.columns = column_names + [f'col_{i}' for i in range(len(column_names), blast_df.shape[1])]
        
        logger.info(f"Found {len(blast_df)} BLAST hits")
        
        # Apply filtering criteria
        filtered_df = blast_df.copy()
        
        # Convert columns to numeric to ensure proper filtering
        for col in ['pident', 'length']:
            if col in filtered_df.columns:
                filtered_df[col] = pd.to_numeric(filtered_df[col], errors='coerce')
        
        # Apply filters
        original_count = len(filtered_df)
        filtered_df = filtered_df[filtered_df['pident'] >= min_identity]
        filtered_df = filtered_df[filtered_df['length'] >= min_alignment]
        
        logger.info(f"Applied filters: ≥{min_identity}% identity, ≥{min_alignment}bp alignment")
        logger.info(f"Filtered from {original_count} to {len(filtered_df)} hits")
        
        # Get all unique ASV IDs from the original results to track no-hits
        all_asv_ids = set(blast_df['qseqid'].unique())
        hit_asv_ids = set(filtered_df['qseqid'].unique())
        no_hit_asv_ids = all_asv_ids - hit_asv_ids
        
        # Parse taxonomy information from sseqid
        result_data = []
        
        # Track best hit per ASV
        best_hits = defaultdict(lambda: {'pident': 0, 'hit': None})
        
        # Process each BLAST hit
        for i, row in filtered_df.iterrows():
            asv_id = row['qseqid']
            pident = float(row['pident'])
            
            # Track the best hit for each ASV
            if pident > best_hits[asv_id]['pident']:
                best_hits[asv_id]['pident'] = pident
                best_hits[asv_id]['hit'] = row
        
        # Process the best hits
        for asv_id, data in best_hits.items():
            row = data['hit']
            
            # Parse taxonomy from sseqid
            taxonomy = parse_nmnh_sseqid(row['sseqid'], marker_type)
            
            # Create the result row with all required columns
            result_row = {
                'qseqid': row['qseqid'],
                'sseqid': row['sseqid'],
                'pident': row['pident'],
                'length': row['length'],
                'qcovs': row.get('qcovs', 0),
                'evalue': row['evalue'],
                'bitscore': row['bitscore'],
                'database': 'NMNH_voucher',
                'specimen_id': taxonomy['specimen_id'],
                'phylum': taxonomy['phylum'],
                'class': taxonomy['class'],
                'order': taxonomy['order'],
                'family': taxonomy['family'],
                'genus': taxonomy['genus'],
                'species': taxonomy['species'],
                'scientific_name': taxonomy['scientific_name'],
                'mismatch': '',  # Not available in NMNH data
                'gapopen': '',   # Not available in NMNH data
                'qstart': '',    # Not available in NMNH data
                'qend': '',      # Not available in NMNH data
                'sstart': '',    # Not available in NMNH data
                'send': '',      # Not available in NMNH data
                'accession': taxonomy['specimen_id'],  # Use specimen_id as accession
                'marker': marker_type
            }
            
            result_data.append(result_row)
        
        # Create result dataframe and save to file
        result_df = pd.DataFrame(result_data)
        result_df.to_csv(output_file, sep='\t', index=False)
        
        logger.info(f"Processed {len(result_df)} ASVs with NMNH taxonomy information")
        logger.info(f"Saved results to {output_file}")
        
        # Save no-hit ASVs if requested
        if no_hits_file and no_hit_asv_ids:
            no_hits_df = pd.DataFrame({'ASV_ID': list(no_hit_asv_ids)})
            no_hits_df.to_csv(no_hits_file, sep='\t', index=False)
            logger.info(f"Saved {len(no_hit_asv_ids)} ASVs with no hits to {no_hits_file}")
        
    except Exception as e:
        logger.error(f"Error processing NMNH BLAST file: {str(e)}")
        raise

def main():
    parser = argparse.ArgumentParser(description="Process NMNH voucher BLAST results and extract taxonomy")
    parser.add_argument("--blast", default="NMNH_blast_arms_results.txt", 
                        help="Input NMNH BLAST results file (default: NMNH_blast_arms_results.txt)")
    parser.add_argument("--output", default="nmnh_taxonomy_results.tsv",
                        help="Output file for taxonomy results (default: nmnh_taxonomy_results.tsv)")
    parser.add_argument("--marker", choices=['COI', '18S'], default='COI',
                        help="Marker type: COI or 18S (default: COI)")
    parser.add_argument("--min-identity", type=float, default=None,
                        help="Minimum percent identity threshold (default: depends on marker)")
    parser.add_argument("--min-alignment", type=int, default=None,
                        help="Minimum alignment length (default: depends on marker)")
    parser.add_argument("--no-hits-file", default=None,
                        help="Output file for ASVs with no hits meeting criteria (optional)")
    
    args = parser.parse_args()
    
    # Set marker-specific settings if not provided
    marker_settings = MARKER_SETTINGS.get(args.marker, MARKER_SETTINGS['COI'])
    min_identity = args.min_identity if args.min_identity is not None else marker_settings['min_identity']
    min_alignment = args.min_alignment if args.min_alignment is not None else marker_settings['min_alignment']
    
    logger.info(f"Processing NMNH data for {args.marker} marker")
    logger.info(f"Using thresholds: min_identity={min_identity}%, min_alignment={min_alignment}bp")
    
    process_nmnh_blast(
        args.blast, 
        args.output, 
        marker_type=args.marker,
        min_identity=min_identity,
        min_alignment=min_alignment,
        no_hits_file=args.no_hits_file
    )
    
    logger.info("NMNH processing complete")

if __name__ == "__main__":
    main()