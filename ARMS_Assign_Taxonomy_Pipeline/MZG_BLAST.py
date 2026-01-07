#!/usr/bin/env python3
"""
Configurable MetaZooGene BLAST script for performing a local blast on downloaded MZG sequences
Author: Modified from Nicole Pittoors' original script
Date: April 15, 2025

This script:
1. Takes ASV sequences (fasta) and the MetaZooGene reference database as input
2. Creates a BLAST database if needed
3. Runs BLAST to find matches for each ASV with configurable parameters
4. Correctly parses taxonomy from MZG headers, handling the genus_species format
5. Selects the best hit for each ASV based on percent identity
6. Includes all ASVs in the output file, even those without BLAST hits
7. Outputs a taxonomy.tsv file with full taxonomy information

Key features:
- Supports both COI and 18S markers with appropriate settings
- Configurable identity threshold and alignment length
- Flexible taxonomy parsing for different marker types
"""

import os
import subprocess
import pandas as pd
from Bio import SeqIO
import argparse
import logging

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

def run_blast(query_fasta, reference_db, output_file, min_identity=0, min_alignment=0, threads=8):
    """
    Run blastn on input query against the reference database with configurable thresholds
    
    Args:
        query_fasta (str): Path to query FASTA file
        reference_db (str): Path to reference database FASTA file
        output_file (str): Path to output BLAST results
        min_identity (float): Minimum percent identity threshold (0 = no threshold)
        min_alignment (int): Minimum alignment length (0 = no threshold)
        threads (int): Number of threads for BLAST
        
    Returns:
        str: Path to BLAST results file
    """
    # Create BLAST database if .nin files don't exist
    if not os.path.exists(f"{reference_db}.nin"):
        logger.info(f"Creating BLAST database from {reference_db}...")
        subprocess.run([
            "makeblastdb", 
            "-in", reference_db, 
            "-dbtype", "nucl"
        ], check=True)
    
    # Build BLAST command
    blast_cmd = [
        "blastn",
        "-query", query_fasta,
        "-db", reference_db,
        "-out", output_file,
        "-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore",
        "-num_threads", str(threads),
        "-max_target_seqs", "5"  # Return 5 hits per query as recommended
    ]
    
    # Add identity threshold if specified
    if min_identity > 0:
        blast_cmd.extend(["-perc_identity", str(min_identity)])
    
    # Run BLAST
    logger.info(f"Running BLAST search with parameters: min_identity={min_identity}, min_alignment={min_alignment}")
    subprocess.run(blast_cmd, check=True)
    
    # Filter for alignment length if specified
    if min_alignment > 0:
        logger.info(f"Filtering results for minimum alignment length of {min_alignment}bp")
        # Create a temporary file
        temp_file = f"{output_file}.temp"
        
        with open(output_file, 'r') as infile, open(temp_file, 'w') as outfile:
            for line in infile:
                fields = line.strip().split('\t')
                if len(fields) >= 4 and int(fields[3]) >= min_alignment:
                    outfile.write(line)
        
        # Replace the original file
        os.replace(temp_file, output_file)
    
    logger.info(f"BLAST results saved to {output_file}")
    return output_file

def parse_taxonomy_from_header(header, marker_type="COI"):
    """
    Parse taxonomy from FASTA header in MetaZooGene format based on marker type
    
    Args:
        header (str): FASTA header string
        marker_type (str): Marker type (COI or 18S)
        
    Returns:
        dict: Taxonomy information
    """
    # Split by semicolons to get taxonomy levels
    taxonomy_parts = header.split(';')
    
    # Initialize empty taxonomy
    taxonomy = {
        "Kingdom": "", "Phylum": "", "Class": "", "Order": "", 
        "Family": "", "Genus": "", "Species": "", "Scientific_Name": ""
    }
    
    # Handle different marker types
    if marker_type == "COI":
        # Standard MZG COI format
        if len(taxonomy_parts) >= 7:
            taxonomy["Kingdom"] = taxonomy_parts[0]
            taxonomy["Phylum"] = taxonomy_parts[1]
            taxonomy["Class"] = taxonomy_parts[2]
            taxonomy["Order"] = taxonomy_parts[3]
            taxonomy["Family"] = taxonomy_parts[4]
            taxonomy["Genus"] = taxonomy_parts[5]
            
            # For species, split at the underscore to get genus and species separately
            genus_species = taxonomy_parts[6]
            try:
                genus, species = genus_species.split('_', 1)  # Split at first underscore only
                taxonomy["Species"] = species
                taxonomy["Scientific_Name"] = f"{genus} {species}"
            except ValueError:
                # Handle case where there's no underscore
                taxonomy["Species"] = ""
                taxonomy["Scientific_Name"] = genus_species
    
    elif marker_type == "18S":
        # Handle 18S format which might be different
        # Adjust these indices based on your 18S reference format
        if len(taxonomy_parts) >= 7:
            taxonomy["Kingdom"] = taxonomy_parts[0]
            taxonomy["Phylum"] = taxonomy_parts[1]
            taxonomy["Class"] = taxonomy_parts[2]
            taxonomy["Order"] = taxonomy_parts[3]
            taxonomy["Family"] = taxonomy_parts[4]
            taxonomy["Genus"] = taxonomy_parts[5]
            
            # For 18S, the format might be different - adjust as needed
            if '_' in taxonomy_parts[6]:
                genus, species = taxonomy_parts[6].split('_', 1)
                taxonomy["Species"] = species
                taxonomy["Scientific_Name"] = f"{genus} {species}"
            else:
                taxonomy["Species"] = ""
                taxonomy["Scientific_Name"] = taxonomy_parts[6]
    
    return taxonomy

def assign_taxonomy(blast_results, reference_fasta, query_fasta, marker_type="COI", min_identity=0, min_alignment=0):
    """
    Assign taxonomy based on BLAST results
    Also ensures every ASV in the query file gets an entry, even if no BLAST hit was found
    
    Args:
        blast_results (str): Path to BLAST results file
        reference_fasta (str): Path to reference FASTA file
        query_fasta (str): Path to query FASTA file
        marker_type (str): Marker type (COI or 18S)
        min_identity (float): Minimum percent identity threshold for filtering
        min_alignment (int): Minimum alignment length threshold for filtering
        
    Returns:
        pandas.DataFrame: Taxonomy assignments
    """
    # Read all query sequences to get a complete list of ASVs
    query_ids = [record.id for record in SeqIO.parse(query_fasta, "fasta")]
    
    # Read BLAST results
    try:
        blast_df = pd.read_csv(blast_results, sep='\t', header=None, 
                            names=["qseqid", "sseqid", "pident", "length", "mismatch", 
                                    "gapopen", "qstart", "qend", "sstart", "send", 
                                    "evalue", "bitscore"])
    except pd.errors.EmptyDataError:
        logger.warning("Warning: BLAST results file is empty. No hits found.")
        blast_df = pd.DataFrame(columns=["qseqid", "sseqid", "pident", "length", "mismatch", 
                                   "gapopen", "qstart", "qend", "sstart", "send", 
                                   "evalue", "bitscore"])
    
    # Load reference sequences headers for taxonomy lookup
    ref_headers = {}
    for record in SeqIO.parse(reference_fasta, "fasta"):
        ref_headers[record.id] = record.description
    
    # Create taxonomy dataframe with all BLAST hits
    all_hits = []
    
    # Track which ASVs have hits
    asvs_with_hits = set()
    
    for _, row in blast_df.iterrows():
        query_id = row["qseqid"]
        subject_id = row["sseqid"]
        percent_identity = row["pident"]
        alignment_length = row["length"]
        
        # Skip hits below identity or alignment threshold (if specified)
        if min_identity > 0 and percent_identity < min_identity:
            continue
        if min_alignment > 0 and alignment_length < min_alignment:
            continue
            
        asvs_with_hits.add(query_id)
        
        # Parse taxonomy from the reference header
        if subject_id in ref_headers:
            taxonomy = parse_taxonomy_from_header(subject_id, marker_type)
            
            all_hits.append({
                "ASV_ID": query_id,
                "Reference_ID": subject_id,
                "Percent_Identity": percent_identity,
                "Length": row["length"],
                "Kingdom": taxonomy["Kingdom"],
                "Phylum": taxonomy["Phylum"],
                "Class": taxonomy["Class"],
                "Order": taxonomy["Order"],
                "Family": taxonomy["Family"],
                "Genus": taxonomy["Genus"],
                "Species": taxonomy["Species"],
                "Scientific_Name": taxonomy["Scientific_Name"],
                "Marker": marker_type
            })
    
    # Convert to dataframe
    all_hits_df = pd.DataFrame(all_hits)
    
    # Select the best hit for each ASV based on percent identity
    taxonomy_results = []
    
    if not all_hits_df.empty:
        # Group by ASV_ID and select the hit with the highest percent identity
        best_hits = all_hits_df.loc[all_hits_df.groupby('ASV_ID')['Percent_Identity'].idxmax()]
        taxonomy_results.extend(best_hits.to_dict('records'))
    
    # Add entries for ASVs with no hits
    no_hits_asvs = set()
    for query_id in query_ids:
        if query_id not in asvs_with_hits:
            # Add a row with empty taxonomy for ASVs with no hits
            taxonomy_results.append({
                "ASV_ID": query_id,
                "Reference_ID": "No_hit",
                "Percent_Identity": 0.0,
                "Length": 0,
                "Kingdom": "",
                "Phylum": "",
                "Class": "",
                "Order": "",
                "Family": "",
                "Genus": "",
                "Species": "",
                "Scientific_Name": "",
                "Marker": marker_type
            })
            no_hits_asvs.add(query_id)
    
    # Convert to dataframe
    taxonomy_df = pd.DataFrame(taxonomy_results)
    
    # Log stats
    logger.info(f"Total ASVs: {len(query_ids)}")
    logger.info(f"ASVs with hits: {len(asvs_with_hits)}")
    logger.info(f"ASVs with no hits: {len(no_hits_asvs)}")
    
    return taxonomy_df, no_hits_asvs

def main():
    parser = argparse.ArgumentParser(description="Assign taxonomy to ASVs using BLAST against MetaZooGene database")
    parser.add_argument("--query", required=True, help="Query FASTA file containing ASVs")
    parser.add_argument("--reference", required=True, help="MetaZooGene reference FASTA file")
    parser.add_argument("--output", required=True, help="Output prefix for results")
    parser.add_argument("--threads", type=int, default=8, help="Number of threads for BLAST")
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
    
    logger.info(f"Running with marker: {args.marker}")
    logger.info(f"Minimum identity threshold: {min_identity}%, minimum alignment length: {min_alignment}bp")
    
    # Run BLAST with no minimum identity threshold initially (filtering will be done later)
    blast_out = f"{args.output}.blast.out"
    run_blast(args.query, args.reference, blast_out, min_identity=0, min_alignment=0, threads=args.threads)
    
    # Assign taxonomy - passing query file to ensure all ASVs are included
    taxonomy_df, no_hits_asvs = assign_taxonomy(
        blast_out, 
        args.reference, 
        args.query, 
        marker_type=args.marker,
        min_identity=min_identity,
        min_alignment=min_alignment
    )
    
    # Save taxonomy results
    taxonomy_df.to_csv(f"{args.output}.tsv", sep='\t', index=False)
    logger.info(f"Taxonomy assignments saved to {args.output}.tsv")
    
    # Save ASVs with no hits if requested
    if args.no_hits_file and no_hits_asvs:
        no_hits_df = pd.DataFrame({'ASV_ID': list(no_hits_asvs)})
        no_hits_df.to_csv(args.no_hits_file, sep='\t', index=False)
        logger.info(f"Saved {len(no_hits_asvs)} ASVs with no hits to {args.no_hits_file}")
    
    return 0

if __name__ == "__main__":
    main()