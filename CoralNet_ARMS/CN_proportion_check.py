#!/usr/bin/env python3
"""
Quick validation script to check proportion calculations
Run this after your main processing script
"""

import pandas as pd
import numpy as np
import os

def validate_proportion_files(counts_file, props_file, exclude_algae=False, description=""):
    """Validate proportion calculations from saved files."""
    
    print("="*50)
    print(f"VALIDATING: {props_file}")
    if description:
        print(f"Description: {description}")
    print("="*50)
    
    # Load files
    counts = pd.read_csv(counts_file, index_col=0)
    props = pd.read_csv(props_file, index_col=0)
    
    # Get numeric columns (skip sitelocality)
    numeric_cols = counts.columns[1:]
    
    # Calculate proportion sums
    proportion_sums = props[numeric_cols].sum(axis=1)
    
    # Check if close to 1.0
    close_to_one = np.abs(proportion_sums - 1.0) < 0.001
    
    print(f"Proportion sum statistics:")
    print(f"  Mean: {proportion_sums.mean():.6f}")
    print(f"  Min: {proportion_sums.min():.6f}")
    print(f"  Max: {proportion_sums.max():.6f}")
    print(f"  Samples summing to ~1.0: {close_to_one.sum()}/{len(proportion_sums)}")
    
    # Show problematic samples
    problems = proportion_sums[~close_to_one]
    if len(problems) > 0:
        print(f"\nProblematic samples:")
        for sample, sum_val in problems.head(5).items():
            print(f"  {sample}: {sum_val:.6f}")
    
    # Manual verification of a few samples
    print(f"\nManual verification (3 random samples):")
    sample_indices = np.random.choice(counts.index, min(3, len(counts)), replace=False)
    
    for sample_id in sample_indices:
        sample_counts = counts.loc[sample_id, numeric_cols]
        sample_props = props.loc[sample_id, numeric_cols]
        
        total_count = sample_counts.sum()
        total_prop = sample_props.sum()
        
        print(f"\n  Sample {sample_id}:")
        print(f"    Total count: {total_count}")
        print(f"    Total proportion: {total_prop:.6f}")
        
        # Check specific categories if present
        categories_to_check = ['Unavailable', 'Chlorophyta', 'Ochrophyta', 'No_recruitment', 'Rhodophyta']
        for category in categories_to_check:
            if category in sample_counts:
                cat_count = sample_counts[category]
                cat_prop = sample_props[category]
                print(f"    {category}: {cat_count} count, {cat_prop:.4f} prop")
    
    return close_to_one.sum() == len(proportion_sums)

def main():
    """Run validation on all available files."""
    
    # Define file pairs to check
    file_pairs = [
        {
            'counts': "CoralNet_phylum_counts.csv",
            'props': "CoralNet_phylum_proportions.csv", 
            'description': "Standard files (all categories included)",
            'exclude_algae': False
        },
        {
            'counts': "CoralNet_phylum_counts_no_algae.csv",
            'props': "CoralNet_phylum_proportions_no_algae.csv",
            'description': "Excluding Chlorophyta, Ochrophyta, and Unavailable (bio-available space, keeps Rhodophyta)",
            'exclude_algae': True
        },
        {
            'counts': "CoralNet_phylum_counts_recruited_only.csv", 
            'props': "CoralNet_phylum_proportions_recruited_only.csv",
            'description': "Recruited organisms only (excluding No_recruitment, Unavailable, possibly algae)",
            'exclude_algae': False  # Depends on flags used
        }
    ]
    
    results = {}
    
    for file_info in file_pairs:
        counts_file = file_info['counts']
        props_file = file_info['props']
        
        if os.path.exists(counts_file) and os.path.exists(props_file):
            print(f"Checking {file_info['description']}...")
            valid = validate_proportion_files(
                counts_file,
                props_file,
                exclude_algae=file_info['exclude_algae'],
                description=file_info['description']
            )
            results[file_info['description']] = valid
            print("\n" + "="*60)
        else:
            print(f"Skipping {file_info['description']} - files not found")
            print(f"  Looking for: {counts_file}, {props_file}")
    
    # Summary
    print("\n" + "="*50)
    print("VALIDATION SUMMARY")
    print("="*50)
    
    for description, is_valid in results.items():
        status = "✅ VALID" if is_valid else "❌ INVALID"
        print(f"{status}: {description}")
    
    if all(results.values()):
        print(f"\n🎉 All proportion files are valid!")
    else:
        print(f"\n⚠️  Some files have proportion issues - check individual results above")
    
    # Additional algae-specific validation
    if os.path.exists("CoralNet_phylum_counts_no_algae.csv"):
        print(f"\n" + "="*50)
        print("ALGAE EXCLUSION VALIDATION")
        print("="*50)
        
        counts_no_algae = pd.read_csv("CoralNet_phylum_counts_no_algae.csv", index_col=0)
        counts_standard = pd.read_csv("CoralNet_phylum_counts.csv", index_col=0)
        
        # Check that Chlorophyta and Ochrophyta are excluded but Rhodophyta is kept
        excluded_algae = ['Chlorophyta', 'Ochrophyta']
        kept_algae = ['Rhodophyta']
        
        print(f"Checking algae exclusion:")
        for algae in excluded_algae:
            if algae in counts_standard.columns and algae not in counts_no_algae.columns:
                print(f"  ✅ {algae} successfully excluded")
            elif algae in counts_no_algae.columns:
                print(f"  ❌ {algae} still present in no-algae file")
            else:
                print(f"  ℹ️  {algae} not found in original data")
        
        for algae in kept_algae:
            if algae in counts_standard.columns and algae in counts_no_algae.columns:
                print(f"  ✅ {algae} correctly retained")
            elif algae in counts_standard.columns and algae not in counts_no_algae.columns:
                print(f"  ❌ {algae} incorrectly excluded")
            else:
                print(f"  ℹ️  {algae} not found in data")
        
        # Calculate how much data was excluded
        if 'Chlorophyta' in counts_standard.columns or 'Ochrophyta' in counts_standard.columns:
            numeric_cols = counts_standard.columns[1:]  # Skip sitelocality
            total_points = counts_standard[numeric_cols].sum().sum()
            
            excluded_points = 0
            for algae in excluded_algae:
                if algae in counts_standard.columns:
                    excluded_points += counts_standard[algae].sum()
            
            excluded_percent = (excluded_points / total_points) * 100
            print(f"  Total algae points excluded: {excluded_points:,} ({excluded_percent:.1f}% of all data)")

if __name__ == "__main__":
    main()
