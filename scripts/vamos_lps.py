#!/usr/bin/env python3
"""
VCF VNTR Longest Repeat Analyzer

This script processes VCF files containing VNTR data and identifies the longest
pure repeat segments for each locus based on the ALTANNO_H1 field.
Can process single files or entire directories of compressed VCFs.
"""

import re
import gzip
import os
import glob
import argparse
from collections import Counter

def open_file(filename):
    """Open file, handling both compressed and uncompressed formats."""
    if filename.endswith('.gz'):
        return gzip.open(filename, 'rt')
    else:
        return open(filename, 'r')

def extract_sample_ids(header_line):
    """Extract sample IDs from VCF header line starting with #CHROM."""
    fields = header_line.strip().split('\t')
    if len(fields) > 9:  # FORMAT is at index 8, samples start at index 9
        return fields[9:]
    return []

def parse_vcf_line(line, sample_ids):
    """Parse a VCF data line and extract relevant fields."""
    if line.startswith('#') or not line.strip():
        return None
    
    fields = line.strip().split('\t')
    if len(fields) < 8:
        return None
    
    chrom = fields[0]
    pos = fields[1]
    info = fields[7]
    
    # Parse INFO field
    info_dict = {}
    for item in info.split(';'):
        if '=' in item:
            key, value = item.split('=', 1)
            info_dict[key] = value
    
    # Extract sample data if available
    sample_data = {}
    if len(fields) > 9 and sample_ids:
        format_field = fields[8] if len(fields) > 8 else ""
        for i, sample_id in enumerate(sample_ids):
            if i + 9 < len(fields):
                sample_data[sample_id] = fields[i + 9]
    
    return {
        'CHROM': chrom,
        'POS': pos,
        'INFO': info_dict,
        'SAMPLES': sample_data
    }

def find_longest_pure_repeat(altanno_sequence, motifs):
    """
    Find the longest pure repeat in the ALTANNO sequence.
    
    Args:
        altanno_sequence: String of motif indices separated by '-'
        motifs: List of motif sequences from RU field
    
    Returns:
        Tuple of (motif_sequence, repeat_count) for longest pure repeat
    """
    if not altanno_sequence or not motifs:
        return None, 0
    
    # Split the sequence into individual motif indices
    indices = altanno_sequence.split('-')
    
    # Convert to integers and handle any parsing errors
    try:
        indices = [int(x) for x in indices if x.strip()]
    except ValueError:
        return None, 0
    
    if not indices:
        return None, 0
    
    # Find longest consecutive runs of the same motif index
    max_count = 0
    max_motif_idx = 0
    current_count = 1
    current_idx = indices[0]
    
    for i in range(1, len(indices)):
        if indices[i] == current_idx:
            current_count += 1
        else:
            if current_count > max_count:
                max_count = current_count
                max_motif_idx = current_idx
            current_count = 1
            current_idx = indices[i]
    
    # Check the last run
    if current_count > max_count:
        max_count = current_count
        max_motif_idx = current_idx
    
    # Get the corresponding motif sequence
    if 0 <= max_motif_idx < len(motifs):
        return motifs[max_motif_idx], max_count
    else:
        return None, 0

def process_vcf_file(filename):
    """
    Process a VCF file and extract longest pure repeats for each VNTR locus.
    
    Args:
        filename: Path to the VCF file
    
    Returns:
        Tuple of (results_list, sample_ids)
    """
    results = []
    sample_ids = []
    
    with open_file(filename) as f:
        for line in f:
            # Extract sample IDs from header
            if line.startswith('#CHROM'):
                sample_ids = extract_sample_ids(line)
                continue
            
            parsed = parse_vcf_line(line, sample_ids)
            if not parsed:
                continue
            
            info = parsed['INFO']
            
            # Check if this is a VNTR entry with required fields
            if 'ALTANNO_H1' not in info or 'RU' not in info or 'END' not in info:
                continue
            
            # Parse motifs from RU field
            ru_field = info['RU']
            motifs = [motif.strip() for motif in ru_field.split(',')]
            
            # Get ALTANNO_H1 sequence
            altanno = info['ALTANNO_H1']
            
            # Find longest pure repeat
            longest_motif, repeat_count = find_longest_pure_repeat(altanno, motifs)
            
            if longest_motif and repeat_count > 0:
                # Process each sample
                for sample_id in sample_ids:
                    result = {
                        'CHROM': parsed['CHROM'],
                        'POS': int(parsed['POS']),
                        'END': int(info['END']),
                        'SAMPLE_ID': sample_id,
                        'LONGEST_MOTIF': longest_motif,
                        'REPEAT_COUNT': repeat_count,
                        'FILENAME': os.path.basename(filename)
                    }
                    results.append(result)
    
    return results, sample_ids

def process_directory(directory_path, pattern="*.vcf.gz"):
    """
    Process all VCF files in a directory matching the given pattern.
    
    Args:
        directory_path: Path to directory containing VCF files
        pattern: File pattern to match (default: *.vcf.gz)
    
    Returns:
        List of all results from all files
    """
    all_results = []
    file_pattern = os.path.join(directory_path, pattern)
    vcf_files = glob.glob(file_pattern)
    
    if not vcf_files:
        # Try uncompressed VCFs as well
        uncompressed_pattern = os.path.join(directory_path, "*.vcf")
        vcf_files = glob.glob(uncompressed_pattern)
    
    if not vcf_files:
        print(f"No VCF files found in {directory_path} matching pattern {pattern}")
        return []
    
    print(f"Found {len(vcf_files)} VCF files to process:")
    for vcf_file in vcf_files:
        print(f"  {os.path.basename(vcf_file)}")
    
    for vcf_file in vcf_files:
        print(f"\nProcessing: {os.path.basename(vcf_file)}")
        try:
            results, sample_ids = process_vcf_file(vcf_file)
            all_results.extend(results)
            print(f"  Found {len(results)} VNTR loci")
        except Exception as e:
            print(f"  Error processing {vcf_file}: {e}")
    
    return all_results

def write_bed_file(results, output_file):
    """
    Write results to a BED format file.
    BED format: chrom start end name score strand
    """
    # Sort results by chromosome and position
    results.sort(key=lambda x: (x['CHROM'], x['POS']))
    
    with open(output_file, 'w') as f:
        # Write header as comment
        f.write("# BED format: chrom start end name score strand\n")
        f.write("# name format: SAMPLE_ID|MOTIF|COUNT\n")
        
        for result in results:
            # BED uses 0-based coordinates, VCF uses 1-based
            start = result['POS'] - 1
            end = result['END']
            name = f"{result['SAMPLE_ID']}|{result['LONGEST_MOTIF']}|{result['REPEAT_COUNT']}"
            score = result['REPEAT_COUNT']  # Use repeat count as score
            strand = "."
            
            f.write(f"{result['CHROM']}\t{start}\t{end}\t{name}\t{score}\t{strand}\n")

def print_summary(results):
    """Print summary statistics."""
    if not results:
        print("No VNTR loci found.")
        return
    
    print(f"\nSummary:")
    print(f"Total VNTR loci processed: {len(results)}")
    
    # Sample statistics
    sample_counts = Counter(result['SAMPLE_ID'] for result in results)
    print(f"Samples found: {len(sample_counts)}")
    for sample, count in sample_counts.items():
        print(f"  {sample}: {count} loci")
    
    # Motif statistics
    motif_counts = Counter(result['LONGEST_MOTIF'] for result in results)
    print(f"\nMost common longest motifs:")
    for motif, count in motif_counts.most_common(10):
        print(f"  {motif}: {count} loci")
    
    # Chromosome distribution
    chrom_counts = Counter(result['CHROM'] for result in results)
    print(f"\nChromosome distribution:")
    for chrom in sorted(chrom_counts.keys(), key=lambda x: (x.replace('chr', ''), x)):
        print(f"  {chrom}: {chrom_counts[chrom]} loci")

def main():
    """Main function to run the analysis."""
    parser = argparse.ArgumentParser(
        description="Analyze VNTR longest repeats from VCF files",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Process single VCF file
  python vntr_analyzer.py -i sample.vcf.gz -o results.bed
  
  # Process all compressed VCFs in directory
  python vntr_analyzer.py -d /path/to/vcfs/ -o all_vntrs.bed
  
  # Process all VCFs (compressed and uncompressed) in directory
  python vntr_analyzer.py -d /path/to/vcfs/ -p "*.vcf*" -o results.bed
        """
    )
    
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-i', '--input', help='Input VCF file')
    group.add_argument('-d', '--directory', help='Directory containing VCF files')
    
    parser.add_argument('-o', '--output', required=True, 
                       help='Output BED file')
    parser.add_argument('-p', '--pattern', default='*.vcf.gz',
                       help='File pattern for directory mode (default: *.vcf.gz)')
    parser.add_argument('--verbose', action='store_true',
                       help='Print detailed progress information')
    
    args = parser.parse_args()
    
    try:
        if args.input:
            # Process single file
            print(f"Processing VCF file: {args.input}")
            results, sample_ids = process_vcf_file(args.input)
        else:
            # Process directory
            print(f"Processing directory: {args.directory}")
            print(f"File pattern: {args.pattern}")
            results = process_directory(args.directory, args.pattern)
        
        if not results:
            print("No VNTR data found to process.")
            return
        
        # Write BED file
        write_bed_file(results, args.output)
        print(f"\nResults written to BED file: {args.output}")
        
        # Print summary
        print_summary(results)
        
    except FileNotFoundError as e:
        print(f"Error: {e}")
    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    main()