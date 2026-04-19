#!/usr/bin/env python3
"""
Simple Pathogenic Allele Detection for vamos VCF Files

This program identifies pathogenic tandem repeat alleles in a single VCF file 
and outputs results to a TSV file.

Usage: python pathogenic_detector.py <disease_loci.tsv> <sample.vcf> <output.tsv>
"""

import sys
import pandas as pd
import logging
from pathlib import Path
from typing import Dict, List, Tuple, Optional
import gzip

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class PathogenicDetector:
    def __init__(self, disease_loci_file: str):
        """Initialize the detector with disease loci reference data."""
        self.disease_loci = self._load_disease_loci(disease_loci_file)
        self.skipped_genes = {'COMP', 'FOXL2', 'DMD', 'MUC1', 'VWA1', 'MIR7-2', 'POLG', 'XYLT1', 'ZIC3', 'NIPA1'}
        
    def _load_disease_loci(self, file_path: str) -> pd.DataFrame:
        """Load and process the disease loci reference file."""
        try:
            df = pd.read_csv(file_path, sep='\t')
            logger.info(f"Loaded {len(df)} disease loci from {file_path}")
            
            # Get the actual chromosome column name (might be '#chrom' or 'chrom')
            chrom_col = None
            for col in df.columns:
                if 'chrom' in col.lower():
                    chrom_col = col
                    break
            
            if chrom_col is None:
                raise ValueError("Could not find chromosome column in the file")
            
            # Rename the chromosome column to standardize it
            if chrom_col != 'chrom':
                df = df.rename(columns={chrom_col: 'chrom'})
            
            return df
        except Exception as e:
            logger.error(f"Error loading disease loci file: {e}")
            raise
    
    def _parse_vcf_line(self, line: str) -> Optional[Dict]:
        """Parse a VCF line and extract relevant information."""
        if line.startswith('#') or not line.strip():
            return None
            
        fields = line.strip().split('\t')
        if len(fields) < 9:
            return None
            
        chrom, pos, vcf_id, ref, alt, qual, filter_field, info, format_field = fields[:9]
        sample_data = fields[9:] if len(fields) > 9 else []
        
        # Parse INFO field
        info_dict = {}
        for item in info.split(';'):
            if '=' in item:
                key, value = item.split('=', 1)
                info_dict[key] = value
        
        return {
            'chrom': chrom,
            'pos': int(pos),
            'id': vcf_id,
            'ref': ref,
            'alt': alt,
            'qual': qual,
            'filter': filter_field,
            'info': info_dict,
            'format': format_field,
            'samples': sample_data
        }
    
    def _extract_repeat_info(self, vcf_record: Dict, sample_idx: int = 0) -> Dict:
        """Extract repeat information from vamos VCF record.

        
        For ALTANNO fields: values are INDICES into the RU (motif) list, not repeat counts.
        """
        try:
            info = vcf_record['info'] or {}
            motifs = []
            allele_lengths = []
            altanno_sequences = []  # Store ALTANNO as index sequences
            
            # Prefer RU from INFO as motif list (VAMOS / STRchive)
            if 'RU' in info and info['RU']:
                motifs = [m for m in info['RU'].split(',') if m != '']
            
            # ALTANNO_H1 / ALTANNO_H2 contain motif INDEX sequences like "0-3-4-4-4-..."
            # Each number is an index into the RU motif list
            for altanno_key in ('ALTANNO_H1', 'ALTANNO_H2'):
                if altanno_key in info and info[altanno_key] and info[altanno_key] != '.':
                    try:
                        # Parse as hyphen-separated indices
                        indices = [int(x) for x in info[altanno_key].split('-') if x != '' and x != '.']
                        if indices:
                            altanno_sequences.append(indices)
                    except ValueError:
                        # non-numeric ALTANNO values: skip
                        pass
            
            # LEN_H1 / LEN_H2 provide allele length in bp (secondary info)
            for len_key in ('LEN_H1', 'LEN_H2'):
                if len_key in info and info[len_key] and info[len_key] != '.':
                    try:
                        allele_lengths.append(int(info[len_key]))
                    except ValueError:
                        pass
            
            # If sample FORMAT contains AL or MC, use those as fallback
            try:
                format_fields = vcf_record['format'].split(':') if vcf_record.get('format') else []
                sample_field = vcf_record['samples'][sample_idx] if vcf_record.get('samples') else ''
                sample_parts = sample_field.split(':') if sample_field else []
                sample_dict = dict(zip(format_fields, sample_parts))
                
                if not motifs and 'MOTIFS' in vcf_record['info']:
                    motifs = [m for m in vcf_record['info'].get('MOTIFS', '').split(',') if m != '']
                
                if 'AL' in sample_dict and sample_dict['AL'] != '.':
                    al_values = sample_dict['AL'].split(',')
                    for x in al_values:
                        try:
                            allele_lengths.append(int(x))
                        except ValueError:
                            pass
                
                # MC field format is different - these are repeat counts, not indices
                # We'll handle this separately if ALTANNO is not available
            except Exception:
                pass
            
            return {
                'motifs': motifs,
                'allele_lengths': allele_lengths,
                'altanno_sequences': altanno_sequences  # List of index sequences
            }
            
        except (ValueError, IndexError, KeyError) as e:
            logger.warning(f"Could not extract repeat information: {e}")
            return {'motifs': [], 'allele_lengths': [], 'altanno_sequences': []}
    
    def _find_matching_disease_locus(self, chrom: str, pos: int) -> Optional[pd.Series]:
        """Find matching disease locus for given chromosome and position."""
        chrom_matches = self.disease_loci[self.disease_loci['chrom'] == chrom]
        position_matches = chrom_matches[
            (chrom_matches['start'] <= pos) & (chrom_matches['stop'] >= pos)
        ]
        
        if len(position_matches) > 0:
            return position_matches.iloc[0]
        return None
    
    def _calculate_copy_numbers(self, repeat_info: Dict, disease_locus: pd.Series) -> Tuple[List[int], List[str]]:
        """Calculate copy numbers from vamos data and return detected pathogenic motifs.
        
        For ALTANNO format: each value is an INDEX into the RU motif list.
        Count how many times pathogenic motif indices appear in the sequence.
        
        Returns:
            Tuple of (copy_numbers, detected_pathogenic_motifs)
        """
        copy_numbers = []
        detected_motifs = []
        
        # Get pathogenic motif information
        pathogenic_motifs_str = disease_locus.get('pathogenic_motif_reference_orientation', '') or ''
        pathogenic_motif_list = [m.strip() for m in pathogenic_motifs_str.split(',')] if ',' in pathogenic_motifs_str else [pathogenic_motifs_str.strip()]
        pathogenic_motif_list = [m for m in pathogenic_motif_list if m and m.upper() != 'NONE']
        
        motifs = repeat_info.get('motifs', []) or []
        altanno_sequences = repeat_info.get('altanno_sequences', [])
        
        # Primary method: Use ALTANNO sequences (indices into RU/motifs list)
        if altanno_sequences and motifs:
            # Build set of pathogenic indices by matching motifs to pathogenic list
            pathogenic_indices = set()
            for i, motif in enumerate(motifs):
                if any(motif.upper() == pm.upper() for pm in pathogenic_motif_list):
                    pathogenic_indices.add(i)
            
            if not pathogenic_indices:
                # No pathogenic motifs found in the motif list
                logger.warning(f"No pathogenic motifs found in RU motif list")
            
            # Count occurrences of pathogenic indices in each ALTANNO sequence
            for idx_sequence in altanno_sequences:
                pathogenic_count = 0
                allele_detected_motifs = set()
                
                for idx in idx_sequence:
                    if idx in pathogenic_indices and idx < len(motifs):
                        pathogenic_count += 1
                        allele_detected_motifs.add(motifs[idx])
                
                copy_numbers.append(pathogenic_count)
                detected_motifs.append(','.join(sorted(allele_detected_motifs)))
            
            return copy_numbers, detected_motifs
        
        # Fallback: Calculate from allele lengths
        if repeat_info.get('allele_lengths'):
            reference_motif = disease_locus.get('reference_motif_reference_orientation', '')
            pathogenic_motif = disease_locus.get('pathogenic_motif_reference_orientation', '')
            motif_length = len(reference_motif) if reference_motif else 1
            
            for al in repeat_info['allele_lengths']:
                try:
                    copy_number = int(al) // motif_length
                    copy_numbers.append(copy_number)
                    # Use pathogenic motif if available, otherwise reference motif
                    detected_motifs.append(pathogenic_motif if pathogenic_motif else reference_motif)
                except Exception:
                    pass
            
            return copy_numbers, detected_motifs
        
        return [], []
    
    def _is_pathogenic(self, repeat_info: Dict, disease_locus: pd.Series) -> Tuple[bool, List[int], List[str]]:
        """Determine if any alleles are pathogenic.
        
        Returns:
            Tuple of (is_pathogenic, pathogenic_allele_counts, detected_motifs)
        """
        try:
            pathogenic_min = disease_locus['pathogenic_min']
            
            if pd.isna(pathogenic_min) or pathogenic_min == 'None':
                return False, [], []
            
            threshold = int(pathogenic_min)
            copy_numbers, detected_motifs = self._calculate_copy_numbers(repeat_info, disease_locus)
            
            if not copy_numbers:
                return False, [], []
            
            pathogenic_alleles = []
            pathogenic_motifs = []
            for i, cn in enumerate(copy_numbers):
                if cn >= threshold:
                    pathogenic_alleles.append(cn)
                    if i < len(detected_motifs):
                        pathogenic_motifs.append(detected_motifs[i])
            
            return len(pathogenic_alleles) > 0, pathogenic_alleles, pathogenic_motifs
            
        except (ValueError, TypeError) as e:
            logger.warning(f"Error determining pathogenicity: {e}")
            return False, [], []
    
    def _open_file(self, file_path: str):
        """Open file, handling compressed and uncompressed formats."""
        if file_path.endswith('.gz'):
            return gzip.open(file_path, 'rt', encoding='utf-8')
        else:
            return open(file_path, 'r', encoding='utf-8')
    
    def process_vcf_files(self, vcf_files: List[str], output_file: str):
        """Process one or more VCF files and output pathogenic findings to TSV."""
        pathogenic_results = []
        
        for vcf_file in vcf_files:
            logger.info(f"Processing VCF file: {vcf_file}")
            
            try:
                with self._open_file(vcf_file) as f:
                    for line in f:
                        vcf_record = self._parse_vcf_line(line)
                        if not vcf_record:
                            continue
                        
                        disease_locus = self._find_matching_disease_locus(vcf_record['chrom'], vcf_record['pos'])
                        if disease_locus is None:
                            continue
                        
                        gene_name = disease_locus['gene']
                        if gene_name in self.skipped_genes:
                            continue
                        
                        # Extract repeat info
                        repeat_info = self._extract_repeat_info(vcf_record)
                        
                        if repeat_info['allele_lengths'] or repeat_info['altanno_sequences']:
                            is_pathogenic, pathogenic_alleles, pathogenic_motifs = self._is_pathogenic(repeat_info, disease_locus)
                            
                            if is_pathogenic:
                                locus = f"{disease_locus['chrom']}:{int(disease_locus['start'])}-{int(disease_locus['stop'])}"
                                
                                for i, cn in enumerate(pathogenic_alleles):
                                    detected_motif = pathogenic_motifs[i] if i < len(pathogenic_motifs) else ''
                                    
                                    result = {
                                        'locus': locus,
                                        'gene': disease_locus['gene'],
                                        'inheritance': disease_locus['inheritance'],
                                        'condition': disease_locus['disease'],
                                        'motif': detected_motif,
                                        'count': int(cn),
                                        'threshold': int(disease_locus['pathogenic_min'])
                                    }
                                    pathogenic_results.append(result)

                                logger.info(f"PATHOGENIC: {gene_name} ({disease_locus['disease']}) - Count: {pathogenic_alleles} - Motifs: {pathogenic_motifs} (threshold: {disease_locus['pathogenic_min']})")
            except FileNotFoundError:
                logger.error(f"VCF file not found: {vcf_file}")
            except Exception as e:
                logger.error(f"Error processing VCF file {vcf_file}: {e}")
        
        # Write results to TSV
        columns = ['locus', 'gene', 'inheritance', 'condition', 'motif', 'count', 'threshold']
        df = pd.DataFrame(pathogenic_results, columns=columns)
        df.to_csv(output_file, sep='\t', index=False)
        
        if len(pathogenic_results) > 0:
            logger.info(f"Found {len(pathogenic_results)} pathogenic alleles")
        else:
            logger.info("No pathogenic alleles found")
        
        logger.info(f"Results written to: {output_file}")

def main():
    # New usage: at least 3 args: disease_loci.tsv, one-or-more vcf files, output.tsv
    if len(sys.argv) < 4:
        print("Usage: python strchive_pathogenic_detector.py <disease_loci.tsv> <sample1.vcf> [sample2.vcf ...] <output.tsv>")
        print("")
        print("Arguments:")
        print("  disease_loci.tsv  - Disease loci reference file")
        print("  sample.vcf        - One or more input VCF files (can be .vcf or .vcf.gz). Two VCFs per sample (hp1/hp2) are supported.")
        print("  output.tsv        - Output TSV file for pathogenic results")
        print("")
        print("Example:")
        print("  python strchive_pathogenic_detector.py disease_loci.tsv GM19238_hp1.vcf GM19238_hp2.vcf pathogenic_results.tsv")
        sys.exit(1)
    
    disease_loci_file = sys.argv[1]
    vcf_files = sys.argv[2:-1]
    output_file = sys.argv[-1]
    
    # Validate input files
    if not Path(disease_loci_file).exists():
        logger.error(f"Disease loci file not found: {disease_loci_file}")
        sys.exit(1)
    
    for v in vcf_files:
        if not Path(v).exists():
            logger.error(f"VCF file not found: {v}")
            sys.exit(1)
    
    # Initialize detector and process
    try:
        detector = PathogenicDetector(disease_loci_file)
        detector.process_vcf_files(vcf_files, output_file)
        print(f"✅ Analysis complete. Results saved to: {output_file}")
    except Exception as e:
        logger.error(f"Failed to process: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()