#!/usr/bin/env python3
"""
Simple Pathogenic Allele Detection for TRGT VCF Files

This program identifies pathogenic tandem repeat alleles in a single VCF file 
and outputs results to a TSV file.

Usage: python trgt_pathogenic_detector_simple.py <disease_loci.tsv> <sample.vcf> <output.tsv>
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

class TRGTPathogenicDetector:
    def __init__(self, disease_loci_file: str):
        """Initialize the detector with disease loci reference data."""
        self.disease_loci = self._load_disease_loci(disease_loci_file)
        self.skipped_genes = {'COMP', 'FOXL2', 'DMD', 'MUC1', 'VWA1', 'pre-MIR7-2', 'POLG'}
        
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
        """Extract repeat information from TRGT VCF record.

        This is more robust to STRchive/VAMOS-style VCFs (RU, ALTANNO_H1, ALTANNO_H2, LEN_H1, LEN_H2)
        while still supporting TRGT-style AL/MC fields in the sample columns.
        """
        try:
            info = vcf_record['info'] or {}
            motifs = []
            allele_lengths = []
            motif_counts = []
            
            # Prefer RU from INFO as motif list (VAMOS / STRchive)
            if 'RU' in info and info['RU']:
                motifs = [m for m in info['RU'].split(',') if m != '']
            
            # ALTANNO_H1 / ALTANNO_H2 contain motif count series like "0-0-1-0-..."
            for altanno_key in ('ALTANNO_H1', 'ALTANNO_H2'):
                if altanno_key in info and info[altanno_key] and info[altanno_key] != '.':
                    parts = info[altanno_key].split('-')
                    # convert to underscore-separated counts to be compatible with existing logic
                    # and keep as string (existing code expects strings like "1_0_2")
                    try:
                        nums = [str(int(x)) for x in parts if x != '' and x != '.']
                        if nums:
                            motif_counts.append('_'.join(nums))
                    except ValueError:
                        # non-numeric ALTANNO values: keep raw joined string as fallback
                        motif_counts.append('_'.join(parts))
            
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
                
                if 'MC' in sample_dict and sample_dict['MC'] != '.':
                    mc_values = sample_dict['MC'].split(',')
                    motif_counts.extend([x for x in mc_values if x != '.' and x != ''])
            except Exception:
                # ignore parsing issues here, we already used INFO-first strategy
                pass
            
            return {
                'motifs': motifs,
                'allele_lengths': allele_lengths,
                'motif_counts': motif_counts
            }
            
        except (ValueError, IndexError, KeyError) as e:
            logger.warning(f"Could not extract repeat information: {e}")
            return {'motifs': [], 'allele_lengths': [], 'motif_counts': []}
    
    def _find_matching_disease_locus(self, chrom: str, pos: int) -> Optional[pd.Series]:
        """Find matching disease locus for given chromosome and position."""
        chrom_matches = self.disease_loci[self.disease_loci['chrom'] == chrom]
        position_matches = chrom_matches[
            (chrom_matches['start'] <= pos) & (chrom_matches['stop'] >= pos)
        ]
        
        if len(position_matches) > 0:
            return position_matches.iloc[0]
        return None
    
    def _calculate_copy_numbers(self, repeat_info: Dict, disease_locus: pd.Series) -> List[int]:
        """Calculate copy numbers from TRGT data."""
        copy_numbers = []
        
        if repeat_info['motif_counts']:
            reference_motif = disease_locus.get('reference_motif_reference_orientation', '')
            pathogenic_motifs = disease_locus.get('pathogenic_motif_reference_orientation', '') or ''
            motifs = repeat_info.get('motifs', []) or []
            
            try:
                # Normalize pathogenic motif list
                pathogenic_motif_list = [m.strip() for m in pathogenic_motifs.split(',')] if ',' in pathogenic_motifs else [pathogenic_motifs.strip()]
                pathogenic_motif_list = [m for m in pathogenic_motif_list if m]
                
                for mc_string in repeat_info['motif_counts']:
                    # parse counts into int list (underscore-separated or single)
                    if '_' in mc_string:
                        try:
                            mc_list = [int(x) for x in mc_string.split('_')]
                        except ValueError:
                            # fallback: skip malformed entries
                            continue
                    else:
                        try:
                            mc_list = [int(mc_string)]
                        except ValueError:
                            continue
                    
                    total_pathogenic_count = 0
                    
                    # If we have multiple pathogenic motifs, sum counts at positions whose mapped motif
                    # matches any pathogenic motif. Map mc_list indices to motifs cyclically when mc_list
                    # is longer than motifs (common for ALTANNO vectors).
                    if pathogenic_motif_list and len(pathogenic_motif_list) > 0:
                        if motifs:
                            for j, count in enumerate(mc_list):
                                if count == 0:
                                    continue
                                motif_idx = j % len(motifs)
                                observed_motif = motifs[motif_idx] if motif_idx < len(motifs) else ''
                                is_pathogenic_motif = any(
                                    observed_motif.upper() == pm.upper() or
                                    observed_motif.upper() in pm.upper() or
                                    pm.upper() in observed_motif.upper()
                                    for pm in pathogenic_motif_list
                                )
                                if is_pathogenic_motif:
                                    total_pathogenic_count += count
                        else:
                            # No motif labels available: treat all counts as potentially pathogenic (conservative)
                            total_pathogenic_count = sum(mc_list)
                    
                    # If no explicit pathogenic motif list (or previous logic didn't compute), attempt to use
                    # reference_motif index mapping (single motif case)
                    if total_pathogenic_count == 0:
                        # try to find reference motif index among motifs
                        motif_index = None
                        if reference_motif and motifs:
                            for i, motif in enumerate(motifs):
                                if motif and reference_motif and motif.upper() == reference_motif.upper():
                                    motif_index = i
                                    break
                            if motif_index is None:
                                for i, motif in enumerate(motifs):
                                    if motif and reference_motif and (reference_motif.upper() in motif.upper() or motif.upper() in reference_motif.upper()):
                                        motif_index = i
                                        break
                        
                        if motif_index is not None and motifs:
                            # Sum counts at positions mapping to motif_index (consider cyclic mapping)
                            for j, count in enumerate(mc_list):
                                if (j % len(motifs)) == motif_index:
                                    total_pathogenic_count += count
                        else:
                            # fallback: take sum of mc_list (conservative)
                            total_pathogenic_count = sum(mc_list)
                    
                    copy_numbers.append(total_pathogenic_count)
                
                return copy_numbers
                
            except (ValueError, IndexError) as e:
                logger.warning(f"Error parsing MC field {repeat_info.get('motif_counts')}: {e}")
        
        # Secondary method: Calculate from AL or LEN_Hx
        if repeat_info.get('allele_lengths'):
            reference_motif = disease_locus.get('reference_motif_reference_orientation', '')
            motif_length = len(reference_motif) if reference_motif else 1
            
            for al in repeat_info['allele_lengths']:
                try:
                    copy_number = int(al) // motif_length
                    copy_numbers.append(copy_number)
                except Exception:
                    pass
            
            return copy_numbers
        
        return []
    
    def _is_pathogenic(self, repeat_info: Dict, disease_locus: pd.Series) -> Tuple[bool, List[int]]:
        """Determine if any alleles are pathogenic."""
        try:
            pathogenic_min = disease_locus['pathogenic_min']
            
            if pd.isna(pathogenic_min) or pathogenic_min == 'None':
                return False, []
            
            threshold = int(pathogenic_min)
            copy_numbers = self._calculate_copy_numbers(repeat_info, disease_locus)
            
            if not copy_numbers:
                return False, []
            
            pathogenic_alleles = [cn for cn in copy_numbers if cn >= threshold]
            return len(pathogenic_alleles) > 0, pathogenic_alleles
            
        except (ValueError, TypeError) as e:
            logger.warning(f"Error determining pathogenicity: {e}")
            return False, []
    
    def _extract_sample_name(self, file_path: str) -> str:
        """Extract sample name from file path. keeps base sample (before underscore)."""
        file_name = Path(file_path).stem
        if file_name.endswith('.vcf'):
            file_name = file_name[:-4]
        return file_name.split('_')[0]
    
    def _open_file(self, file_path: str):
        """Open file, handling compressed and uncompressed formats."""
        if file_path.endswith('.gz'):
            return gzip.open(file_path, 'rt', encoding='utf-8')
        else:
            return open(file_path, 'r', encoding='utf-8')
    
    def process_vcf_files(self, vcf_files: List[str], output_file: str):
        """Process one or more VCF files (e.g. two per sample) and output pathogenic findings to TSV."""
        pathogenic_results = []
        
        for vcf_file in vcf_files:
            sample_name = self._extract_sample_name(vcf_file)
            logger.info(f"Processing VCF file: {vcf_file} (Sample: {sample_name})")
            
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
                        
                        # Extract repeat info, prefer INFO fields (ALTANNO_H1/H2, RU, LEN_Hx), fallback to sample FORMAT
                        repeat_info = self._extract_repeat_info(vcf_record)
                        
                        if repeat_info['allele_lengths'] or repeat_info['motif_counts']:
                            is_pathogenic, pathogenic_alleles = self._is_pathogenic(repeat_info, disease_locus)
                            
                            if is_pathogenic:
                                copy_numbers = self._calculate_copy_numbers(repeat_info, disease_locus)
                                
                                result = {
                                    'chromosome': vcf_record['chrom'],
                                    'position': vcf_record['pos'],
                                    'gene': disease_locus['gene'],
                                    'disease': disease_locus['disease'],
                                    'inheritance': disease_locus['inheritance'],
                                    'reference_motif': disease_locus['reference_motif_reference_orientation'],
                                    'observed_motifs': ','.join(repeat_info['motifs']),
                                    'pathogenic_threshold': disease_locus['pathogenic_min'],
                                    'allele_lengths_bp': ','.join(map(str, repeat_info['allele_lengths'])),
                                    'calculated_copy_numbers': ','.join(map(str, copy_numbers)),
                                    'filter': vcf_record['filter']
                                }
                                
                                pathogenic_results.append(result)
                                
                                logger.info(f"PATHOGENIC: {gene_name} ({disease_locus['disease']}) - "
                                          f"Copy numbers: {pathogenic_alleles} (threshold: {disease_locus['pathogenic_min']})")
            except FileNotFoundError:
                logger.error(f"VCF file not found: {vcf_file}")
            except Exception as e:
                logger.error(f"Error processing VCF file {vcf_file}: {e}")
        
        # Write results to TSV (single combined output)
        if pathogenic_results:
            df = pd.DataFrame(pathogenic_results)
            df.to_csv(output_file, sep='\t', index=False)
            logger.info(f"Found {len(pathogenic_results)} pathogenic alleles")
            logger.info(f"Results written to: {output_file}")
        else:
            # Create empty file with headers
            columns = ['sample', 'vcf_file', 'chromosome', 'position', 'gene', 'disease', 'inheritance',
                      'reference_motif', 'pathogenic_threshold', 'observed_motifs', 
                      'allele_lengths_bp', 'calculated_copy_numbers', 'pathogenic_alleles',
                      'vcf_id', 'quality', 'filter']
            df = pd.DataFrame(columns=columns)
            df.to_csv(output_file, sep='\t', index=False)
            logger.info("No pathogenic alleles found")
            logger.info(f"Empty results file created: {output_file}")

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
        detector = TRGTPathogenicDetector(disease_loci_file)
        detector.process_vcf_files(vcf_files, output_file)
        print(f"✅ Analysis complete. Results saved to: {output_file}")
    except Exception as e:
        logger.error(f"Failed to process: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
