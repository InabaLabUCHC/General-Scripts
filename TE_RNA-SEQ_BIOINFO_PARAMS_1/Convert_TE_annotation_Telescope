#!/usr/bin/env python3
"""
Convert RepeatMasker GFF3 output to Telescope-compatible GTF format
Usage: python3 repeatmasker_to_gtf.py genome.fa.out genome.fa
"""

import sys
import re
from collections import defaultdict

def parse_repeatmasker_out(out_file):
    """Parse RepeatMasker .out file"""
    tes = []
    
    with open(out_file, 'r') as f:
        # Skip header lines
        for _ in range(3):
            next(f)
        
        locus_counter = defaultdict(int)
        
        for line in f:
            line = line.strip()
            if not line:
                continue
            
            # RepeatMasker .out format (space-delimited)
            parts = line.split()
            if len(parts) < 15:
                continue
            
            # Extract fields
            score = parts[0]
            perc_div = parts[1]
            perc_del = parts[2]
            perc_ins = parts[3]
            query_seq = parts[4]
            query_start = int(parts[5])
            query_end = int(parts[6])
            query_left = parts[7]
            strand = parts[8]
            repeat_name = parts[9]
            repeat_class = parts[10]
            
            # Convert strand
            if strand == 'C':
                strand = '-'
            else:
                strand = '+'
            
            # Create locus ID (format: repeat_name{}counter)
            locus_counter[repeat_name] += 1
            locus_id = f"{repeat_name}{'{}'}{locus_counter[repeat_name]}"
            
            # Generate FBti-style ID (pseudo)
            fbti_id = f"TE{abs(hash(locus_id)) % 10000000:07d}"
            
            tes.append({
                'chr': query_seq,
                'start': query_start,
                'end': query_end,
                'strand': strand,
                'family': repeat_name,
                'class': repeat_class,
                'gene_id': fbti_id,
                'transcript_id': fbti_id,
                'locus_id': locus_id,
                'length': query_end - query_start + 1,
                'perc_div': perc_div,
                'score': score
            })
    
    return tes

def write_gtf(tes, output_file):
    """Write GTF file in Telescope-compatible format"""
    
    with open(output_file, 'w') as f:
        # Write header
        f.write("# GTF annotation of Drosophila melanogaster transposable elements\n")
        f.write("# Generated from RepeatMasker output\n")
        f.write("# Format: GTF (compatible with Telescope)\n")
        f.write("# Overlap mode: union (handles nested TEs)\n")
        f.write("#\n")
        
        # Write TE entries
        for te in tes:
            # GTF format: chr source feature start end score strand frame attributes
            attributes = (
                f'gene_id "{te["gene_id"]}"; '
                f'transcript_id "{te["transcript_id"]}"; '
                f'family_id "{te["family"]}"; '
                f'locus_id "{te["locus_id"]}"; '
                f'length "{te["length"]}"; '
                f'repeat_class "{te["class"]}"; '
                f'divergence "{te["perc_div"]}";'
            )
            
            gtf_line = "\t".join([
                te['chr'],
                'RepeatMasker',
                'exon',
                str(te['start']),
                str(te['end']),
                te['score'],
                te['strand'],
                '.',
                attributes
            ])
            
            f.write(gtf_line + "\n")

def main():
    if len(sys.argv) != 3:
        print("Usage: python3 repeatmasker_to_gtf.py genome.fa.out output.gtf")
        sys.exit(1)
    
    out_file = sys.argv[1]
    output_gtf = sys.argv[2]
    
    print(f"Parsing RepeatMasker output: {out_file}")
    tes = parse_repeatmasker_out(out_file)
    
    print(f"Found {len(tes)} TE annotations")
    
    # Count families
    families = {}
    for te in tes:
        families[te['family']] = families.get(te['family'], 0) + 1
    
    print(f"Found {len(families)} unique TE families")
    print("\nTop 10 families:")
    for family, count in sorted(families.items(), key=lambda x: x[1], reverse=True)[:10]:
        print(f"  {family}: {count}")
    
    print(f"\nWriting GTF to: {output_gtf}")
    write_gtf(tes, output_gtf)
    
    print("Done!")

if __name__ == "__main__":
    main()
