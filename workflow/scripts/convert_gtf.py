#!/usr/bin/env python3
import argparse

def parse_attributes(attr_str):
    """Parse the attributes column (9th column) into a dictionary."""
    attrs = {}
    for attr in attr_str.strip().split('; '):
        if attr:
            key, value = attr.split(' ', 1)
            attrs[key] = value.strip('"')
    return attrs

def process_isoform(isoform_str):
    """Process isoform string to get exon coordinates.
    Returns list of (start, end) tuples."""
    # Clean the isoform string by removing quotes and semicolons
    isoform_str = isoform_str.strip('"').strip(';')
    parts = isoform_str.split(',')
    
    exons = []
    for part in parts:
        # Clean each part of any remaining quotes
        part = part.strip('"')
        start, end = map(int, part.split('-'))
        # GTF is 1-based, so we don't need to adjust the coordinates
        exons.append((start, end))
    
    return exons

def format_gtf_attributes(attrs_dict):
    """Format attributes for GTF format."""
    formatted = []
    for key, value in attrs_dict.items():
        formatted.append(f'{key} "{value}"')
    return "; ".join(formatted)

def convert_to_gtf(input_file, output_file):
    """Convert CIRI-long .info file to GTF3 format."""
    with open(input_file, 'r') as fin, open(output_file, 'w') as fout:
        for line in fin:
            if line.startswith('#') or not line.strip():
                continue
                
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
                
            chrom = fields[0]
            start = int(fields[3])  # Start position
            end = int(fields[4])    # End position
            score = fields[5]  # Score/depth
            strand = fields[6] if fields[6] != 'None' else '.'
            
            # Parse attributes
            attrs = parse_attributes(fields[8])
            circ_id = attrs.get('circ_id', '.').strip('"')
            gene_id = attrs.get('gene_id', 'novel_gene').strip('"')
            gene_name = attrs.get('gene_name', 'novel_gene').strip('"')
            
            # Process isoform information
            isoform = attrs.get('isoform', '')
            
            if '|' in isoform:
                # Handle multiple isoforms
                isoform = isoform.strip('"').strip(';')
                isoform_list = isoform.split('|')
            else:
                isoform_list = [isoform] if isoform else [(f"{start}-{end}")]

            # Process each isoform
            for i, iso in enumerate(isoform_list):
                transcript_id = f"{circ_id}_iso{i+1}"
                
                # Write transcript entry
                transcript_attrs = {
                    'gene_id': gene_id,
                    'transcript_id': transcript_id,
                    'gene_name': gene_name,
                    'transcript_biotype': 'circRNA'
                }
                fout.write(f'{chrom}\tCIRI-long\ttranscript\t{start}\t{end}\t{score}\t{strand}\t.\t{format_gtf_attributes(transcript_attrs)}\n')
                
                # Process and write exon entries
                exons = process_isoform(iso)
                for j, (exon_start, exon_end) in enumerate(exons, 1):
                    exon_attrs = {
                        'gene_id': gene_id,
                        'transcript_id': transcript_id,
                        'gene_name': gene_name,
                        'exon_number': str(j),
                        'exon_id': f"{transcript_id}.{j}"
                    }
                    fout.write(f'{chrom}\tCIRI-long\texon\t{exon_start}\t{exon_end}\t{score}\t{strand}\t.\t{format_gtf_attributes(exon_attrs)}\n')

def main():
    parser = argparse.ArgumentParser(
        description='Convert CIRI-long .info file to GTF3 format',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument(
        '-i', '--input',
        required=True,
        help='Input CIRI-long .info file'
    )
    
    parser.add_argument(
        '-o', '--output',
        required=True,
        help='Output GTF3 file'
    )
    
    args = parser.parse_args()
    convert_to_gtf(args.input, args.output)

if __name__ == '__main__':
    main()
