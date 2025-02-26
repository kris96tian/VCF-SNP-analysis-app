import os
import re
import gzip
import logging
import tempfile

logger = logging.getLogger(__name__)

def detect_file_format(filepath):
    is_gzipped = filepath.endswith('.gz')
    
    if filepath.endswith('.vcf.gz') or filepath.endswith('.vcf'):
        return is_gzipped, 'vcf'
    elif filepath.endswith('.bcf'):
        return False, 'bcf'
    
    try:
        if is_gzipped:
            with gzip.open(filepath, 'rt') as f:
                first_line = f.readline()
        else:
            with open(filepath, 'r') as f:
                first_line = f.readline()
        
        if first_line.startswith('##fileformat=VCF'):
            return is_gzipped, 'vcf'
        elif first_line.startswith('BCF'):
            return False, 'bcf'
        else:
            return is_gzipped, 'unknown'
    except Exception as e:
        logger.warning(f"Error detecting file format: {str(e)}")
        return is_gzipped, 'unknown'


def fix_vcf_header(filepath, output_filepath=None):

    if output_filepath is None:
        output_filepath = f"{filepath}.fixed"
    
    is_compressed, format_type = detect_file_format(filepath)
    
    if format_type != 'vcf':
        raise ValueError(f"Cannot fix non-VCF file format: {format_type}")
    
    try:
        infile = gzip.open(filepath, 'rt') if is_compressed else open(filepath, 'r')
        outfile = open(output_filepath, 'w')
        for line in infile:
            if line.startswith('#CHROM'):
                parts = re.split(r'\s+', line.strip())
                fixed_line = '\t'.join(parts) + '\n'
                outfile.write(fixed_line)
            elif line.startswith('##'):
                if not re.match(r'^##\w+=', line):
                    fixed_line = re.sub(r'^##\s*(\w+)\s*=', r'##\1=', line)
                    outfile.write(fixed_line)
                else:
                    outfile.write(line)
            else:
                if ' ' in line and '\t' not in line:
                    parts = re.split(r'\s+', line.strip())
                    fixed_line = '\t'.join(parts) + '\n'
                    outfile.write(fixed_line)
                else:
                    outfile.write(line)        
        infile.close()
        outfile.close()
        
        return output_filepath
    except Exception as e:
        logger.error(f"Error fixing VCF file: {str(e)}")
        if infile:
            infile.close()
        if outfile:
            outfile.close()
        raise
