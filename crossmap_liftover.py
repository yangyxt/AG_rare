import argparse
import os
import subprocess

def liftover_vcf(input_vcf, chain_file, ref_fasta, output_vcf):
    cmd = [
        'CrossMap.py', 'vcf',
        chain_file,
        input_vcf,
        ref_fasta,
        output_vcf
    ]
    subprocess.run(cmd, check=True)
    print(f"Liftover complete: {output_vcf}")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Liftover VCF from hg19/GRCh37 to hg38/GRCh38.')
    parser.add_argument('rare_vcf', help='Path to rare variants VCF (hg19)')
    parser.add_argument('unfiltered_vcf', help='Path to unfiltered variants VCF (hg19)')
    parser.add_argument('--chain_file', default='hg19ToHg38.over.chain.gz', help='Path to chain file')
    parser.add_argument('--ref_fasta', default='hg38.fa', help='Path to hg38 FASTA')
    args = parser.parse_args()

    # Liftover rare VCF
    rare_output = os.path.splitext(args.rare_vcf)[0] + '_hg38.vcf'
    liftover_vcf(args.rare_vcf, args.chain_file, args.ref_fasta, rare_output)

    # Liftover unfiltered VCF
    unfiltered_output = os.path.splitext(args.unfiltered_vcf)[0] + '_hg38.vcf'
    liftover_vcf(args.unfiltered_vcf, args.chain_file, args.ref_fasta, unfiltered_output)

    print("Use the _hg38.vcf files as input for subsequent scripts.")
