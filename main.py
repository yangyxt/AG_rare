import argparse
import os
import subprocess

def run_command(cmd):
    subprocess.run(cmd, check=True)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Main pipeline for AlphaGenome variant analysis.')
    parser.add_argument('assembly', choices=['hg19', 'GRCh37', 'hg38', 'GRCh38'], help='VCF assembly (liftover if hg19/GRCh37)')
    parser.add_argument('rare_vcf', help='Path to rare variants VCF')
    parser.add_argument('unfiltered_vcf', help='Path to unfiltered variants VCF')
    parser.add_argument('ontology_terms', nargs='+', help='Ontology terms (e.g., UBERON:0001157)')
    parser.add_argument('--api_key', required=True, help='AlphaGenome API key')
    parser.add_argument('--modalities', nargs='+', default=['RNA_SEQ'], help='Modalities (e.g., RNA_SEQ)')
    parser.add_argument('--af_field', default='AF', help='AF INFO field')
    parser.add_argument('--threshold', type=float, default=0.1, help='Effect score threshold')
    parser.add_argument('--chain_file', default='hg19ToHg38.over.chain.gz', help='Liftover chain file')
    parser.add_argument('--ref_fasta', default='hg38.fa', help='hg38 reference FASTA')
    args = parser.parse_args()

    rare_vcf = args.rare_vcf
    unfiltered_vcf = args.unfiltered_vcf

    # Step 1: Liftover if needed
    if args.assembly in ['hg19', 'GRCh37']:
        print("Performing liftover to hg38...")
        rare_output = os.path.splitext(rare_vcf)[0] + '_hg38.vcf'
        unfiltered_output = os.path.splitext(unfiltered_vcf)[0] + '_hg38.vcf'
        run_command(['python', 'liftover.py', rare_vcf, unfiltered_vcf, '--chain_file', args.chain_file, '--ref_fasta', args.ref_fasta])
        rare_vcf = rare_output
        unfiltered_vcf = unfiltered_output

    # Step 2: Process VCFs to generate CSVs
    print("Processing VCFs...")
    run_command(['python', 'vcf_process.py', rare_vcf, unfiltered_vcf, '--af_field', args.af_field])

    # Step 3: Call AlphaGenome API
    print("Processing with AlphaGenome API...")
    ontology_str = ' '.join(['--ontology_terms'] + args.ontology_terms)
    outputs_str = ' '.join(['--requested_outputs'] + args.modalities)
    run_command(['python', 'AG_connect.py', 'variants_with_intervals.csv', '--api_key', args.api_key] + ontology_str.split() + outputs_str.split())

    # Step 4: Post-process and visualize
    print("Visualizing results...")
    modalities_str = ' '.join(['--modalities'] + args.modalities)
    run_command(['python', 'post_process.py', 'results', modalities_str, '--threshold', str(args.threshold)])

    print("Pipeline complete. Check visualizations/ for plots and scores.")
