import argparse
import os
import pandas as pd

# Imports from modular scripts
from crossmap_liftover import liftover_vcf
from load_interval_variants import load_rare_variants, cluster_rare_variants, generate_intervals_and_variants
from AG_connect import setup_api_client, process_variants, save_results
from vis_AG_results import load_results, filter_and_visualize

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
        liftover_vcf(rare_vcf, args.chain_file, args.ref_fasta, rare_output)
        unfiltered_output = os.path.splitext(unfiltered_vcf)[0] + '_hg38.vcf'
        liftover_vcf(unfiltered_vcf, args.chain_file, args.ref_fasta, unfiltered_output)
        rare_vcf = rare_output
        unfiltered_vcf = unfiltered_output

    # Step 2: Process VCFs to generate CSVs
    print("Processing VCFs...")
    rare_df = load_rare_variants(rare_vcf, args.af_field)
    if rare_df.empty:
        print("No rare variants found with the specified AF field.")
    else:
        clusters = cluster_rare_variants(rare_df)
        variants_df, intervals_df = generate_intervals_and_variants(clusters, unfiltered_vcf)
        variants_df.to_csv('variants_with_intervals.csv', index=False)
        intervals_df.to_csv('intervals.csv', index=False)

    # Step 3: Call AlphaGenome API
    print("Processing with AlphaGenome API...")
    variants_df = pd.read_csv('variants_with_intervals.csv')  # Reload if needed
    model = setup_api_client(args.api_key)
    results = process_variants(variants_df, model, args.ontology_terms, args.modalities)
    save_results(results, 'results')

    # Step 4: Post-process and visualize
    print("Visualizing results...")
    results = load_results('results')
    filter_and_visualize(results, args.modalities, args.threshold, 'visualizations')

    print("Pipeline complete. Check visualizations/ for plots and scores.")
