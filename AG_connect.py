import argparse
import os
import pandas as pd
from alphagenome.data import genome
from alphagenome.models import dna_client
import pickle


def setup_api_client(api_key):
    return dna_client.create(api_key)


def process_variants(variants_df, model, ontology_terms, requested_outputs):
    results = {}
    for interval_str, group in variants_df.groupby('interval'):
        chr_, start_end = interval_str.split(':')
        start, end = map(int, start_end.split('-'))
        interval = genome.Interval(chromosome=chr_, start=start, end=end)

        for _, row in group.iterrows():
            variant = genome.Variant(
                chromosome=row['chr'],
                position=row['pos'],
                reference_bases=row['ref'],
                alternate_bases=row['alt']
            )
            key = f"{row['chr']}_{row['pos']}_{row['ref']}_{row['alt']}"
            try:
                outputs = model.predict_variant(
                    interval=interval,
                    variant=variant,
                    ontology_terms=ontology_terms,
                    requested_outputs=requested_outputs
                )
                results[key] = outputs
            except Exception as e:
                print(f"Error processing variant {key}: {e}")
    return results


def save_results(results, output_dir='results'):
    os.makedirs(output_dir, exist_ok=True)
    for key, outputs in results.items():
        with open(os.path.join(output_dir, f"{key}.pkl"), 'wb') as f:
            pickle.dump(outputs, f)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process variants with AlphaGenome API.')
    parser.add_argument('variants_csv', help='Path to variants_with_intervals.csv')
    parser.add_argument('--api_key', required=True, help='AlphaGenome API key')
    parser.add_argument('--ontology_terms', nargs='+', default=['UBERON:0001157'], help='Ontology terms (e.g., tissues)')
    parser.add_argument('--requested_outputs', nargs='+', default=['RNA_SEQ'], help='Output types (e.g., RNA_SEQ)')
    parser.add_argument('--output_dir', default='results', help='Directory to save pickled results')
    args = parser.parse_args()

    variants_df = pd.read_csv(args.variants_csv)
    model = setup_api_client(args.api_key)
    results = process_variants(variants_df, model, args.ontology_terms, args.requested_outputs)
    save_results(results, args.output_dir)

    print(f"Results saved to {args.output_dir}")
