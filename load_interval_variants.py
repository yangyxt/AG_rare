import argparse
import pandas as pd
import pysam


def load_rare_variants(rare_vcf_path, af_field='AF'):
    rare_vcf = pysam.VariantFile(rare_vcf_path)
    rare_records = []
    for rec in rare_vcf:
        af = rec.info.get(af_field)
        if af is None:
            continue  # Skip if no AF field
        af = af[0] if isinstance(af, tuple) else af
        alts = ','.join(rec.alts)  # Handle multi-allelic
        rare_records.append({'chr': rec.chrom, 'pos': rec.pos, 'ref': rec.ref, 'alt': alts, 'af': af})
    rare_df = pd.DataFrame(rare_records)
    return rare_df


def cluster_rare_variants(rare_df):
    clusters = []
    for chr_, chr_group in rare_df.groupby('chr'):
        if chr_group.empty:
            continue
        chr_group = chr_group.sort_values('pos')
        current_vars = [chr_group.iloc[0]]
        current_min = chr_group.iloc[0]['pos']
        for i in range(1, len(chr_group)):
            next_pos = chr_group.iloc[i]['pos']
            if next_pos - current_min <= 1000000:
                current_vars.append(chr_group.iloc[i])
            else:
                clusters.append(pd.DataFrame(current_vars))
                current_min = next_pos
                current_vars = [chr_group.iloc[i]]
        if current_vars:
            clusters.append(pd.DataFrame(current_vars))
    return clusters


def generate_intervals_and_variants(clusters, unfiltered_vcf_path):
    unfiltered_vcf = pysam.VariantFile(unfiltered_vcf_path)
    variant_data = []
    interval_set = set()
    for cluster in clusters:
        if len(cluster) == 0:
            continue
        min_pos = cluster['pos'].min()
        max_pos = cluster['pos'].max()
        chr_ = cluster['chr'].iloc[0]
        span = max_pos - min_pos
        if span > 1000000:
            print(f"Warning: Cluster on {chr_} spans {span} bp > 1M bp, but proceeding as is (may need splitting).")

        rarest_idx = cluster['af'].idxmin()
        rarest_pos = cluster.loc[rarest_idx, 'pos']

        desired_start = rarest_pos - 500000
        start_low = max_pos - 1000000
        start_high = min_pos
        start = max(start_low, min(start_high, desired_start))
        end = start + 1000000

        interval_str = f"{chr_}:{int(start)}-{int(end)}"
        interval_set.add(interval_str)

        for rec in unfiltered_vcf.fetch(chr_, int(start) - 1, int(end) + 1):
            if rec.pos >= start and rec.pos <= end:
                alts = ','.join(rec.alts)
                variant_data.append({
                    'chr': rec.chrom,
                    'pos': rec.pos,
                    'ref': rec.ref,
                    'alt': alts,
                    'interval': interval_str
                })
    return pd.DataFrame(variant_data), pd.DataFrame({'interval': list(interval_set)})


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process VCF files to generate intervals and variants for AlphaGenome.')
    parser.add_argument('rare_vcf', help='Path to the rare variants VCF file')
    parser.add_argument('unfiltered_vcf', help='Path to the unfiltered variants VCF file')
    parser.add_argument('--af_field', default='AF', help='INFO field name for population allele frequency (default: AF)')
    args = parser.parse_args()

    rare_df = load_rare_variants(args.rare_vcf, args.af_field)
    if rare_df.empty:
        print("No rare variants found with the specified AF field.")
    else:
        clusters = cluster_rare_variants(rare_df)
        variants_df, intervals_df = generate_intervals_and_variants(clusters, args.unfiltered_vcf)

        variants_df.to_csv('variants_with_intervals.csv', index=False)
        intervals_df.to_csv('intervals.csv', index=False)

        print("Generated files: variants_with_intervals.csv and intervals.csv")
