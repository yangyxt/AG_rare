import argparse
import os
import pickle
import pandas as pd
import numpy as np
from alphagenome.visualization import plot_components
import matplotlib.pyplot as plt
from alphagenome.data import genome

def load_results(results_dir):
    results = {}
    for filename in os.listdir(results_dir):
        if filename.endswith('.pkl'):
            key = filename[:-4]
            with open(os.path.join(results_dir, filename), 'rb') as f:
                results[key] = pickle.load(f)
    return results

def compute_effect_scores(outputs, modalities):
    scores = {}
    for modality in modalities:
        ref_track = getattr(outputs.reference, modality.lower(), None)
        alt_track = getattr(outputs.alternate, modality.lower(), None)
        if ref_track is not None and alt_track is not None:
            diff = np.abs(ref_track - alt_track)
            scores[modality] = np.mean(diff)  # Mean abs diff as effect score
    return scores

def filter_and_visualize(results, modalities, threshold, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    score_data = []
    for key, outputs in results.items():
        scores = compute_effect_scores(outputs, modalities)
        max_score = max(scores.values()) if scores else 0
        if max_score > threshold:
            chr_, pos, ref, alt = key.split('_')
            variant = genome.Variant(chromosome=chr_, position=int(pos), reference_bases=ref, alternate_bases=alt)
            
            # Visualize (example: overlay tracks for each modality)
            for modality in modalities:
                ref_track = getattr(outputs.reference, modality.lower())
                alt_track = getattr(outputs.alternate, modality.lower())
                plot_components.plot(
                    [plot_components.OverlaidTracks(
                        tdata={'REF': ref_track, 'ALT': alt_track},
                        colors={'REF': 'dimgrey', 'ALT': 'red'}
                    )],
                    interval=ref_track.interval.resize(2**15),  # Adjust zoom as needed
                    annotations=[plot_components.VariantAnnotation([variant], alpha=0.8)]
                )
                plot_file = os.path.join(output_dir, f"{key}_{modality}.png")
                plt.savefig(plot_file)
                plt.close()
                print(f"Saved visualization: {plot_file}")
            
            # Collect scores
            score_data.append({'variant_key': key, **scores})
    
    # Save scores CSV
    scores_df = pd.DataFrame(score_data)
    scores_df.to_csv(os.path.join(output_dir, 'effect_scores.csv'), index=False)
    print(f"Saved effect scores: {os.path.join(output_dir, 'effect_scores.csv')}")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Post-process and visualize AlphaGenome results.')
    parser.add_argument('results_dir', help='Directory with pickled results from Script 2')
    parser.add_argument('--modalities', nargs='+', default=['RNA_SEQ'], help='Modalities to process (e.g., RNA_SEQ CHROMATIN_ACCESSIBILITY)')
    parser.add_argument('--threshold', type=float, default=0.1, help='Effect score threshold for filtering/visualization')
    parser.add_argument('--output_dir', default='visualizations', help='Directory for plots and scores CSV')
    args = parser.parse_args()

    results = load_results(args.results_dir)
    filter_and_visualize(results, args.modalities, args.threshold, args.output_dir)
