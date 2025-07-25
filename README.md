# AG_rare: AlphaGenome Pipeline for Rare Variant Analysis

[![GitHub license](https://img.shields.io/github/license/yangyxt/AG_rare)](https://github.com/yangyxt/AG_rare/blob/main/LICENSE)
[![Python Version](https://img.shields.io/badge/python-3.11%2B-blue)](https://www.python.org/)

## Overview

This repository provides a modular Python pipeline for analyzing rare variants from patient VCF files using Google DeepMind's AlphaGenome API. The pipeline:
- Handles liftover from GRCh37/hg19 to GRCh38/hg38 if needed.
- Clusters rare variants (<0.1% PAF) into 1M bp intervals.
- Extracts all variants within those intervals from an unfiltered VCF.
- Queries the AlphaGenome API for variant effect predictions on regulatory elements (e.g., tissue-specific disruptions).
- Computes effect scores, filters significant disruptions, and generates visualizations.

AlphaGenome predicts multimodal genomic features (e.g., RNA expression, chromatin accessibility) and scores variant impacts, focusing on non-coding regions like exon-adjacent areas. This is suitable for small- to medium-scale analyses due to API quotas.

**Note:** AlphaGenome is free for non-commercial use but requires an API key from https://deepmind.google/technologies/alphagenome/. Respect usage limits.

## Prerequisites

- Python 3.11+ (via Miniconda recommended for CentOS 7 compatibility).
- Install dependencies:
  ```
  pip install pandas pysam alphagenome CrossMap matplotlib numpy
  ```
- For liftover: Download hg19ToHg38 chain file (`wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz`) and hg38 FASTA (`wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz` and unzip).
- AlphaGenome API key (register at https://deepmind.google/technologies/alphagenome/).

## Installation

Clone the repo:
```
git clone https://github.com/yangyxt/AG_rare.git
cd AG_rare
```

Install packages as above.

## Usage

Run the main pipeline script with required arguments:

```
python main.py <assembly> <rare_vcf> <unfiltered_vcf> <ontology_terms> --api_key YOUR_API_KEY [optional args]
```

### Arguments
- `assembly`: VCF genome build (hg19, GRCh37, hg38, or GRCh38). Liftover applied if hg19/GRCh37.
- `rare_vcf`: Path to VCF with rare variants (<0.1% PAF, filtered externally).
- `unfiltered_vcf`: Path to full (unfiltered) VCF from the same patient.
- `ontology_terms`: Space-separated list of ontology terms (e.g., UBERON:0001157 for colon).
- `--api_key`: Required AlphaGenome API key.
- `--modalities` (optional): Space-separated modalities (default: RNA_SEQ; others: CHROMATIN_ACCESSIBILITY, etc.).
- `--af_field` (optional): VCF INFO field for allele frequency (default: AF).
- `--threshold` (optional): Effect score filter for visualization (default: 0.1).
- `--chain_file` (optional): Liftover chain file path (default: hg19ToHg38.over.chain.gz).
- `--ref_fasta` (optional): hg38 FASTA path (default: hg38.fa).

### Example
```
python main.py hg19 rare.vcf unfiltered.vcf UBERON:0001157 UBERON:0001114 --api_key abc123 --modalities RNA_SEQ CHROMATIN_ACCESSIBILITY --threshold 0.05
```

### Output
- `variants_with_intervals.csv` & `intervals.csv`: Intermediate files with variants and intervals.
- `results/`: Pickled AlphaGenome API outputs.
- `visualizations/`: PNG plots of REF vs. ALT tracks for filtered variants, plus `effect_scores.csv`.

## Scripts Overview

- `main.py`: Orchestrates the full pipeline via function imports.
- `crossmap_liftover.py`: Lifts over VCFs from hg19/GRCh37 to hg38 using CrossMap.
- `load_interval_variants.py`: Parses VCFs, clusters rare variants into 1M bp intervals, extracts all variants within them.
- `AG_connect.py`: Calls AlphaGenome API for variant predictions, saves results.
- `vis_AG_results.py`: Loads results, computes mean absolute differences as effect scores, filters, and visualizes tracks using AlphaGenome's plotting tools.

## Ontology Terms

Ontology terms specify biological contexts (e.g., tissues, cell types) for AlphaGenome predictions, drawn from UBERON (anatomy), CL (cell types), CLO (cell lines), and EFO (experimental factors). Only terms in AlphaGenome's training data (e.g., ENCODE biosamples) are supported; invalid terms raise errors.

For the full list of supported terms, run this official Colab notebook and export the DataFrame: https://colab.research.google.com/github/google-deepmind/alphagenome/blob/main/colabs/tissue_ontology_mapping.ipynb.

Examples:
- UBERON:0001157 (colon)
- UBERON:0001114 (liver)
- CL:0000127 (astrocyte)

Validate terms via EBI OLS: https://www.ebi.ac.uk/ols4/.

## Limitations
- Designed for rare variant analysis; not optimized for whole-genome scale (filter variants first).
- API quotas apply—test on small datasets.
- Assumes VCFs use "chr" prefixes and are from the same patient.

## Contributing
Pull requests welcome. For issues, open a GitHub issue.

## License
MIT License (see [LICENSE](LICENSE)).
