# OMHI: Oral Microbiome Health Index

<p align="center">
  <img src="images/omhi_workflow.png" width="850">
</p>

<p align="center">
  <img src="https://img.shields.io/badge/license-MIT-green">
  <img src="https://img.shields.io/badge/platform-Linux%20%7C%20HPC-blue">
  <img src="https://img.shields.io/badge/input-paired--end%20FASTQ-orange">
  <img src="https://img.shields.io/badge/profiler-MetaPhlAn4-purple">
  <img src="https://img.shields.io/badge/status-beta-yellow">
</p>

## Description

OMHI, the **Oral Microbiome Health Index**, is a metagenome-based score designed to quantify whether an oral microbial profile is closer to a health-associated or disease-associated oral community state.

This repository provides a full OMHI inference pipeline starting from paired-end shotgun metagenomic FASTQ files, as well as a lightweight scoring script for users who already have MetaPhlAn4 species-level abundance profiles.

The full pipeline performs:

```text
FASTQ
→ quality control
→ host read removal
→ MetaPhlAn4 taxonomic profiling
→ species-level abundance extraction
→ alignment to the fixed OMHI species universe
→ CLR transformation
→ OMHI score calculation

Installation

Clone the repository:

git clone https://github.com/minjunhoo2001/OMHI.git
cd OMHI

Create the conda environment:

conda env create -f environment.yml
conda activate omhi_env

Check installation:

bash workflow/test_installation.sh

If all required commands and model files are detected, the installation is ready.

Required databases

The full OMHI pipeline requires two external databases:

A human Bowtie2 index for host read removal.
A MetaPhlAn4 database for taxonomic profiling.

These databases are not included in this repository.

Example paths:

HOST_INDEX=/path/to/human_bowtie2_index/human_ref
MPA_DB=/path/to/metaphlan4_database
MPA_INDEX=mpa_vJun23_CHOCOPhlAnSGB_202403

HOST_INDEX should be the Bowtie2 index prefix, not the folder name.

For example, if the index files are:

/path/to/human_ref.1.bt2
/path/to/human_ref.2.bt2
/path/to/human_ref.3.bt2
/path/to/human_ref.4.bt2

then use:

HOST_INDEX=/path/to/human_ref
Input format

The full pipeline requires a tab-delimited sample sheet:

sample_id	r1	r2
sample01	/path/to/sample01_R1.fastq.gz	/path/to/sample01_R2.fastq.gz
sample02	/path/to/sample02_R1.fastq.gz	/path/to/sample02_R2.fastq.gz

The header must be exactly:

sample_id	r1	r2
Run the full OMHI pipeline
THREADS=8 \
HOST_INDEX=/path/to/human_ref \
MPA_DB=/path/to/metaphlan4_database \
MPA_INDEX=mpa_vJun23_CHOCOPhlAnSGB_202403 \
bash workflow/run_omhi_pipeline.sh samples.tsv results

The final OMHI scores will be written to:

results/05.omhi/omhi_scores.tsv
Output

The output file contains:

sample_id	OMHI	predicted_state	matched_universe_taxa	total_universe_taxa
Column	Description
sample_id	Sample identifier
OMHI	Oral Microbiome Health Index score
predicted_state	health_associated if OMHI > 0; otherwise disease_associated
matched_universe_taxa	Number of OMHI universe species matched in the input
total_universe_taxa	Number of species in the fixed OMHI species universe

Interpretation:

OMHI > 0    more health-associated oral microbial profile
OMHI < 0    more disease-associated oral microbial profile
Calculate OMHI from an existing species-level table

If MetaPhlAn4 has already been run, users can calculate OMHI directly from a species-level abundance table.

Input format:

feature	sample01	sample02
Neisseria_subflava	0.012	0.004
Streptococcus_sanguinis	0.003	0.001
Tannerella_forsythia	0.000	0.005

Run:

Rscript scripts/calc_omhi.R \
  --input taxonomy_species.tsv \
  --coef model/omhi_species_coefficients.tsv \
  --universe model/omhi_species_universe.tsv \
  --output omhi_scores.tsv

Important note:

OMHI is calculated after aligning MetaPhlAn4 species profiles to the fixed OMHI species universe.
Species not detected in a sample are assigned zero abundance before CLR transformation.
CLR transformation is always performed over this fixed OMHI species universe.
Calculate OMHI from a MetaPhlAn4 merged table

If users already have a merged MetaPhlAn4 table, first extract species-level abundance:

python scripts/metaphlan_level_extractor.py \
  metaphlan4_merged.tsv \
  level

Then calculate OMHI:

Rscript scripts/calc_omhi.R \
  --input level/species.tsv \
  --coef model/omhi_species_coefficients.tsv \
  --universe model/omhi_species_universe.tsv \
  --output omhi_scores.tsv
HPC / LSF example

An example LSF submission script is provided:

bsub < submit/submit_omhi_test.lsf

Users should modify the following variables in the submission script:

PROJECT_DIR="/path/to/OMHI"
HOST_INDEX="/path/to/human_ref"
MPA_DB="/path/to/metaphlan4_database"
MPA_INDEX="mpa_vJun23_CHOCOPhlAnSGB_202403"
SAMPLE_SHEET="samples.tsv"
OUTDIR="results"
Model files

The released OMHI model files are stored in model/:

model/
├── omhi_species_coefficients.tsv
├── omhi_species_coefficients_full.tsv
├── omhi_species_coefficients_nonzero.tsv
├── omhi_species_universe.tsv
└── omhi_species_model_config.tsv

omhi_species_universe.tsv defines the fixed species set used for CLR transformation.

omhi_species_coefficients.tsv contains the coefficients used to calculate OMHI.

Citation

If you use OMHI, please cite:

[Your OMHI manuscript citation here]
License

This project is released under the MIT License.