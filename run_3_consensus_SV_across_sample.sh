#!/bin/bash

# Exit on error
set -e

# Directory containing the VCF files
VCF_DIR="/u/scratch/g/guyiqian/lung/out/consensus_vcf/"

# Check that directory is not empty and exists
if [ -z "$VCF_DIR" ] || [ ! -d "$VCF_DIR" ]; then
  echo "Usage: $0 /path/to/vcf_dir/"
  echo "Error: Directory not found: $VCF_DIR"
  exit 1
fi

# Extract sample names from method_consensus.vcf files
sample_list=($(find "$VCF_DIR" -name "*.method_consensus.vcf" -exec basename {} .method_consensus.vcf  \; | sort))

echo "Samples with method consensus VCFs: ${sample_list[@]}"

# Loop through matching samples and run consensus
echo "Generating consensus for $sample..."
python 3_consensus_SV_across_sample.py \
--input_dir ${VCF_DIR} \
--output_tsv ${VCF_DIR}/sample_consensus.tsv \
--min_support 2 \
--max_distance 100
