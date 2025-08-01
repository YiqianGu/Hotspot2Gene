#!/bin/bash

# Exit on error
set -e

# Directory containing the VCF files
VCF_DIR="/u/scratch/g/guyiqian/lung/out"

# Check that directory is not empty and exists
if [ -z "$VCF_DIR" ] || [ ! -d "$VCF_DIR" ]; then
  echo "Usage: $0 /path/to/vcf_dir/"
  echo "Error: Directory not found: $VCF_DIR"
  exit 1
fi

# Extract sample names from cuteSV and sniffles2 VCF files
list1=($(find "$VCF_DIR" -name "*.cuteSV.vcf" -exec basename {} .cuteSV.vcf \; | sort))
list2=($(find "$VCF_DIR" -name "*.sniffles2.vcf" -exec basename {} .sniffles2.vcf \; | sort))

# Find intersection: sample names present in both lists
sample_list=($(comm -12 <(printf "%s\n" "${list1[@]}") <(printf "%s\n" "${list2[@]}")))

echo "Samples with both VCFs: ${sample_list[@]}"

# Loop through matching samples and run consensus
for sample in "${sample_list[@]}"; do
  echo "Generating consensus for $sample..."

  python /u/scratch/g/guyiqian/lung/scripts/2_consensus_SV_by_method.py \
    --sniffles_vcf "${VCF_DIR}/${sample}.sniffles2.vcf" \
    --cutesv_vcf "${VCF_DIR}/${sample}.cuteSV.vcf" \
    --output_vcf "${VCF_DIR}" \
    --sample_name "$sample"
done

