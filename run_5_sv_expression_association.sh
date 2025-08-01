nput directory containing SV hotspot CSV files
SV_DIR="/u/scratch/g/guyiqian/lung/out/consensus_vcf/hotspot"

# Expression matrices
EXPR="/u/scratch/g/guyiqian/lung/rna/combined_exp_cancer_control.csv"

# Comma-separated sample names
CANCER_LIST="lung1,lung2,lung3,lung4,lung5" 
CONTROL_LIST="HG00597,HG00558,HG00544,HG00609,HG00673" 

# Annotation file and settings
GTF="/u/scratch/g/guyiqian/ref/refdata-gex-GRCh38-2020-A/genes/genes.gtf"
GTF_MODE="gene"
WINDOW_SIZE=250000
N_JOBS=8

# Output directory
OUT_DIR="/u/scratch/g/guyiqian/lung/out/consensus_vcf/association_2"
mkdir -p "$OUT_DIR"

# Loop through each SV hotspot file
for SV_FILE in "$SV_DIR"/*_hotspot_filtered.csv; do
    # Get filename without path and extension
    BASENAME=$(basename "$SV_FILE" .csv)

    # Construct output filename
    OUT_FILE="$OUT_DIR/${BASENAME}_gene_association.csv"

    # Run the new association script
    echo "Processing $SV_FILE..."
    python /u/scratch/g/guyiqian/lung/scripts/5_sv_expression_association_2.py \
        --expr "$EXPR" \
        --cancer_list "$CANCER_LIST" \
        --control_list "$CONTROL_LIST" \
        --sv "$SV_FILE" \
        --gtf "$GTF" \
        --gtf_mode "$GTF_MODE" \
        --window "$WINDOW_SIZE" \
        --n_jobs "$N_JOBS" \
        --out "$OUT_FILE"
done

echo "All jobs finished."

