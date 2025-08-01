python 4_detect_hotspot_v2.py \
    -i /u/scratch/g/guyiqian/lung/out/consensus_vcf/sample_consensus.tsv \
    -o /u/scratch/g/guyiqian/lung/out/consensus_vcf/hotspot/ \
    --bandwidth 500 \
    --window-size 1000 \
    --step-size 1000 \
    --n-shuffle 100 \
    --n-jobs 8

