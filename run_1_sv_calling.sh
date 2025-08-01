python 1_sv_calling.py \
  --input_dir /u/scratch/r/rmodlin/wqe/data \
  --output_dir /u/scratch/g/guyiqian/lung/out \
  --reference /u/scratch/g/guyiqian/ref/refdata-gex-GRCh38-2020-A/fasta/genome.fa \
  --data ont \
  --sv_callers sniffles2 cuteSV \
  --sniffles_dir /u/home/g/guyiqian/miniconda3/envs/wqe3/bin \
  --cuteSV_dir /u/home/g/guyiqian/miniconda3/envs/wqe3/bin \
  -n 8
