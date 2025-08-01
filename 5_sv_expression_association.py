import argparse
import pandas as pd
import numpy as np
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests
from joblib import Parallel, delayed
import gzip

# ------------------------- Utility Functions -------------------------

def counts_to_cpm(count_df):
    """
    Convert raw counts to CPM (Counts Per Million)
    """
    lib_sizes = count_df.sum(axis=0)
    cpm = count_df.divide(lib_sizes, axis=1) * 1e6
    return cpm


def read_gtf_gene_positions(gtf_path, mode='gene'):
    """
    Parse GTF to extract gene or transcript positions.
    Returns: dict {gene_name: {chrom, start, end, strand}}
    """
    gene_info = {}
    open_func = gzip.open if gtf_path.endswith('.gz') else open
    with open_func(gtf_path, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if fields[2] != mode:
                continue
            chrom, start, end, strand, attr = fields[0], int(fields[3]), int(fields[4]), fields[6], fields[8]
            
            gene_name = ''
            for item in attr.strip().split(';'):
                item = item.strip()
                if mode == 'gene' and 'gene_name' in item:
                    gene_name = item.split('"')[1]
                elif mode == 'transcript' and 'transcript_name' in item:
                    gene_name = item.split('"')[1]
            
            if gene_name and gene_name not in gene_info:
                gene_info[gene_name] = {
                    'chrom': chrom,
                    'start': start,
                    'end': end,
                    'strand': strand
                }
    return gene_info


def differential_expression_near_hotspot(expr_df, cancer_samples, control_samples,
                                        sv_file, gtf_info, window=500_000):
    """
    Compute differential expression of genes near SV hotspots.
    Returns a DataFrame with log2FC and adjusted p-values.
    """
    sv_df = pd.read_csv(sv_file)
    results = []

    for _, sv in sv_df.iterrows():
        chrom, pos, svtype = sv['chrom'], sv['position'], sv['svtype']

        # Find nearby genes
        for gene, info in gtf_info.items():
            if info['chrom'] != chrom:
                continue
            if abs(info['start'] - pos) > window and abs(info['end'] - pos) > window:
                continue

            if gene not in expr_df.index:
                continue

            # Extract expression
            c_vals = expr_df.loc[gene, cancer_samples].values
            n_vals = expr_df.loc[gene, control_samples].values

            # Compute t-test and log2FC
            t_stat, p_val = ttest_ind(c_vals, n_vals, equal_var=False)
            log2fc = np.log2(np.mean(c_vals) + 1) - np.log2(np.mean(n_vals) + 1)

            results.append([chrom, pos, svtype, gene, info['strand'], log2fc, p_val])

    # Multiple testing correction
    results_df = pd.DataFrame(results, columns=[
        'chrom', 'sv_pos', 'svtype', 'gene', 'strand', 'log2FC', 'p_value'
    ])
    if not results_df.empty:
        results_df['adj_p_value'] = multipletests(results_df['p_value'], method='fdr_bh')[1]
    
    return results_df


# ------------------------- Main -------------------------

def main():
    parser = argparse.ArgumentParser(description="Associate SV hotspots with differential gene expression")
    parser.add_argument('--expr', required=True, help='Raw count expression matrix CSV (rows=genes, columns=samples)')
    parser.add_argument('--cancer_list', required=True, help='Comma-separated list of cancer sample names')
    parser.add_argument('--control_list', required=True, help='Comma-separated list of control sample names')
    parser.add_argument('--gtf', required=True, help='GTF file')
    parser.add_argument('--gtf_mode', default='gene', choices=['gene', 'transcript'])
    parser.add_argument('--sv', required=True, help='SV hotspot CSV file')
    parser.add_argument('--window', type=int, default=500_000, help='Window size for nearby genes')
    parser.add_argument('--n_jobs', type=int, default=1)
    parser.add_argument('--out', required=True, help='Output CSV file path')

    args = parser.parse_args()

    # Load expression and normalize to CPM
    expr_df = pd.read_csv(args.expr, index_col=0)
    expr_df = counts_to_cpm(expr_df)

    # Parse sample groups
    cancer_samples = args.cancer_list.split(',')
    control_samples = args.control_list.split(',')

    # Parse GTF
    gtf_info = read_gtf_gene_positions(args.gtf, mode=args.gtf_mode)

    # Run analysis
    results_df = differential_expression_near_hotspot(
        expr_df, cancer_samples, control_samples,
        args.sv, gtf_info,
        window=args.window
    )

    # Save
    results_df.to_csv(args.out, index=False)


if __name__ == "__main__":
    main()

