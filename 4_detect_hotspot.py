#!/usr/bin/env python3

import numpy as np
import pandas as pd
from scipy.stats import norm
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from tqdm import tqdm
import argparse
import os
from numpy import unique

# ---------- Core Functions ----------

def smooth_sv_density(sv_data, chrom, bandwidth=500, window_size=1000, step_size=100, fixed_range=None):
    breakpoints = []
    svtype_contrib = []

    for sv in sv_data:
        if sv['chrom'] != chrom:
            continue
        svtype = sv['svtype']
        pos = sv['pos']
        end = sv['end']
        if svtype == 'INS':
            breakpoints.append(pos)
            svtype_contrib.append((pos, svtype))
        elif svtype in {'BND', 'INV'}:
            breakpoints.extend([pos, end])
            svtype_contrib.extend([(pos, svtype), (end, svtype)])
        elif svtype in {'DUP', 'DEL'}:
            region_points = list(range(pos, end + 1, step_size))
            breakpoints.extend(region_points)
            svtype_contrib.extend([(p, svtype) for p in region_points])

    if not breakpoints:
        return np.array([]), np.array([]), [], []

    if fixed_range is None:
        min_pos = max(0, min(breakpoints) - 5_000)
        max_pos = max(breakpoints) + 5_000
    else:
        min_pos, max_pos = fixed_range

    num_windows = int((max_pos - min_pos) // step_size)
    positions = np.linspace(min_pos, min_pos + step_size * num_windows, num_windows)
    densities = np.zeros_like(positions, dtype=float)
    svtype_matrix = [defaultdict(float) for _ in range(len(positions))]

    for bp, svtype in svtype_contrib:
        weights = norm.pdf((positions - bp), scale=bandwidth)
        densities += weights
        for i, w in enumerate(weights):
            svtype_matrix[i][svtype] += w

    svtype_summary = [max(counts.items(), key=lambda x: x[1])[0] if counts else 'NA' for counts in svtype_matrix]
    chrom_list = [chrom] * len(positions)
    return positions, densities, svtype_summary, chrom_list

def generate_shuffled_sv(chrom_sv, chrom, min_pos, max_pos):
    shuffled_sv = []
    for sv in chrom_sv:
        svtype = sv['svtype']
        length = sv['end'] - sv['pos']
        new_pos = np.random.randint(min_pos, max_pos - max(1, length))
        new_end = new_pos if svtype in ['INS', 'BND'] else new_pos + length
        shuffled_sv.append({
            'chrom': chrom,
            'pos': new_pos,
            'end': new_end,
            'svtype': svtype
        })
    return shuffled_sv

def shuffle_and_smooth(chrom_sv, chrom, bandwidth, window_size, step_size, fixed_range, seed):
    np.random.seed(seed)
    shuffled_sv = generate_shuffled_sv(chrom_sv, chrom, *fixed_range)
    _, density, _, _ = smooth_sv_density(shuffled_sv, chrom, bandwidth, window_size, step_size, fixed_range=fixed_range)
    return density

def detect_hotspots_by_shuffling_parallel(sv_data, chrom, bandwidth=500, window_size=1000,
                                          step_size=100, n_shuffle=1000, seed=42, n_jobs=1):
    np.random.seed(seed)
    positions, observed_density, svtype_summary, chrom_list = smooth_sv_density(
        sv_data, chrom, bandwidth, window_size, step_size)

    chrom_sv = [sv for sv in sv_data if sv['chrom'] == chrom]
    fixed_range = (min(positions), max(positions))

    seeds = np.random.randint(0, 1_000_000, size=n_shuffle)
    null_densities = []

    with ProcessPoolExecutor(max_workers=n_jobs) as executor:
        futures = [
            executor.submit(shuffle_and_smooth, chrom_sv, chrom, bandwidth, window_size, step_size, fixed_range, s)
            for s in seeds
        ]
        for f in tqdm(futures, desc="Shuffling in parallel"):
            null_densities.append(f.result())

    null_densities = np.array(null_densities)
    mean_null = null_densities.mean(axis=0)
    std_null = null_densities.std(axis=0) + 1e-6
    z_scores = (observed_density - mean_null) / std_null
    p_values = np.mean(null_densities >= observed_density, axis=0)

    results = pd.DataFrame({
        'chrom': chrom_list,
        'position': positions,
        'density': observed_density,
        'z_score': z_scores,
        'p_value': p_values,
        'svtype': svtype_summary
    })

    return results

# ---------- CLI ----------

def load_sv_file(path):
    df = pd.read_csv(path, sep='\t')
    return df.to_dict(orient='records')

def save_results(results, out_file):
    results.to_csv(out_file, index=False)

def main():
    parser = argparse.ArgumentParser(description="Detect SV hotspots using Gaussian smoothing + shuffling.")
    parser.add_argument('-i', '--input', required=True, help="Path to SV input TSV file.")
    parser.add_argument('-o', '--output', required=True, help="Path to output directory.")
    parser.add_argument('--bandwidth', type=int, default=500, help="Gaussian bandwidth in bp (default: 500)")
    parser.add_argument('--window-size', type=int, default=1000, help="Window size (default: 1000)")
    parser.add_argument('--step-size', type=int, default=100, help="Step size (default: 100)")
    parser.add_argument('--n-shuffle', type=int, default=1000, help="Number of shuffles (default: 1000)")
    parser.add_argument('--n-jobs', type=int, default=4, help="Number of parallel jobs (default: 4)")
    parser.add_argument('--seed', type=int, default=42, help="Random seed (default: 42)")

    args = parser.parse_args()

    sv_data = load_sv_file(args.input)
    chrom_list = list(set(row["chrom"] for row in sv_data))
    print(chrom_list)

    os.makedirs(args.output, exist_ok=True)

    for chrom in chrom_list:
        results = detect_hotspots_by_shuffling_parallel(
            sv_data, chrom=chrom,
            bandwidth=args.bandwidth,
            window_size=args.window_size,
            step_size=args.step_size,
            n_shuffle=args.n_shuffle,
            seed=args.seed,
            n_jobs=args.n_jobs
        )
        output_path = os.path.join(args.output, f"{chrom}_hotspot.csv")
        save_results(results, output_path)
        print(f"Done. Results written to {output_path}")

if __name__ == '__main__':
    main()
