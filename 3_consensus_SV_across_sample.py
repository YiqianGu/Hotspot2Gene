#!/usr/bin/env python3

import os
import argparse
import vcf
from intervaltree import IntervalTree
from collections import defaultdict
from Levenshtein import distance as lev_distance

def parse_sv(record):
    return {
        'chrom': record.CHROM,
        'pos': record.POS,
        'end': int(record.INFO.get('END', record.POS)),
        'alt': record.ALT,
        'svtype': record.INFO['SVTYPE'],
        'svlen': int(record.INFO.get('SVLEN', 0)),
        'id': record.ID,
        'info': record.INFO,
        'sample': record.INFO.get('SOURCE', 'unknown')
    }

class CrossSampleConsensusSV:
    def __init__(self, max_distance=100, min_support=2):
        self.max_distance = max_distance
        self.min_support = min_support
        self.sv_trees = defaultdict(IntervalTree)
        self.cluster_id = 0
        self.similarity_threshold = 0.5

    def add_vcf_file(self, vcf_path):
        reader = vcf.Reader(filename=vcf_path)

        for record in reader:
            sv = parse_sv(record)
            if not sv:
                continue

            key = (sv['chrom'], sv['svtype'])

            if sv['svtype'] == 'INS':
                raw_matches = self.sv_trees[key].overlap(sv['pos'] - self.max_distance, sv['end'] + self.max_distance)
                matches = []
                for interval in raw_matches:
                    candidate_sv = interval.data
                    if sv['alt'] and candidate_sv['svs'][0]['alt']:
                        seq1 = str(sv['alt'][0])
                        seq2 = str(candidate_sv['svs'][0]['alt'][0])
                        similarity = 1 - lev_distance(seq1, seq2) / max(len(seq1), len(seq2))
                        if similarity >= self.similarity_threshold:
                            matches.append(interval)
            else:
                matches = self.sv_trees[key].overlap(sv['pos'] - self.max_distance, sv['end'] + self.max_distance)

            if matches:
                matched = sorted(matches, key=lambda m: abs(m.begin - sv['pos']))[0]
                if sv['sample'] not in matched.data['samples']:
                    matched.data['samples'].append(sv['sample'])
                    matched.data['svs'].append(sv)
            else:
                start = sv['pos']
                end = sv['end'] if sv['end'] > sv['pos'] else sv['pos'] + 1
                self.sv_trees[key][start:end] = {
                    "cluster_id": self.cluster_id,
                    "svs": [sv],
                    "samples": [sv['sample']]
                }
                self.cluster_id += 1

    def finalize_consensus(self):
        final_consensus = []
        for tree in self.sv_trees.values():
            for interval in tree:
                data = interval.data
                if len(data["samples"]) >= self.min_support:
                    representative = max(data["svs"], key=lambda sv: abs(sv['svlen']))
                    final_consensus.append({
                        "chrom": representative['chrom'],
                        "pos": representative['pos'],
                        "end": representative['end'],
                        "svtype": representative['svtype'],
                        "svlen": representative['svlen'],
                        "samples": data["samples"],
                        "support": len(data["samples"])
                    })
        return final_consensus

def write_consensus_to_tsv(consensus_list, output_file):
    with open(output_file, 'w') as f:
        f.write("chrom\tpos\tend\tsvtype\tsvlen\tsupport\tsamples\n")
        for sv in consensus_list:
            f.write(f"{sv['chrom']}\t{sv['pos']}\t{sv['end']}\t{sv['svtype']}\t{sv['svlen']}\t{sv['support']}\t{sv['samples']}\n")

def main():
    parser = argparse.ArgumentParser(description="Cross-sample SV consensus tool")
    parser.add_argument("--input_dir", required=True, help="Directory containing VCF files")
    parser.add_argument("--output_tsv", required=True, help="Output TSV file for consensus SVs")
    parser.add_argument("--min_support", type=int, default=2, help="Minimum sample support to report consensus")
    parser.add_argument("--max_distance", type=int, default=100, help="Max breakpoint distance for clustering")
    args = parser.parse_args()

    consensus_builder = CrossSampleConsensusSV(max_distance=args.max_distance, min_support=args.min_support)

    for filename in os.listdir(args.input_dir):
        if filename.endswith(".vcf"):
            vcf_path = os.path.join(args.input_dir, filename)
            print(f"Processing {vcf_path}")
            consensus_builder.add_vcf_file(vcf_path)

    consensus = consensus_builder.finalize_consensus()
    print(f"Total consensus SVs: {len(consensus)}")

    write_consensus_to_tsv(consensus, args.output_tsv)
    print(f"Written consensus SVs to: {args.output_tsv}")

if __name__ == "__main__":
    main()

