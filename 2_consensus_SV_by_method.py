import os
import argparse
from intervaltree import Interval, IntervalTree
import vcf

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
        'record': record
    }

def build_interval_trees(sv_list, flank=100):
    trees = {}
    for sv in sv_list:
        key = (sv['chrom'], sv['svtype'])
        if key not in trees:
            trees[key] = IntervalTree()
        trees[key][sv['pos'] - flank : sv['pos'] + flank] = sv
    return trees

def generate_consensus_sv_fast(sniffles_vcf_path, cutesv_vcf_path, output_vcf_path, sample_name, max_end_diff=100):
    sniffles_reader = vcf.Reader(filename=sniffles_vcf_path)
    cutesv_reader = vcf.Reader(filename=cutesv_vcf_path)

    consensus_dir = f"{args.output_vcf}/consensus_vcf"
    os.makedirs(consensus_dir, exist_ok=True)

    # Prepare writer based on Sniffles2 header
    vcf_writer = vcf.Writer(open(f"{consensus_dir}/{sample_name}.method_consensus.vcf", 'w'), sniffles_reader)

    sniffles_svs = [parse_sv(r) for r in sniffles_reader]
    cutesv_svs = [parse_sv(r) for r in cutesv_reader]
    cutesv_tree = build_interval_trees(cutesv_svs)

    matched_cutesv_ids = set()

    for sv1 in sniffles_svs:
        key = (sv1['chrom'], sv1['svtype'])
        
        start = sv1['pos'] - max_end_diff
        end = sv1['pos'] + max_end_diff
        
        hits = cutesv_tree.get(key, IntervalTree()).overlap(start, end)

        best_match = None
        for hit in hits:
            sv2 = hit.data
            if abs(sv1['end'] - sv2['end']) <= max_end_diff:
                best_match = sv2
                break

        if best_match:
            matched_cutesv_ids.add(best_match['id'])

            # Build a new consensus record
            consensus_record = sv1['record']
            consensus_record.POS = int((sv1['pos'] + best_match['pos']) / 2)
            consensus_record.INFO['SVTYPE'] = sv1['svtype']
            consensus_record.INFO['END'] = int((sv1['end'] + best_match['end']) / 2)
            consensus_record.INFO['SVLEN'] = int((abs(sv1['svlen']) + abs(best_match['svlen'])) / 2)
            consensus_record.INFO['SOURCE'] = sample_name
            consensus_record.ID = f"consensus.{sv1['svtype']}.{sv1['chrom']}.{consensus_record.POS}"

            # Optional: adjust FORMAT field if needed
            consensus_record.samples = []

            # Write to VCF
            vcf_writer.write_record(consensus_record)

    vcf_writer.close()
    print(f"Consensus SV VCF {sample_name} written to: {output_vcf_path}/{sample_name}_method_consensus.vcf")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate consensus SV calls from Sniffles2 and cuteSV.")
    parser.add_argument("--sniffles_vcf", required=True, help="VCF file from Sniffles2")
    parser.add_argument("--cutesv_vcf", required=True, help="VCF file from cuteSV")
    parser.add_argument("--output_vcf", required=True, help="Path to output consensus VCF")
    parser.add_argument("--sample_name", required=True, help="Sample name to annotate SOURCE in INFO field")
    parser.add_argument("--max_end_diff", type=int, default=100, help="Maximum position difference to consider SVs matching")

    args = parser.parse_args()

    generate_consensus_sv_fast(
        sniffles_vcf_path=args.sniffles_vcf,
        cutesv_vcf_path=args.cutesv_vcf,
        output_vcf_path=args.output_vcf,
        sample_name=args.sample_name,
        max_end_diff=args.max_end_diff
    )
