import os
import glob
import argparse

def run_sv_callers(input_dir, output_dir, reference, data, sv_callers = ['sniffles2', 'cuteSV'], n=1, sniffles_dir=None, cuteSV_dir=None):
    # input_dir: the path to bam files
    # output_dir: the path to store output files
    # reference: the reference genome
    # sv_callers: a list of sv callers to use. defalut: sniffles2 and cuteSV
    # data: choose from 'pacbio_hifi' and 'ont'
    # n: number of cores to use
    # sniffles_dir: the path to sniffles2
    # cuteSV_dir: the path to cuteSV
    
    bam_files = glob.glob(os.path.join(input_dir, "*.bam"))
    print(f"Found BAM files: {bam_files}")

    for bam in bam_files:
        sample_name = os.path.basename(bam).removesuffix(".bam")
        for caller in sv_callers:
            if caller == 'sniffles2':
                output_vcf = os.path.join(output_dir, f"{sample_name}.sniffles2.vcf")
                cmd = f"{sniffles_dir}/sniffles -i {bam} -v {output_vcf} --reference {reference} --threads {n}"
            elif caller == 'cuteSV':
                if data == "pacbio_hifi":
                    output_vcf = os.path.join(output_dir, f"{sample_name}.cuteSV.vcf")
                    cmd = f"{cuteSV_dir}/cuteSV {bam} {reference} {output_vcf} {output_dir} --threads {n} --max_cluster_bias_INS 1000 --diff_ratio_merging_INS 0.9 --max_cluster_bias_DEL 1000 --diff_ratio_merging_DEL 0.5"
                elif data == "ont":
                    output_vcf = os.path.join(output_dir, f"{sample_name}.cuteSV.vcf")
                    cmd = f"{cuteSV_dir}/cuteSV {bam} {reference} {output_vcf} {output_dir} --threads {n} --max_cluster_bias_INS 1000 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 100 --diff_ratio_merging_DEL 0.3"
            else:
                print(f"Unsupported caller: {caller}")
                continue

            print(f"Running: {cmd}")
            os.system(cmd)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run Sniffles2 and/or cuteSV on BAM files.")
    parser.add_argument("--input_dir", required=True, help="Directory containing BAM files")
    parser.add_argument("--output_dir", required=True, help="Directory to store output VCFs")
    parser.add_argument("--reference", required=True, help="Reference genome file (FASTA)")
    parser.add_argument("--data", choices=["pacbio_hifi", "ont"], required=True, help="Sequencing data type")
    parser.add_argument("--sv_callers", nargs="+", default=["sniffles2", "cuteSV"], help="SV callers to use")
    parser.add_argument("-n", type=int, default=1, help="Number of threads to use")
    parser.add_argument("--sniffles_dir", help="Path to Sniffles2 binary directory")
    parser.add_argument("--cuteSV_dir", help="Path to cuteSV binary directory")

    args = parser.parse_args()

    run_sv_callers(
        input_dir=args.input_dir,
        output_dir=args.output_dir,
        reference=args.reference,
        data=args.data,
        sv_callers=args.sv_callers,
        n=args.n,
        sniffles_dir=args.sniffles_dir,
        cuteSV_dir=args.cuteSV_dir
    )


