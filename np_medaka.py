from pathlib import Path
import re
from argparse import ArgumentParser

from utils import Cmds
from configs import SOFTWARE, REF_FA


class Medaka:

    def __init__(self):
        ap = ArgumentParser()
        ap.add_argument("--in_dir", type=str, required=True, help="trimmed_bams/*.trimmed.fq")
        ap.add_argument("--out_dir", type=str, required=True, help="consensus_vcfs/*.consensus.vcf")
        ap.add_argument("--cpu", type=int, default=1)
        self.args = ap.parse_args()
        self.cpu = int(self.args.cpu)
        self.in_dir = Path(self.args.in_dir).absolute()
        self.out_dir = Path(self.args.out_dir).absolute()
        self.in_files = list(
            self.in_dir.glob("*.trimmed.fq")
        )
        self.samples = [
            re.sub(r"\.trimmed\.fq$", "", i.name)
            for i in self.in_files
        ]
        self.temp_dirs = [
            self.out_dir.joinpath(i)
            for i in self.samples
        ]
        self.out_consensus_fas = [
            i.joinpath("consensus.fasta")
            for i in self.temp_dirs
        ]
        self.out_hdfs = [
            i.joinpath("consensus_probs.hdf")
            for i in self.temp_dirs
        ]
        self.out_bams = [
            i.joinpath("calls_to_draft.bam")
            for i in self.temp_dirs
        ]
        self.raw_vcfs = [
            i.joinpath("raw.vcf")
            for i in self.temp_dirs
        ]
        self.annotated_vcfs = [
            i.joinpath("annotated.vcf")
            for i in self.temp_dirs
        ]
        self.longshot_vcfs = [
            i.joinpath("longshot.vcf")
            for i in self.temp_dirs
        ]
        self.consensus_vcfs = [
            self.out_dir.joinpath(f"{s}.consensus.vcf")
            for s in self.samples
        ]
        self.sh_dir = self.out_dir.joinpath("work_sh")
        self.sh_dir.mkdir(parents=True, exist_ok=True)
        self.shell_files = [
            self.sh_dir.joinpath(f"{s}.sh")
            for s in self.samples
        ]

    def medaka(self):
        Cmds(
            [
                f"export PATH={SOFTWARE.joinpath('bin')}:$PATH && "
                f"medaka_consensus -i {in_fq} -d {REF_FA} -t 8 -m r941_min_hac_variant_g507 -o {temp_dir} && "
                f"medaka variant --verbose {REF_FA} {out_hdf} {raw_vcf} && "
                f"medaka tools annotate {raw_vcf} {REF_FA} {out_bam} {annotated_vcf} && "
                f"longshot -P 0 -F -A --no_haps --bam {out_bam} --ref {REF_FA} "
                f"--out {longshot_vcf} --potential_variants {annotated_vcf} && "
                f"bcftools view -i \"FILTER=='PASS'\" {longshot_vcf} > {consensus_vcf}\n"
                for in_fq, temp_dir, out_hdf, raw_vcf, out_bam, annotated_vcf, longshot_vcf, consensus_vcf
                in zip(
                        self.in_files, self.temp_dirs, self.out_hdfs,
                        self.raw_vcfs, self.out_bams, self.annotated_vcfs,
                        self.longshot_vcfs, self.consensus_vcfs
                    )
            ],
            self.shell_files
        ).multi_run(self.cpu)
        return self


if __name__ == "__main__":
    Medaka().medaka()
