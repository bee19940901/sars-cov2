from pathlib import Path
import re

from utils import Cmds
from configs import DATABASE, SOFTWARE


class BcftoolsCall:

    def __init__(self, args):

        self.cpu = int(args.cpu)
        self.in_dir = Path(args.in_dir).absolute()
        self.out_dir = Path(args.out_dir).absolute()

        self.sh_dir = self.out_dir.joinpath("work_sh")
        self.sh_dir.mkdir(parents=True, exist_ok=True)
        self.sorted_bams = list(self.in_dir.glob("*.sorted.bam"))
        self.samples = [re.sub(r"\.sorted\.bam$", "", i.name) for i in self.sorted_bams]
        self.raw_vcfs = [self.out_dir.joinpath(f"{i}.raw.vcf") for i in self.samples]
        self.good_vcfs = [self.out_dir.joinpath(f"{i}.good.vcf") for i in self.samples]

    def bcftools_call(self):
        Cmds(
            [
                f"{SOFTWARE.joinpath('bin', 'bcftools')} mpileup --annotate FORMAT/AD --max-depth 1000000 --threads 8 "
                f"-f {DATABASE.joinpath('cons_ref', 'sars-cov2.ref.fa')} {sorted_bam} | "
                f"{SOFTWARE.joinpath('bin', 'bcftools')} call --multiallelic-caller  --variants-only "
                f"--ploidy 1 --threads 8  -o {raw_vcf} "
                f"&& {SOFTWARE.joinpath('bin', 'bcftools')} view -i 'QUAL >= 10 && INFO/DP >= 10' "
                f"{raw_vcf} > {good_vcf}\n"
                for sorted_bam, raw_vcf, good_vcf in zip(self.sorted_bams, self.raw_vcfs, self.good_vcfs)
            ],
            [
                self.sh_dir.joinpath(f"{sample}.sh")
                for sample in self.samples
            ]
        ).multi_run(self.cpu)


def main(args):
    BcftoolsCall(args).bcftools_call()
