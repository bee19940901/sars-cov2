from argparse import ArgumentParser
from pathlib import Path
from re import sub, search

import pandas as pd

from utils import Cmds, Results
from configs import NGMLR, SNIFFLES, SAMTOOLS109, IVAR, REF_FA, PRIMERS_BED


class NpSniffles:

    def __init__(self):

        ap = ArgumentParser()
        ap.add_argument("--in_fqs_dir", type=str, required=True)
        ap.add_argument("--out_dir", type=str, required=True)
        ap.add_argument("--cpu", type=int, default=1)

        self.args = ap.parse_args()
        self.in_fqs_dir = Path(self.args.in_fqs_dir).absolute()
        self.out_dir = Path(self.args.out_dir).absolute()
        self.sh_dir = self.out_dir.joinpath("work_sh")
        self.sh_dir.mkdir(parents=True, exist_ok=True)
        self.in_fqs = list(
            self.in_fqs_dir.glob("*.fq.gz")
        )
        self.samples = [
            sub(r"(?:\.clean)?\.fq.gz$", "", i.name)
            for i in self.in_fqs
        ]
        self.out_bams = [
            self.out_dir.joinpath(f"{s}.bam")
            for s in self.samples
        ]
        self.out_vcfs = [
            self.out_dir.joinpath(f"{s}.vcf")
            for s in self.samples
        ]

    def run(self):
        Cmds(
            [
                f"{NGMLR} -t 8 -r {REF_FA} -q {in_fq} -x ont | "
                f"{SAMTOOLS109} view -@ 8 -F 4 -h | "
                f"{IVAR} trim -b {PRIMERS_BED} -e -q 0 | "
                f"{SAMTOOLS109} sort -@ 8 -o {out_bam} && "
                f"{SNIFFLES} -m {out_bam} -v {out_vcf} -l 10 -s 20\n"
                for in_fq, out_bam, out_vcf in zip(self.in_fqs, self.out_bams, self.out_vcfs)
            ],
            [
                self.sh_dir.joinpath(f"{s}.sh")
                for s in self.samples
            ]
        ).multi_run(self.args.cpu)
        (
            Results(self.out_bams + self.out_vcfs)
            .check_exists()
            .check_empty()
        )
        return self

    def merge(self):
        pd.concat(
            [
                (
                    pd.read_table(
                        out_vcf_file,
                        sep="\t",
                        comment="#",
                        header=None,
                        usecols=[2, 1, 7],
                        names=["START", "ID", "INFO"])
                    .assign(SAMPLE=sample)
                    .assign(RE=lambda x: x["INFO"].apply(lambda y: search(r"RE=(\d+)$", y).group(1)))
                    .assign(SVTYPE=lambda x: x["INFO"].apply(lambda y: search(r";SVTYPE=(\w+);", y).group(1)))
                    .assign(END=lambda x: x["INFO"].apply(lambda y: search(r";END=(\d+);", y).group(1)))
                    .assign(SVLEN=lambda x: x["INFO"].apply(lambda y: search(r";SVLEN=-?(\d+);", y).group(1)))
                    .filter(["SAMPLE", "ID", "START", "END", "SVLEN", "SVTYPE", "RE"])
                )
                for sample, out_vcf_file in zip(self.samples, self.out_vcfs)
            ]
        ).sort_values(by="SAMPLE").to_csv(self.out_dir.joinpath("all.merged.xls"), sep="\t", header=True, index=False)
        return self


if __name__ == "__main__":
    (
        NpSniffles()
        .run()
        .merge()
    )
