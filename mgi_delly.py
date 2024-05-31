from argparse import ArgumentParser
from pathlib import Path
from re import sub, search

import pandas as pd

from configs import DELLY, REF_FA, PLATFORM3, PLATFORM2
from utils import Cmds, Results


class MgiDelly:

    def __init__(self):

        ap = ArgumentParser()
        ap.add_argument("--in_bam_dir", type=str, required=True)
        ap.add_argument("--out_vcf_dir", type=str, required=True)
        ap.add_argument("--cpu", type=int, default=1)
        ap.add_argument("--platform", type=str, required=True, choices=[PLATFORM2, PLATFORM3])

        self.args = ap.parse_args()
        self.in_bam_dir = Path(self.args.in_bam_dir).absolute()
        self.in_bam_files = list(
            self.in_bam_dir.glob("*.bam")
        )
        self.samples = [
            sub(r"(?:\.sorted)?\.bam$", "", i.name)
            for i in self.in_bam_files
        ]
        self.out_vcf_dir = Path(self.args.out_vcf_dir).absolute()
        self.out_vcf_files = [
            self.out_vcf_dir.joinpath(f"{s}.vcf")
            for s in self.samples
        ]
        self.out_sh_dir = self.out_vcf_dir.joinpath("work_sh")
        self.out_sh_dir.mkdir(parents=True, exist_ok=True)

    def run(self):
        (
            Results(self.in_bam_files)
            .check_exists()
            .check_empty()
        )
        Cmds(
            [
                f"{DELLY} call -g {REF_FA} {in_bam_file} > {out_vcf_file}\n"
                if {self.args.platform} == PLATFORM2
                else
                f"{DELLY} lr -y ont -g {REF_FA} {in_bam_file} > {out_vcf_file}\n"
                for in_bam_file, out_vcf_file
                in zip(self.in_bam_files, self.out_vcf_files)
            ],
            [
                self.out_sh_dir.joinpath(f"{s}.sh")
                for s in self.samples
            ]
        ).multi_run(self.args.cpu)
        (
            Results(self.out_vcf_files)
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
                    .assign(SVTYPE=lambda x: x["INFO"].apply(lambda y: search(r"SVTYPE=(\w+);", y).group(1)))
                    .assign(END=lambda x: x["INFO"].apply(lambda y: search(r";END=(\d+);", y).group(1)))
                    .filter(["SAMPLE", "ID", "START", "END", "SVTYPE"])
                )
                for sample, out_vcf_file in zip(self.samples, self.out_vcf_files)
            ]
        ).sort_values(by="SAMPLE").to_csv(
            self.out_vcf_dir.joinpath("all.merged.xls"),
            sep="\t",
            header=True,
            index=False
        )
        return self


if __name__ == "__main__":
    (
        MgiDelly()
        .run()
        .merge()
    )
