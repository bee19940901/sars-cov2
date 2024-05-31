from argparse import ArgumentParser
from pathlib import Path
import re

import pandas as pd
import numpy as np

from utils import Cmds
from configs import BCFTOOLS,  REF_FA, BLACKLIST


class BcftoolsCall:

    def __init__(self):
        ap = ArgumentParser()
        ap.add_argument("--in_dir", type=str, required=True, help="trimmed_bams/*.sorted.bam")
        ap.add_argument("--out_dir", type=str, required=True, help="consensus_vcfs/*.consensus.vcf")
        ap.add_argument("--cpu", type=int, default=8)
        self.args = ap.parse_args()
        self.in_dir = Path(self.args.in_dir).absolute()
        self.out_dir = Path(self.args.out_dir).absolute()
        self.sh_dir = self.out_dir.joinpath("work_sh")
        self.sh_dir.mkdir(parents=True, exist_ok=True)
        self.sorted_bams = list(
            self.in_dir.glob("*.sorted.bam")
        )
        self.samples = [
            re.sub(r"\.sorted\.bam$", "", i.name)
            for i in self.sorted_bams
        ]
        self.temp_dirs = [
            self.out_dir.joinpath(i)
            for i in self.samples
        ]
        self.raw_vcfs = [
            i.joinpath("raw.vcf")
            for i in self.temp_dirs
        ]
        self.good_vcfs = [
            i.joinpath("good.vcf")
            for i in self.temp_dirs
        ]
        self.consensus_vcfs = [
            i.joinpath("consensus.vcf")
            for i in self.temp_dirs
        ]
        for i in self.temp_dirs:
            i.mkdir(parents=True, exist_ok=True)

    def bcftools_call(self):
        Cmds(
            [
                f"{BCFTOOLS} mpileup --no-BAQ --min-BQ 0 --annotate FORMAT/AD --threads 8 "
                f"-f {REF_FA} {sorted_bam} | "
                f"{BCFTOOLS} call -c  --variants-only "
                f"--ploidy 1 --threads 8  -o {raw_vcf} "
                f"&& {BCFTOOLS} view -i 'QUAL >= 20 && INFO/DP >= 30' "
                f"{raw_vcf} > {good_vcf}\n"
                for sorted_bam, raw_vcf, good_vcf
                in zip(self.sorted_bams, self.raw_vcfs, self.good_vcfs)
            ],
            [
                self.sh_dir.joinpath(f"{sample}.sh")
                for sample in self.samples
            ]
        ).multi_run(self.args.cpu)
        return self

    def merge(self):

        def func(x):
            if "," in x["ALT"]:
                return "FAIL"
            matches = re.search(r";DP4=(\d+),(\d+),(\d+),(\d+);", x["INFO"])
            rate = (
                (int(matches.group(3)) + int(matches.group(4))) /
                (int(matches.group(1)) + int(matches.group(2)) + int(matches.group(3)) + int(matches.group(4)))
            )
            if rate >= 0.8:
                return "PASS"
            else:
                return "FAIL"

        black_sites = list()
        dfs = list()
        with open(BLACKLIST, "r") as fr:
            for line in fr:
                if line := line.strip():
                    items = line.split("\t")
                    black_sites += list(range(int(items[1]), int(items[2]) + 1))
        for sample, pass_vcf, consensus_vcf in zip(self.samples, self.good_vcfs, self.consensus_vcfs):
            df = (
                pd
                .read_csv(
                    pass_vcf, sep="\t", comment="#", header=None,
                    usecols=[1, 3, 4, 7], names=["POS", "REF", "ALT", "INFO"]
                )
                .assign(
                    BLACKED=lambda x: x["POS"].apply(lambda y: "YES" if int(y) in black_sites else "NO"),
                    SAMPLE=sample,
                    FILTER=lambda x: x.apply(func, axis=1)
                )
            )
            df['TYPE'] = np.where((df['REF'].str.len() == 1) & (df['ALT'].str.len() == 1), 'SNV', 'INDEL')
            df.to_csv(consensus_vcf, sep="\t", header=True, index=False)
            dfs.append(df.query("FILTER == 'PASS'").drop(columns=["FILTER", "INFO"]))
        pd.concat(dfs).to_csv(self.out_dir.joinpath("all.merged.consensus.vcf"), sep="\t", header=True, index=False)
        return self


if __name__ == "__main__":
    BcftoolsCall().bcftools_call().merge()
