import re
from pathlib import Path
import pandas as pd

from configs import DATABASE, SOFTWARE
from utils import Cmds, Results


class Coverage:

    def __init__(self, args):

        self.cpu = int(args.cpu)
        self.in_dir = Path(args.in_dir).absolute()
        self.out_dir = Path(args.out_dir).absolute()

        self.merged_file = self.out_dir.joinpath("all.coverage.merged.tsv")
        self.sh_dir = self.out_dir.joinpath("work_sh")
        self.sh_dir.mkdir(parents=True, exist_ok=True)
        self.in_bams = list(self.in_dir.glob("*.sorted.bam"))
        self.samples = [re.sub(r"\.sorted\.bam$", "", i.name) for i in self.in_bams]
        self.out_coverages = [self.out_dir.joinpath(f"{i}.coverage") for i in self.samples]

    def coverage(self):
        Cmds(
            [
                f"{SOFTWARE.joinpath('bin','samtools')} coverage -d 1000000 {in_bam} > {out_coverage}\n"
                for in_bam, out_coverage in zip(self.in_bams, self.out_coverages)
            ],
            [
                self.sh_dir.joinpath(f"{sample}.sh")
                for sample in self.samples
            ]
        ).multi_run(self.cpu)
        return self

    def merge(self):
        dfs = list()
        for sample, i in zip(self.samples, self.out_coverages):
            df = pd.read_csv(i, sep=r"\s+")
            df["sample"] = sample
            dfs.append(df)
        all_df = pd.concat(dfs, axis=0)
        all_df.sort_values(by="coverage", ascending=False)\
            .to_csv(self.merged_file, sep="\t",header=True, index=False)
        return self


def main(args):
    Coverage(args).coverage().merge()

