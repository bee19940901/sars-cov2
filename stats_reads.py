"""
统计过滤得到的新冠reads信息
"""
from pathlib import Path
import re
import pandas as pd

from utils import Cmds, Results
from configs import SOFTWARE


class StatsReads:

    def __init__(self, args):
        """
        :param args: 命令行参数
        """
        self.cpu = int(args.cpu)
        self.in_dir = Path(args.in_dir).absolute()
        self.out_dir = Path(args.out_dir).absolute()

        self.sh_dir = self.out_dir.joinpath("works_sh")
        self.sh_dir.mkdir(parents=True, exist_ok=True)

        self.merged_file = self.out_dir.joinpath("all.stats.merged.tsv")
        self.in_fqs = list(self.in_dir.glob("*.fq.gz"))
        self.samples = [re.sub(r"\.fq\.gz$", "", i.name) for i in self.in_fqs]
        self.out_stats = [self.out_dir.joinpath(f"{s}.stat") for s in self.samples]

    def stats_reads(self):
        Cmds(
            [
                f"{SOFTWARE.joinpath('bin', 'seqkit')} stats  -a -j 8 {in_fq} > {out_stat} "
                for in_fq, out_stat in zip(self.in_fqs, self.out_stats)
            ],
            [
                self.sh_dir.joinpath(f"{sample}.sh")
                for sample in self.samples
            ]
        ).multi_run(self.cpu)
        Results(self.out_stats).check_exists().check_empty()
        return self

    def merge_stats(self):
        """
        将所有的统计结果进行合并
        :return:
        """
        stat_dfs = [pd.read_csv(i, sep=r"\s+") for i in self.out_stats]
        all_stats_df = pd.concat(stat_dfs, axis=0)
        all_stats_df.to_csv(self.merged_file, sep="\t", header=True, index=False)
        return self


def main(args):
    """
    :param args: 命令哈参数
    :return:
    """
    StatsReads(args).stats_reads().merge_stats()
