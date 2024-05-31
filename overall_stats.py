"""
统计二三代结果的基本信息
"""
from pathlib import Path
from argparse import ArgumentParser
from collections import defaultdict, namedtuple

import pandas as pd

from configs import PLATFORM2, PLATFORM3

class OverallStats:

    def __init__(self):

        ap = ArgumentParser()
        ap.add_argument("-mcf", "--mgi_coverage_file", type=str, required=True)
        ap.add_argument("-ncf", "--np_coverage_file", type=str, required=True)
        ap.add_argument("-msf", "--mgi_stats_file", type=str, required=True)
        ap.add_argument("-nsf", "--np_stats_file", type=str, required=True)
        ap.add_argument("-od", "--out_dir", type=str, required=True)
        ag = ap.parse_args()

        self.mgi_coverages_file = Path(ag.mgi_coverage_file).absolute()
        self.np_coverage_file = Path(ag.np_coverage_file).absolute()
        self.mgi_stats_file = Path(ag.mgi_stats_file).absolute()
        self.np_stats_file = Path(ag.np_stats_file).absolute()
        self.out_dir = Path(ag.out_dir).absolute()

        self.out_dir.mkdir(parents=True, exist_ok=True)
        self.overall_stats_file = self.out_dir.joinpath("overall.stats.xls")

    def overall_stats(self):
        Task = namedtuple("Task", ["platform", "stats_file", "coverage_file"])
        tasks = [
            Task(PLATFORM2, self.mgi_stats_file, self.mgi_coverages_file),
            Task(PLATFORM3, self.np_stats_file, self.np_coverage_file)
        ]
        overall_stats_dict = defaultdict(list)
        for t in tasks:
            stats_df = pd.read_table(t.stats_file, thousands=",")
            coverage_df = pd.read_table(t.coverage_file, thousands=",")
            overall_stats_dict["platform"].append(t.platform)
            overall_stats_dict["n"].append(coverage_df.shape[0])
            overall_stats_dict["data_obtained"].append(f"{stats_df['sum_len'].sum() / 1000000:.2f}MB")
            overall_stats_dict["average_coverage"].append(f"{coverage_df['coverage'].mean():.2f}%")
            overall_stats_dict["passing_gisaid_qc"].append(
                f"{((coverage_df['coverage'] >= 80).sum() / (coverage_df.shape[0])) * 100:.2f}%"
            )
            overall_stats_dict["average_q30_score"].append(f"{stats_df['Q30(%)'].mean():.2f}%")
            overall_stats_dict["average_q20_score"].append(f"{stats_df['Q20(%)'].mean():.2f}%")
            overall_stats_dict["average_depth"].append(f"{coverage_df['meandepth'].mean():.2f}")
            overall_stats_dict["total_reads"].append(f"{stats_df['num_seqs'].sum():.2f}")
            overall_stats_dict["average_reads"].append(f"{stats_df['num_seqs'].mean():.2f}")
            overall_stats_dict["average_base"].append(f"{stats_df['sum_len'].mean() / 1000000:.2f}MB")
            overall_stats_dict["average_N50"].append(f"{stats_df['N50'].mean():.2f}")
        pd.DataFrame(overall_stats_dict).to_csv(self.overall_stats_file, sep="\t", header=True, index=False)
        return self


if __name__ == "__main__":
    OverallStats().overall_stats()
