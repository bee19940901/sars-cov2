"""
获取BMC文献中table1中的内容
"""
import json
from pathlib import Path

import pandas as pd


class Overall:
    def __init__(self, args):

        self.coverage_file = Path(args.coverage_file).absolute()
        self.stat_file = Path(args.stat_file).absolute()
        self.overall_file = Path(args.overall_file).absolute()

        self.out_dir = self.overall_file.parent
        self.out_dir.mkdir(parents=True, exist_ok=True)

    def conclude(self):
        coverage_df = pd.read_table(self.coverage_file)
        stat_df = pd.read_table(self.stat_file, thousands=",")
        print(stat_df['sum_len'].sum())
        out_dict = {
            "data_obtained": f"{(stat_df['sum_len'].sum() / 1000000):.1f}MB",
            "average_coverage": f"{coverage_df['coverage'].mean():.2f}%",
            "passing_gisaid_qc": f"{((coverage_df['coverage'] >= 80).sum() / (coverage_df.shape[0])) * 100:.2f}%",
            "q30_score": f"{stat_df['Q30(%)'].mean():.1f}%",
        }
        with open(self.overall_file, "w", encoding="utf-8") as fw:
            json.dump(out_dict, fw)
        return self


def main(args):
    Overall(args).conclude()
