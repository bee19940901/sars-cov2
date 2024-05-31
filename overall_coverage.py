"""
覆盖度分析
"""
from pathlib import Path
from argparse import ArgumentParser
from collections import defaultdict

import pandas as pd
from scipy.stats import mannwhitneyu
from plotnine import *

from configs import PLATFORM2, PLATFORM3


class OverallCoverage:

    def __init__(self):

        ap = ArgumentParser()
        ap.add_argument("-mcf", "--mgi_coverage_file", type=str, required=True)
        ap.add_argument("-ncf", "--np_coverage_file", type=str, required=True)
        ap.add_argument("-od", "--out_dir", type=str, required=True)
        ag = ap.parse_args()

        self.mgi_coverage_file = Path(ag.mgi_coverage_file).absolute()
        self.np_coverage_file = Path(ag.np_coverage_file).absolute()
        self.out_dir = Path(ag.out_dir).absolute()
        self.out_dir.mkdir(parents=True, exist_ok=True)

        self.data_df = pd.DataFrame()
        self.data_xls = self.out_dir.joinpath("data_xls")

        self.stats_df = pd.DataFrame()
        self.stats_xls = self.out_dir.joinpath("stats_xls")

        self.coverage_png = self.out_dir.joinpath("coverages.png")
        self.coverage_pdf = self.out_dir.joinpath("coverages.pdf")

    def parse(self):
        self.data_df = (
            pd.concat(
                [
                    pd.read_table(self.mgi_coverage_file).assign(platform=PLATFORM2),
                    pd.read_table(self.np_coverage_file).assign(platform=PLATFORM3)
                ]
            ).filter(["sample", "platform", "coverage"])
        )
        self.data_df.to_csv(self.data_xls, sep="\t", header=True, index=False)
        return self

    def stats(self):
        stats_dict = defaultdict(list)
        s, p = mannwhitneyu(
            self.data_df.query("platform == @PLATFORM2")["coverage"],
            self.data_df.query("platform == @PLATFORM3")["coverage"],
        )
        mgi_n = self.data_df.query("platform == @PLATFORM2").shape[0]
        np_n = self.data_df.query("platform == @PLATFORM3").shape[0]
        stats_dict["TYPE"].append("coverage")
        stats_dict["mgi_n"].append(mgi_n)
        stats_dict["np_n"].append(np_n)
        stats_dict["statistic"].append(s)
        stats_dict["p-value"].append(p)
        self.stats_df = pd.DataFrame(stats_dict)
        self.stats_df.to_csv(self.stats_xls, sep="\t", header=True, index=False)
        return self

    def plots(self):
        p = (
            ggplot(**{"data": self.data_df, "mapping": aes(**{"x": "platform", "y": "coverage", "fill": "platform"})}) +
            geom_violin() +
            geom_boxplot(
                width=0.1, position=position_dodge(0.9), outlier_size=3, outlier_shape="*", outlier_color="orange"
            ) +
            ylab("Coverage") +
            scale_y_continuous(
                limits=(0, 105),
                breaks=range(0, 110, 10),
                labels=[f"{i}%" for i in range(0, 110, 10)]
            ) +
            theme_classic() +
            theme(
                legend_position="none",
                axis_title_x=element_blank(),
                text=element_text(family="Times New Roman"),
                title=element_text(family="Times New Roman")
            )
        )
        p.save(self.coverage_pdf, width=8, height=6)
        p.save(self.coverage_png, width=8, height=6, dpi=400)
        return self


if __name__ == "__main__":
    (
        OverallCoverage()
        .parse()
        .stats()
        .plots()
    )
