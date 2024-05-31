"""
统计二三代CT值和reads数量 coverage之间的关系
"""
from argparse import ArgumentParser
from pathlib import Path
from collections import defaultdict

import pandas as pd
from plotnine import *
from scipy.stats import linregress

from configs import PLATFORM2, PLATFORM3


class OverallCT:

    def __init__(self):

        ap = ArgumentParser()
        ap.add_argument("-mcf", "--mgi_coverage_file", type=str, required=True)
        ap.add_argument("-ncf", "--np_coverage_file", type=str, required=True)
        ap.add_argument("-csf", "--ct_score_file", type=str, required=True)
        ap.add_argument("-od", "--out_dir", type=str, required=True)
        ag = ap.parse_args()

        self.mgi_coverage_file = Path(ag.mgi_coverage_file).absolute()
        self.np_coverage_file = Path(ag.np_coverage_file).absolute()
        self.ct_score_file = Path(ag.ct_score_file).absolute()
        self.out_dir = Path(ag.out_dir).absolute()
        self.out_dir.mkdir(parents=True, exist_ok=True)
        self.mgi_data_df = pd.DataFrame()
        self.np_data_df = pd.DataFrame()
        self.stat_df = pd.DataFrame()
        self.data_df = pd.DataFrame()

    def data(self):
        ct_score_df = pd.read_csv(
            self.ct_score_file,
            sep="\t"
        )
        self.mgi_data_df = (
            pd
            .read_csv(self.mgi_coverage_file, sep="\t")
            .assign(platform=PLATFORM2)
            .loc[:, ["sample", "platform", "numreads", "coverage"]]
            .merge(ct_score_df, on="sample", how="inner")
        )
        self.np_data_df = (
            pd
            .read_csv(self.np_coverage_file, sep="\t")
            .assign(platform=PLATFORM3)
            .loc[:, ["sample", "platform", "numreads", "coverage"]]
            .merge(ct_score_df, on="sample", how="inner")
        )
        self.data_df = pd.concat(
            [
                self.mgi_data_df,
                self.np_data_df
            ]
        )
        self.data_df.to_csv(
            self.out_dir.joinpath("data.xls"),
            sep="\t",
            header=True,
            index=False
        )
        return self

    def stat(self):
        stat_dict = defaultdict(list)
        for p, d in zip([PLATFORM2, PLATFORM3], [self.mgi_data_df, self.np_data_df]):
            for x in ["ORF1ab", "N"]:
                for y in ["numreads", "coverage"]:
                    res = linregress(d[x], d[y])
                    stat_dict["platform"].append(p)
                    stat_dict["x"].append(x)
                    stat_dict["y"].append(y)
                    stat_dict["r-value"].append(res.rvalue)
                    stat_dict["p-value"].append(res.pvalue)
        self.stat_df = pd.DataFrame(stat_dict)
        self.stat_df.to_csv(self.out_dir.joinpath("stat.xls"), sep="\t", header=True, index=False)
        return self

    def foo(self, ct, value):
        p = (
                ggplot(**{"data": self.data_df, "mapping": aes(**{"x": ct, "y": value, "color": "platform"})}) +
                geom_point(size=2, shape='•') +
                geom_smooth(method='lm', se=True, linetype="dotted") +
                theme_classic() +
                theme(
                    legend_title=element_blank(),
                    text=element_text(family="Times New Roman"),
                    title=element_text(family="Times New Roman")
                )
        )
        if value == "coverage":
            p += ylab("Coverage")
            p += scale_y_continuous(
                limits=(0, 105),
                breaks=range(0, 110, 10),
                labels=[f"{i}%" for i in range(0, 110, 10)]
            )
        else:
            p += ylab("Reads Number")
            p += scale_y_continuous(
                limits=(0, None),
                breaks=range(0, self.data_df[value].max() + 1000000, 1000000),
                labels=[f"{i/1000000:.0f}M" for i in range(0, self.data_df[value].max() + 1000000, 1000000)]
            )
        p.save(self.out_dir.joinpath(f"{ct}-{value}.png"), width=8, height=6, dpi=400)
        p.save(self.out_dir.joinpath(f"{ct}-{value}.pdf"), width=8, height=6)

    def plot(self):
        self.foo("ORF1ab", "numreads",)
        self.foo("N", "numreads")
        self.foo("ORF1ab", "coverage",)
        self.foo("N", "coverage",)
        return self


if __name__ == "__main__":
    (
        OverallCT()
        .data()
        .stat()
        .plot()
    )
