from argparse import ArgumentParser
from pathlib import Path
from collections import defaultdict
import re
import math

import pandas as pd
from scipy.stats import mannwhitneyu
from plotnine import *

from configs import DATABASE, PLATFORM2, PLATFORM3


class GenesCoverage:

    BED = DATABASE.joinpath("cons_ref", "sars.gene.bed")

    def __init__(self):
        ap = ArgumentParser()
        ap.add_argument("-mdd", "--mgi_depth_dir", type=str, required=True)
        ap.add_argument("-ndd", "--np_depth_dir", type=str, required=True)
        ap.add_argument("-od", "--out_dir", type=str, required=True)
        ap.add_argument("-mmd", "--mgi_min_depth", type=int, default=30)
        ap.add_argument("-nmd", "--np_min_depth", type=int, default=30)
        ag = ap.parse_args()
        self.mgi_min_depth = int(ag.mgi_min_depth)
        self.np_min_depth = int(ag.np_min_depth)
        self.out_dir = Path(ag.out_dir).absolute()
        self.mgi_depth_dir = Path(ag.mgi_depth_dir).absolute()
        self.np_depth_dir = Path(ag.np_depth_dir).absolute()
        self.data_df = pd.DataFrame()
        self.genes_dict = defaultdict(tuple)
        self.stats_dict = defaultdict(list)
        self.stats_df = pd.DataFrame()
        self.mgi_depths = list(self.mgi_depth_dir.glob("*.depth"))
        self.mgi_samples = [re.sub(r"\.depth$", "", i.name) for i in self.mgi_depths]
        self.np_depths = list(self.np_depth_dir.glob("*.depth"))
        self.np_samples = [re.sub(r"\.depth$", "", i.name) for i in self.np_depths]
        self.out_dir.mkdir(parents=True, exist_ok=True)

    def parse(self):
        with open(self.BED, "r", encoding="utf-8") as fr:
            for line in fr.readlines():
                if line := line.strip():
                    items = line.split("\t")
                    self.genes_dict[items[3]] = (int(items[1]), int(items[2]))
        data_dict = defaultdict(list)
        for g, r in self.genes_dict.items():
            gene_len = r[1] - r[0] + 1
            for s, d in zip(self.mgi_samples, self.mgi_depths):
                with open(d, 'r', encoding="utf-8") as fd:
                    init_n = 0
                    init_d = 0
                    for line in fd:
                        if line := line.strip():
                            its = line.split("\t")
                            if r[0] <= int(its[1]) <= r[1]:
                                if int(its[2]) > 0:
                                    init_d += int(its[2])
                                    if int(its[2]) >= self.mgi_min_depth:
                                        init_n += 1
                                else:
                                    continue
                    coverage = init_n / gene_len
                    meandepth = init_d / gene_len
                    data_dict["sample"].append(s)
                    data_dict["platform"].append(PLATFORM2)
                    data_dict["gene"].append(g)
                    data_dict["coverage"].append(coverage)
                    data_dict["meandepth"].append(meandepth)
            for s, d in zip(self.np_samples, self.np_depths):
                with open(d, 'r', encoding="utf-8") as fd:
                    init_n = 0
                    init_d = 0
                    for line in fd:
                        if line := line.strip():
                            its = line.split("\t")
                            if r[0] <= int(its[1]) <= r[1]:
                                if int(its[2]) > 0:
                                    init_d += int(its[2])
                                    if int(its[2]) >= self.np_min_depth:
                                        init_n += 1
                                else:
                                    continue
                    coverage = init_n / gene_len
                    meandepth = init_d / gene_len
                    data_dict["sample"].append(s)
                    data_dict["platform"].append(PLATFORM3)
                    data_dict["gene"].append(g)
                    data_dict["coverage"].append(coverage)
                    data_dict["meandepth"].append(meandepth)
        self.data_df = pd.DataFrame(data_dict)
        self.data_df.to_csv(self.out_dir.joinpath("data.xls"), sep="\t", header=True, index=False)
        return self

    def stats(self):
        for gene in self.genes_dict:
            coverage_s, coverage_p = (
                mannwhitneyu(
                    self.data_df.query("gene == @gene and platform == @PLATFORM2")["coverage"],
                    self.data_df.query("gene == @gene and platform == @PLATFORM3")["coverage"],
                    alternative="two-sided"
                )
            )
            meandepth_s, meandepth_p = (
                mannwhitneyu(
                    self.data_df.query("gene == @gene and platform == @PLATFORM2")["meandepth"],
                    self.data_df.query("gene == @gene and platform == @PLATFORM3")["meandepth"],
                    alternative="two-sided"
                )
            )
            mgi_n = self.data_df.query("gene == @gene and platform == @PLATFORM2").shape[0]
            np_n = self.data_df.query("gene == @gene and platform == @PLATFORM3").shape[0]
            self.stats_dict["gene"].append(gene)
            self.stats_dict[f"{PLATFORM2}_n"].append(mgi_n)
            self.stats_dict[f"{PLATFORM3}_n"].append(np_n)
            self.stats_dict["coverage_statistic"].append(coverage_s)
            self.stats_dict["coverage_p-value"].append(coverage_p)
            self.stats_dict["meandepth_statistic"].append(meandepth_s)
            self.stats_dict["meandepth_p-value"].append(meandepth_p)
        self.stats_df = pd.DataFrame(self.stats_dict)
        self.stats_df.to_csv(self.out_dir.joinpath("stats.xls"), sep="\t", header=True, index=False)
        return self

    def plot_depth(self):
        data = self.data_df.assign(meandepth=self.data_df["meandepth"].apply(lambda x: math.log2(x + 1)))
        p = (
                ggplot(**{"data": data, "mapping": aes(**{"x": "gene", "y": "meandepth", "fill": "platform"})}) +
                geom_violin() +
                geom_boxplot(
                    width=0.1, position=position_dodge(0.9), outlier_shape="*", outlier_color="orange"
                ) +
                ylab("log2(Meandepth+1)") +
                scale_x_discrete(limits=list(self.genes_dict.keys())) +
                scale_y_continuous(
                    limits=(0, data.meandepth.max() * 1.1),
                ) +
                theme_classic() +
                theme(
                    axis_title_x=element_blank(),
                    legend_title=element_blank(),
                    text=element_text(family="Times New Roman"),
                    title=element_text(family="Times New Roman")
                )
        )
        p.save(self.out_dir.joinpath(f"meandepth.png"), width=16, height=6, dpi=400)
        p.save(self.out_dir.joinpath(f"meandepth.pdf"), width=16, height=6)
        return self

    def plot_coverage(self):
        data = (
            self.data_df
            .groupby(by=["platform", "gene"], as_index=False)
            .agg({"meandepth": "mean", "coverage": "mean"})
            .assign(coverage=lambda x: x["coverage"]*100)
        )
        p = (
                ggplot(**{"data": data, "mapping": aes(**{"x": "gene", "y": "coverage", "fill": "platform"})}) +
                geom_col(position=position_dodge(0.75)) +
                scale_x_discrete(limits=list(self.genes_dict.keys())) +
                scale_y_continuous(
                    limits=(0, 105),
                    breaks=range(0, 110, 10),
                    labels=[f"{i}%" for i in range(0, 110, 10)]
                ) +
                ylab("Coverage") +
                theme_classic() +
                theme(
                    axis_title_x=element_blank(),
                    legend_title=element_blank(),
                    text=element_text(family="Times New Roman"),
                    title=element_text(family="Times New Roman")
                )
        )
        p.save(self.out_dir.joinpath(f"coverage.png"), width=16, height=6, dpi=400)
        p.save(self.out_dir.joinpath(f"coverage.pdf"), width=16, height=6)
        return self


if __name__ == "__main__":
    (
        GenesCoverage()
        .parse()
        .stats()
        .plot_coverage()
        .plot_depth()
    )
