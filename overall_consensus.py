"""
对nextclade基因组评估结果进行整理
"""
from argparse import ArgumentParser
from pathlib import Path
from collections import defaultdict
from math import log2

import pandas as pd
import plotnine as pn
from scipy.stats import mannwhitneyu


from configs import PLATFORM2, PLATFORM3


class OverallConsensus:

    def __init__(self):

        ap = ArgumentParser()
        ap.add_argument("-mnf", "--mgi_nextclade_csv", type=str, required=True)
        ap.add_argument("-nnf", "--np_nextclade_csv", type=str, required=True)
        ap.add_argument("-od", "--out_dir", type=str, required=True)
        ag = ap.parse_args()

        self.mgi_nextclade_csv = Path(ag.mgi_nextclade_csv).absolute()
        self.np_nextclade_csv = Path(ag.np_nextclade_csv).absolute()
        self.out_dir = Path(ag.out_dir).absolute()
        self.out_dir.mkdir(parents=True, exist_ok=True)

        self.data_df = pd.DataFrame()
        self.stats_df = pd.DataFrame()

        self.data_file = self.out_dir.joinpath("consensus.xls")
        self.stats_file = self.out_dir.joinpath("stats.xls")
        self.png_file = self.out_dir.joinpath("consensus.png")
        self.pdf_file = self.out_dir.joinpath("consensus.pdf")

    def parse(self):
        mgi_nextclade_df = (
            pd
            .read_csv(
                self.mgi_nextclade_csv,
                sep=";",
                usecols=["seqName", "qc.overallScore", "qc.overallStatus"]
            )
            .rename(columns={"qc.overallScore": "score", "qc.overallStatus": "status"})
            .assign(platform=PLATFORM2)
        )
        np_nextclade_df = (
            pd
            .read_csv(
                self.np_nextclade_csv,
                sep=";",
                usecols=["seqName", "qc.overallScore", "qc.overallStatus"]
            )
            .rename(columns={"qc.overallScore": "score", "qc.overallStatus": "status"})
            .assign(platform=PLATFORM3)
        )
        total_df = pd.concat([mgi_nextclade_df, np_nextclade_df]).assign(status="total")
        self.data_df = pd.concat([mgi_nextclade_df, np_nextclade_df, total_df])

        self.data_df.to_csv(self.data_file, sep="\t", header=True, index=False)
        return self

    def stats(self):
        stats_dict = defaultdict(list)
        for i in ["bad", "mediocre", "good", "total"]:
            s, p = mannwhitneyu(
                self.data_df.query("platform == @PLATFORM2 and status == @i")["score"],
                self.data_df.query("platform == @PLATFORM3 and status == @i")["score"],
            )
            n2 = self.data_df.query("platform == @PLATFORM2 and status == @i")["score"].shape[0]
            n3 = self.data_df.query("platform == @PLATFORM3 and status == @i")["score"].shape[0]
            stats_dict["status"].append(i)
            stats_dict[f"{PLATFORM2}_n"].append(n2)
            stats_dict[f"{PLATFORM3}_n"].append(n3)
            stats_dict["statistic"].append(s)
            stats_dict["p-value"].append(p)
        self.stats_df = pd.DataFrame(stats_dict)
        self.data_df["score"] = self.data_df["score"].apply(lambda x: log2(x + 1))
        self.stats_df.to_csv(self.stats_file, sep="\t", header=True, index=False)
        return self

    def plot(self):
        p = (
            pn.ggplot(
                **{
                    "data": self.data_df,
                    "mapping": pn.aes(**{"x": "status", "y": "score", "fill": "platform"})
                }
            ) +
            pn.geom_violin() +
            pn.geom_boxplot(
                width=0.1,
                position=pn.position_dodge(0.9),
                outlier_size=3,
                outlier_shape="*",
                outlier_color="orange"
            ) +
            pn.ylab("log2(Score+1)") +
            pn.scale_x_discrete(
                limits=["total", "good", "mediocre", "bad"],
                labels=["Total", "Good", "Mediocre", "Bad"]
            ) +
            pn.scale_y_continuous(
                limits=(0, int(self.data_df.score.max()) * 1.1)
            ) +
            pn.theme_classic() +
            pn.theme(
                axis_title_x=pn.element_blank(),
                legend_title=pn.element_blank(),
                text=pn.element_text(family="Times New Roman"),
                title=pn.element_text(family="Times New Roman")
            )
        )
        p.save(self.pdf_file, width=8, height=6)
        p.save(self.png_file, width=8, height=6, dpi=300)
        return self


if __name__ == "__main__":
    (
        OverallConsensus()
        .parse()
        .stats()
        .plot()
    )
