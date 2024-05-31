"""
对二三代数据的突变信息进行分析
"""
import re
from argparse import ArgumentParser
from pathlib import Path
from collections import defaultdict

import pandas as pd
import plotnine as pn
import matplotlib as mpl
from scipy.stats import mannwhitneyu

from configs import PLATFORM2, PLATFORM3


class OverallMutation:

    ap = ArgumentParser()
    ap.add_argument("-nvd", "--np_vcf_dir", type=str, required=True)
    ap.add_argument("-mvd", "--mgi_vcf_dir", type=str, required=True)
    ap.add_argument("-od", "--out_dir", type=str, required=True)
    args = ap.parse_args()

    def __init__(self):

        self.np_vcf_dir = Path(self.args.np_vcf_dir).absolute()
        self.mgi_vcf_dir = Path(self.args.mgi_vcf_dir).absolute()
        self.out_dir = Path(self.args.out_dir).absolute()
        self.out_dir.mkdir(parents=True, exist_ok=True)
        self.np_good_vcfs = list(self.np_vcf_dir.glob("*.good.vcf"))
        self.np_samples = [re.sub(r"\.good\.vcf$", "", i.name) for i in self.np_good_vcfs]
        self.mgi_good_vcfs = list(self.mgi_vcf_dir.glob("*.good.vcf"))
        self.mgi_samples = [re.sub(r"\.good\.vcf$", "", i.name)for i in self.mgi_good_vcfs]
        self.stats_df = pd.DataFrame()
        self.stats_xls = self.out_dir.joinpath("stats.xls")
        self.data_df = pd.DataFrame()
        self.data_xls = self.out_dir.joinpath("data.xls")

    def __add__(self, other):
        if callable(other):
            return other(self)
        else:
            return self

    @staticmethod
    def parser(platform, good_vcfs, samples):

        def func(x):
            if len(x["ALT"]) == 1 and len(x["REF"]) == 1:
                return "snp"
            elif len(x["ALT"]) == len(x["REF"]):
                return "mnp"
            elif len(x["ALT"]) < len(x["REF"]):
                return "deletion"
            else:
                return "insertion"

        counts_dict = defaultdict(list)
        for i, s in zip(good_vcfs, samples):
            df = pd.read_csv(
                i, sep="\t", comment='#', header=None, usecols=[0, 1, 3, 4, 5, 7],
                names=["CHROM", "POS", "REF", "ALT", "QUAL", "FORMAT"])
            if df.shape[0] == 0:
                counts_dict["sample"].extend([s] * 6)
                counts_dict["platform"].extend([platform] * 6)
                counts_dict["mutation"].extend(["snp", "mnp", "deletion", "insertion", "substitution", "total"])
                counts_dict["count"].extend([0] * 6)
            else:
                df["TYPE"] = df.apply(func, axis=1)
                snp = df["TYPE"].value_counts().get("snp", 0)
                mnp = df["TYPE"].value_counts().get("mnp", 0)
                substitution = snp + mnp
                deletion = df["TYPE"].value_counts().get("deletion", 0)
                insertion = df["TYPE"].value_counts().get("insertion", 0)
                total = deletion + insertion + substitution
                counts_dict["sample"].extend([s] * 6)
                counts_dict["platform"].extend([platform]*6)
                counts_dict["mutation"].extend(["snp", "mnp", "deletion", "insertion", "substitution", "total"])
                counts_dict["count"].extend([snp, mnp, deletion, insertion, substitution, total])
        return pd.DataFrame(counts_dict)

    def parse(self):
        self.data_df = pd.concat(
            [
                self.parser(PLATFORM2, self.mgi_good_vcfs, self.mgi_samples),
                self.parser(PLATFORM3, self.np_good_vcfs, self.np_samples),
            ]
        )
        self.data_df.to_csv(self.data_xls, sep="\t", header=True, index=False)
        return self

    def stats(self):
        stats_dict = defaultdict(list)
        for i in ["total", "deletion", "insertion", "substitution", "snp", "mnp"]:
            df = self.data_df.loc[self.data_df["mutation"] == i, :]
            s, p = mannwhitneyu(
                df.loc[df["platform"] == PLATFORM2, ]["count"],
                df.loc[df["platform"] == PLATFORM3, ]["count"],
                alternative="two-sided"
            )
            stats_dict["mutation"].append(i)
            stats_dict[f"{PLATFORM2}_n"].append(df.loc[df["platform"] == PLATFORM2, ].shape[0])
            stats_dict[f"{PLATFORM3}_n"].append(df.loc[df["platform"] == PLATFORM3, ].shape[0])
            stats_dict["statistic"].append(s)
            stats_dict["p-value"].append(p)
            df.to_csv(self.out_dir.joinpath(f"{i}.xls"), sep="\t", header=True, index=False)
        self.stats_df = pd.DataFrame(stats_dict)
        self.stats_df.to_csv(self.stats_xls, sep="\t", header=True, index=False)
        return self

    @staticmethod
    def plot(data, t, png_file, pdf_file, x_order=None):
        mpl.use("cairo")
        breaks = list(range(0, data["count"].max() + 25, 25))
        p = (
            pn.ggplot(
                **{
                    "data": data,
                    "mapping": pn.aes(
                        **{"x": "mutation", "y": "count", "fill": "platform"}
                    )
                }
            ) +
            pn.geom_violin() +
            pn.geom_boxplot(
                width=0.1,
                position=pn.position_dodge(0.9),
                outlier_shape="*",
                outlier_color="orange"
            ) +
            pn.ylab(
                "Number of Mutations"
            ) +
            pn.scale_y_continuous(
                limits=(0, breaks[-1] * 1.05),
                breaks=breaks
            ) +
            pn.theme_classic() +
            pn.theme(
                axis_title_x=pn.element_blank(),
                legend_title=pn.element_blank(),
                text=pn.element_text(family="Times New Roman"),
                title=pn.element_text(family="Times New Roman")
            )
        )
        if t == "mutation":
            p += pn.scale_x_discrete(
                limits=x_order,
                labels=lambda x: [i.capitalize() for i in x]
            )
        elif t == "snp":
            p += pn.scale_x_discrete(labels=["SNP"])
        else:
            p += pn.scale_x_discrete(labels=["Total"])
        p.save(png_file, width=8, height=6, dpi=400)
        p.save(pdf_file, width=8, height=6)

    def plots(self):
        self.plot(
            self.data_df.query("mutation == 'total'"),
            'total',
            self.out_dir.joinpath("total.png"),
            self.out_dir.joinpath("total.pdf"),
        )
        self.plot(
            self.data_df.query("mutation == 'snp'"),
            'snp',
            self.out_dir.joinpath("snp.png"),
            self.out_dir.joinpath("snp.pdf"),
        )
        indels_list = ["deletion", "insertion", "substitution"]
        self.plot(
            self.data_df.query("mutation in @indels_list"),
            'mutation',
            self.out_dir.joinpath("mutation.png"),
            self.out_dir.joinpath("mutation.pdf"),
            indels_list
        )
        return self


if __name__ == "__main__":
    (
        # 初始化
        OverallMutation()
        # 输入文件解析
        + OverallMutation.parse
        # 统计分析
        + OverallMutation.stats
        # 数据可视化
        + OverallMutation.plots
    )
