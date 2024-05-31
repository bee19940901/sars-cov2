"""
获取二三代分型不同的基因组并建树
"""
from pathlib import Path
from argparse import ArgumentParser

import pandas as pd
from utils import Cmds
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from configs import SOFTWARE, PLATFORM3, PLATFORM2


class OverTree:

    MAFFT = SOFTWARE.joinpath("bin", "mafft")
    IQTREE2 = SOFTWARE.joinpath("bin", "iqtree")

    def __init__(self):

        ap = ArgumentParser()
        ap.add_argument("-mmf", "--mgi_merged_fas", type=str, required=True)
        ap.add_argument("-nmf", "--np_merged_fas", type=str, required=True)
        ap.add_argument("-mnc", "--mgi_nextclade_csv", type=str, required=True)
        ap.add_argument("-nnc", "--np_nextclade_csv", type=str, required=True)
        ap.add_argument("-od", "--out_dir", type=str, required=True)
        ap.add_argument("-c", "--cpu", type=int, default=16)
        ag = ap.parse_args()

        self.cpu = int(ag.cpu)
        self.mgi_merged_fas = Path(ag.mgi_merged_fas).absolute()
        self.np_merged_fas = Path(ag.np_merged_fas).absolute()
        self.mgi_nextclade_csv = Path(ag.mgi_nextclade_csv).absolute()
        self.np_nextclade_csv = Path(ag.np_nextclade_csv).absolute()
        self.out_dir = Path(ag.out_dir).absolute()
        self.sh_dir = self.out_dir.joinpath("work_sh")
        self.sh_dir.mkdir(parents=True, exist_ok=True)
        self.data_xls = self.out_dir.joinpath("data.xls")
        self.data_df = pd.DataFrame()

    @staticmethod
    def clade_func(x):
        if x[f"clade_{PLATFORM2}"] == "-" or x[f"clade_{PLATFORM3}"] == "-":
            return "-"
        if x[f"clade_{PLATFORM2}"] != "-" and x[f"clade_{PLATFORM3}"] == "-":
            return f"{PLATFORM2}_ONLY"
        if x[f"clade_{PLATFORM2}"] == "-" and x[f"clade_{PLATFORM3}"] != "-":
            return f"{PLATFORM3}_ONLY"
        if x[f"clade_{PLATFORM2}"] == x[f"clade_{PLATFORM3}"]:
            return "YES"
        return "NO"

    @staticmethod
    def pango_func(x):
        if x[f"Nextclade_pango_{PLATFORM2}"] == "-" or x[f"Nextclade_pango_{PLATFORM3}"] == "-":
            return "-"
        if x[f"Nextclade_pango_{PLATFORM2}"] != "-" and x[f"Nextclade_pango_{PLATFORM3}"] == "-":
            return f"{PLATFORM2}_ONLY"
        if x[f"Nextclade_pango_{PLATFORM2}"] == "-" and x[f"Nextclade_pango_{PLATFORM3}"] != "-":
            return f"{PLATFORM3}_ONLY"
        if x[f"Nextclade_pango_{PLATFORM2}"] == x[f"Nextclade_pango_{PLATFORM3}"]:
            return "YES"
        return "NO"

    def parse(self):
        """
        判断二三代分型结果是否一致
        :return:
        """
        # 三代数据基因组评估结果
        np_df = pd.read_csv(self.np_nextclade_csv, sep=";", usecols=[1, 2, 3])
        # 二代数据基因组评估结果
        mgi_df = pd.read_csv(self.mgi_nextclade_csv, sep=";", usecols=[1, 2, 3])
        self.data_df = (
            pd
            # 将两个结果合并
            .merge(mgi_df, np_df, on="seqName", how="outer", suffixes=(f"_{PLATFORM2}", f"_{PLATFORM3}"))
            # 使用 - 填充没有评估结果的值，也就是没有分型结果的样本
            .fillna("-")
            # 判断 clade 是否一致
            .assign(clade_judge=lambda x: x.apply(self.clade_func, axis=1))
            # 判断分型是否一致
            .assign(pango_judge=lambda x: x.apply(self.pango_func, axis=1))
            # 修改样本名
            .assign(seqName=lambda x: x["seqName"].str.replace(r"^XG", "", regex=True))
        )
        # 保存结果
        self.data_df.to_csv(self.data_xls, sep="\t", header=True, index=False)
        return self

    def stat(self, title, judge=None):
        """
        将二三代分型不同的基因组合并
        :return:
        """
        # 使用所有的序列
        if title == "all" and not judge:
            lis = (
                self.data_df
                .seqName
                .tolist()
            )
        # 使用二三代clade结果不同的序列（包含只有二代结果的序列）
        elif title == "only":
            lis = (
                self.data_df
                .query("clade_judge == @judge or clade_judge == 'NO'")
                .seqName
                .tolist()
            )
        # 使用二三代 clade 结果一致或者不一致的序列
        else:
            lis = (
                self.data_df
                .query("clade_judge == @judge")
                .seqName
                .tolist()
            )
        out_records = list()
        # 挑选序列
        with open(self.np_merged_fas, "r") as fn:
            for r in SeqIO.parse(fn, "fasta"):
                if r.id in lis:
                    out_records.append(
                        SeqRecord(seq=r.seq, id=f"{PLATFORM3}_{r.id}", description=f"{PLATFORM3}_{r.description}")
                    )
        with open(self.mgi_merged_fas, "r") as fm:
            for r in SeqIO.parse(fm, "fasta"):
                if r.id in lis:
                    out_records.append(
                        SeqRecord(seq=r.seq, id=f"{PLATFORM2}_{r.id}", description=f"{PLATFORM2}_{r.description}")
                    )
        # 保存序列
        with open(self.out_dir.joinpath(f"{title}.fas"), "w") as fw:
            SeqIO.write(out_records, fw, "fasta")

    def stats(self):
        self.stat("all")
        self.stat("diff", "NO")
        self.stat("same", "YES")
        self.stat("only", f"{PLATFORM2}_ONLY")
        return self

    def tree(self, title):
        """
        建树
        :return:
        """
        Cmds(
            [
                f"{self.MAFFT} --thread {self.cpu} {self.out_dir.joinpath(f'{title}.fas')} "
                f"> {self.out_dir.joinpath(f'{title}.aln')} && "
                f"{self.IQTREE2} -s {self.out_dir.joinpath(f'{title}.aln')} "
                f"--prefix {self.out_dir.joinpath(f'{title}_tree')} -T {self.cpu} -bb 1000 --redo\n"
            ],
            [
                self.sh_dir.joinpath(f"{title}.sh")
            ]
        ).multi_run(1)
        return self

    def trees(self):
        self.tree("diff")
        self.tree("same")
        self.tree("all")
        self.tree("only")


if __name__ == "__main__":
    (
        OverTree()
        .parse()
        .stats()
        .trees()
    )
