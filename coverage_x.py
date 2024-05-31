from argparse import ArgumentParser
from pathlib import Path
from re import sub

import pandas as pd


REF_LEN = 29903


class Coverage30X:

    def __init__(self):

        ap = ArgumentParser()
        ap.add_argument("-cd", "--coverage_dir", type=str, required=True)
        ap.add_argument("-dd", "--depth_dir", type=str, required=True)
        ap.add_argument("-od", "--out_dir", type=str, required=True)
        ap.add_argument("-md", "--min_depth", type=int, default=30)

        self.args = ap.parse_args()
        self.coverage_dir = Path(self.args.coverage_dir).absolute()
        self.depth_dir = Path(self.args.depth_dir).absolute()
        self.out_dir = Path(self.args.out_dir).absolute()
        self.out_dir.mkdir(parents=True, exist_ok=True)
        self.depth_files = list(self.depth_dir.glob("*.depth"))
        self.samples = [
            sub(r"\.depth$", "", i.name)
            for i in self.depth_files
        ]
        self.coverage_files = [
            self.coverage_dir.joinpath(f"{s}.coverage")
            for s in self.samples
        ]
        self.out_merged_file = self.out_dir.joinpath("all.coverage.merged.xls")

    def data(self):
        """
        数据处理
        :return:
        """
        out_dfs = list()
        for sample, coverage_file, depth_file in zip(self.samples, self.coverage_files, self.depth_files):
            coverage_df = pd.read_table(coverage_file)
            depth_df = pd.read_table(depth_file, header=None, names=["chrom", "pos", "depth"])
            depth_30x_df = depth_df.query("depth >= @self.args.min_depth")
            covbases = depth_30x_df.shape[0]
            coverage = covbases / REF_LEN
            mean_depth = depth_30x_df["depth"].mean()
            out_df = (
                pd
                .read_table(coverage_file)
                .assign(covbases=covbases)
                .assign(coverage=round(coverage*100, 4))
                .assign(sample=sample)
            )
            out_df.to_csv(self.out_dir.joinpath(f"{sample}.coverage"), sep="\t", header=True, index=False)
            depth_30x_df.to_csv(self.out_dir.joinpath(f"{sample}.depth"), sep="\t", header=False, index=False)
            out_dfs.append(out_df)
        (
            pd
            .concat(out_dfs)
            .sort_values(by="coverage", ascending=False)
            .to_csv(self.out_merged_file, sep="\t", header=True, index=False)
        )
        return self


if __name__ == "__main__":
    Coverage30X().data()
