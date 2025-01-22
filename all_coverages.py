from argparse import ArgumentParser
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
from re import sub

import pandas as pd

from configs import BLACKLIST, PLATFORM2, PLATFORM3, REF_LEN


# 获取所有的覆盖度
class Worker:

    def __init__(self):
        ap = ArgumentParser()
        ap.add_argument("--in_dir", type=str, required=True, help="depths/*.depth")
        ap.add_argument("--out_dir", type=str, required=True, help="coverage_minus_blacklist")
        ap.add_argument("--platform", type=str, required=True, choices=[PLATFORM2, PLATFORM3])
        ap.add_argument("--cutoff", type=int, default=30)
        ap.add_argument("--threads", type=int, default=16)
        # 命令行参数
        self.args = ap.parse_args()
        # 输入目录，包含所有的覆盖深度文件
        self.in_dir = Path(self.args.in_dir).absolute()
        # 输出文件，包含所有的覆盖度文件
        self.out_dir = Path(self.args.out_dir).absolute()
        self.out_dir.mkdir(parents=True, exist_ok=True)
        # 所有的深度文件
        self.in_depths = list(
            self.in_dir.glob("*.depth")
        )
        # 所有的样本
        self.samples = [
            sub(r"\.depth$", "", i.name)
            for i in self.in_depths
        ]
        # 所有的样本临时目录
        self.temp_dirs = [
            self.out_dir.joinpath(i)
            for i in self.samples
        ]
        # 所有样本的覆盖度文件，五列，样本名，总覆盖度，去掉黑名单位点后的覆盖度，30X覆盖度，30去掉黑名单位点后的覆盖度
        self.out_coverages = [
            self.out_dir.joinpath(s, f"{s}.coverages")
            for s in self.samples
        ]
        # 所有样本结果的汇总
        self.merge_xls = self.out_dir.joinpath(f"{self.args.platform}.coverages.merged.xls")
        # 所有的黑名单位点
        self.black_sites = list()

    def mkdirs(self):
        """
        创建结果目录
        :return:
        """
        self.out_dir.mkdir(parents=True, exist_ok=True)
        for temp_dir in self.temp_dirs:
            temp_dir.mkdir(parents=True, exist_ok=True)
        return self

    def parse_blacklist(self):
        """
        解析黑名单位点
        :return:
        """
        with open(BLACKLIST, "r") as fr:
            for line in fr:
                if line := line.strip():
                    items = line.split("\t")
                    self.black_sites += list(range(int(items[1]), int(items[2]) + 1))
        return self

    def run(self, sample: str, depth: Path, coverage: Path) -> None:
        """
        将测序深度文件转换为覆盖度文件
        :param sample: 样本名
        :param depth: 深度文件
        :param coverage: 覆盖度文件
        :return: None
        """
        total = REF_LEN
        total_minus = REF_LEN - len(self.black_sites)
        cov, cov_x, cov_minus, cov_x_minus = 0, 0, 0, 0
        with open(depth, "r") as fr, open(coverage, "w") as fw:
            fw.write("SAMPLE\tCOV\tCOV_X\tCOV_MINUS\tCOV_X_MINUS\n")
            for line in fr:
                if line := line.strip():
                    _, site, dep = line.split("\t")
                    site = int(site)
                    dep = int(dep)
                    if site > REF_LEN:
                        continue
                    if dep == 0:
                        continue
                    if dep > 0:
                        cov += 1
                        if site not in self.black_sites:
                            cov_minus += 1
                    if dep > self.args.cutoff:
                        cov_x += 1
                        if site not in self.black_sites:
                            cov_x_minus += 1
            fw.write(
                f"{sample}\t"
                f"{(cov/total)*100:.4f}\t"
                f"{(cov_x/total)*100:.4f}\t"
                f"{(cov_minus/total_minus)*100:.4f}\t"
                f"{(cov_x_minus/total_minus)*100:.4f}\n"
            )

    def multi_run(self):
        """
        多线程运行
        :return:
        """
        with ProcessPoolExecutor(self.args.threads) as pool:
            for sample, in_depth, out_coverage, in zip(self.samples, self.in_depths, self.out_coverages):
                pool.submit(self.run, sample, in_depth, out_coverage)
        return self

    def merge(self):
        """
        合并结果
        :return:
        """
        pd.concat(
            [pd.read_table(out_coverage) for out_coverage in self.out_coverages]
        ).to_csv(self.merge_xls, sep="\t", header=True, index=False)
        return self


if __name__ == "__main__":
    (
        Worker()
        .mkdirs()
        .parse_blacklist()
        .multi_run()
        .merge()
    )
