from argparse import ArgumentParser
from pathlib import Path
from re import sub, search
from collections import defaultdict

import pandas as pd

from configs import MINIMAP2, BWA, SAMTOOLS,  REF_FA, HOMO_FA, PLATFORM2, PLATFORM3
from utils import Cmds, Results


class Rate:

    ap = ArgumentParser()
    ap.add_argument("--raw_fq_dir", type=str, required=True)
    ap.add_argument("--out_dir", type=str, required=True)
    ap.add_argument("--platform", type=str, required=True, choices=[PLATFORM2, PLATFORM3])
    ap.add_argument("--cpu", type=int, default=1)
    args = ap.parse_args()

    def __init__(self):
        self.raw_fq_dir = Path(self.args.raw_fq_dir).absolute()
        self.out_dir = Path(self.args.out_dir).absolute()
        self.sh_dir = self.out_dir.joinpath("work_sh")
        self.sh_dir.mkdir(parents=True, exist_ok=True)
        self.raw_fq_files = list(
            self.raw_fq_dir.glob("*.fq.gz")
        )
        self.samples = [
            sub(r"\.fq\.gz$", "", i.name)
            for i in self.raw_fq_files
        ]
        self.out_home_stat_files = [
            self.out_dir.joinpath(f"{i}.home.stats")
            for i in self.samples
        ]
        self.out_sars_stat_files = [
            self.out_dir.joinpath(f"{i}.sars.stats")
            for i in self.samples
        ]
        self.sh_files = [
            self.sh_dir.joinpath(f"{i}.sh")
            for i in self.samples
        ]
        self.merged_df = pd.DataFrame()

    def __or__(self, other):
        if callable(other):
            return other(self)
        else:
            return self

    def run(self):
        Cmds(
            [
                f"{BWA} -t 16 {REF_FA} {raw_fq_file} | "
                f"{SAMTOOLS} flagstats -@ 16 - > {out_sars_stat_file} && "
                f"{BWA} -t 16 {HOMO_FA} {raw_fq_file} | "
                f"{SAMTOOLS}  flagstats -@ 16 - > {out_homo_stat_file} \n"
                if self.args.platform == PLATFORM2
                else
                f"{MINIMAP2} -ax map-ont -t 16 {REF_FA} {raw_fq_file} | "
                f"{SAMTOOLS} flagstats -@ 16 - > {out_sars_stat_file} && "
                f"{MINIMAP2} -ax map-ont -t 16 {HOMO_FA} {raw_fq_file} | "
                f"{SAMTOOLS} flagstats -@ 16 - > {out_homo_stat_file} \n"
                for raw_fq_file, out_sars_stat_file, out_homo_stat_file
                in zip(self.raw_fq_files, self.out_sars_stat_files, self.out_home_stat_files)
            ],
            self.sh_files
        ).multi_run(self.args.cpu)
        Results(
            self.out_sars_stat_files +
            self.out_home_stat_files
        ).check_exists().check_empty()
        return self

    def merge(self):
        out_merged_dict = defaultdict(list)
        for sample, sars_file, homo_file in zip(self.samples, self.out_sars_stat_files,  self.out_home_stat_files):
            out_merged_dict["SAMPLE"].append(sample)
            with open(sars_file, "r", encoding="utf-8") as fr:
                text = fr.read()
                total_reads = int(search(r"(\d+)\s+\+\s+\d+\s+primary\n", text).group(1))
                sars_mapped_reads = int(search(r"(\d+)\s+\+\s+\d+\s+primary\s+mapped\s+", text).group(1))
                out_merged_dict["SARS_TOTAL"].append(int(total_reads))
                out_merged_dict["SARS"].append(int(sars_mapped_reads))
            with open(homo_file, "r", encoding="utf-8") as fr:
                text = fr.read()
                total_reads = int(search(r"(\d+)\s+\+\s+\d+\s+primary\n", text).group(1))
                homo_mapped_reads = int(search(r"(\d+)\s+\+\s+\d+\s+primary\s+mapped\s+", text).group(1))
                out_merged_dict["HOMO_TOTAL"].append(int(total_reads))
                out_merged_dict["HOMO"].append(int(homo_mapped_reads))
        self.merged_df = (
            pd
            .DataFrame(out_merged_dict)
            .assign(
                HOMO_RATE=lambda x: x['HOMO'] / x['HOMO_TOTAL'],
                SARS_RATE=lambda x: x['SARS'] / x['SARS_TOTAL']
            )
            .sort_values(by="HOMO_RATE")

        )
        return self

    def save(self):
        self.merged_df.to_csv(
            self.out_dir.joinpath("merged.xls"),
            sep="\t",
            header=True,
            index=False
        )
        self.merged_df.filter(
            ["SAMPLE", "SARS_TOTAL", "SARS", "SARS_RATE"]
        ).to_csv(
            self.out_dir.joinpath("sars.xls"),
            sep="\t",
            header=True,
            index=False
        )
        self.merged_df.filter(
            ["SAMPLE", "HOMO_TOTAL", "HOMO", "HOMO_RATE"]
        ).to_csv(
            self.out_dir.joinpath("homo.xls"),
            sep="\t",
            header=True,
            index=False
        )
        self.merged_df.filter(
            ["SAMPLE", "HOMO_TOTAL", "SARS", "HOMO", "SARS_RATE", "HOMO_RATE"]
        ).rename(
            columns={"HOMO_TOTAL": "TOTAL"}
        ).to_csv(
            self.out_dir.joinpath("result.xls"),
            sep="\t",
            header=True,
            index=False
        )
        return self


if __name__ == "__main__":
    (
        Rate()
        # 运行命令行
        | Rate.run
        # 合并结果
        | Rate.merge
        # 保存结果
        | Rate.save
    )
