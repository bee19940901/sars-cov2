from argparse import ArgumentParser
from pathlib import Path
from dataclasses import dataclass

import pandas as pd

from configs import BLACKLIST, PLATFORM2, PLATFORM3, REF_LEN


# 用于储存分析结果
@dataclass
class Record:
    analysis: str = None
    cases_analysed: int = None
    genome_coverage: float = None
    negative_positions: int = None
    mgi_variants: int = None
    qnome_variants: int = None
    tps: int = None
    fns: int = None
    fps: int = None
    sensitivity: float = None
    precision: float = None
    jaccard_similarity: float = None
    perfect_concordance: str = None
    mgi_snvs: int = None
    qnome_snvs: int = None
    tps_snvs: int = None
    fns_snvs: int = None
    fps_snvs: int = None
    sensitivity_snvs: float = None
    precision_snvs: float = None
    jaccard_similarity_snvs: float = None
    perfect_concordance_snvs: str = None


class Pipeline:

    def __init__(self):
        """
        初始化
        """
        ap = ArgumentParser()
        ap.add_argument("--mgi_vcf_file", type=str, required=True)
        ap.add_argument("--np_vcf_file", type=str, required=True)
        ap.add_argument("--mgi_coverage_file", type=str, required=True)
        ap.add_argument("--np_coverage_file", type=str, required=True)
        ap.add_argument("--out_file", type=str, required=True)
        self.args = ap.parse_args()
        # 二代突变文件
        self.mgi_vcf_file = Path(self.args.mgi_vcf_file).absolute()
        # 三代突变文件
        self.np_vcf_file = Path(self.args.np_vcf_file).absolute()
        # 二代覆盖度信息汇总文件
        self.mgi_coverage_file = Path(self.args.mgi_coverage_file).absolute()
        # 三代覆盖度信息汇总文件
        self.np_coverage_file = Path(self.args.np_coverage_file).absolute()
        # 输出的结果文件
        self.out_file = Path(self.args.out_file).absolute()
        # 输出目录
        self.out_dir = Path(self.out_file).parent
        # 创建输出目录
        self.out_dir.mkdir(parents=True, exist_ok=True)
        # medaka 结果
        self.medaka_record = Record()
        self.medaka_record.analysis = "Medaka"
        # medaka minus blacklist 结果
        self.medaka_minus_record = Record()
        self.medaka_minus_record.analysis = "Medaka minus blacklist"
        # 总的覆盖度表
        self.all_coverage_df = pd.DataFrame()
        # 总的变异信息表
        self.all_df = pd.DataFrame()
        # 需要分析的样本
        self.samples = list()
        # 所有的黑名单位点
        self.black_sites = list()
        # 二代的变异信息
        self.mgi_vcf_df = pd.DataFrame()
        # 三代的变异信息
        self.np_vcf_df = pd.DataFrame()
        # 总的snv变异信息的表
        self.snv_df = pd.DataFrame()
        # 总的snv blacked变异信息的表
        self.snv_minus_df = pd.DataFrame()
        # 总的覆盖度表 minus black list
        self.all_minus_df = pd.DataFrame()

    @staticmethod
    def func(x: pd.Series):
        """
        辅助函数，用来判断结果是真阳、假阳还是假阴
        :param x:
        :return:
        """
        # 只有一种情况是真阳，就是 REF 和 ALT 都相等的情况
        if x["ALT_QNome"] == x["ALT_MGI"] and x["REF_QNome"] == x["REF_MGI"]:
            return "TP"
        # 如果这个位点二代有变异三代没变异则是假阴
        if pd.isna(x["ALT_QNome"]):
            return "FN"
        # 如果这个位点三代有变异二代没变异则为假阳
        if pd.isna(x["ALT_MGI"]):
            return "FP"
        # 二代和三代都有变异但不一致则为假阳
        return "FP"

    @staticmethod
    def func1(df: pd.DataFrame, record: Record) -> Record:
        """
        辅助函数，用来计算真阳数，假阳数和假阴数等结果
        :param df: 输入的数据框
        :param record: 用于更新的结果信息
        :return:
        """
        value_counts = df["JUDGE"].value_counts()
        record.tps = value_counts.get("TP", 0)
        record.fps = value_counts.get("FP", 0)
        record.fns = value_counts.get("FN", 0)
        record.sensitivity = (record.tps / (record.tps + record.fns)) * 100
        record.precision = (record.tps / (record.tps + record.fps)) * 100
        record.jaccard_similarity = (record.tps / (value_counts.sum())) * 100
        return record

    @staticmethod
    def func2(samples_n: int, df: pd.DataFrame, record: Record) -> Record:
        """
        辅助函数，检查样本二三代数据是否完全一致
        :return:
        """
        # 完全一致样本的数量
        n = 0
        for _, g in df.groupby(by="SAMPLE", as_index=False):
            if all(r["JUDGE"] == "TP" for _, r in g.iterrows()):
                n += 1
        record.perfect_concordance = f"{n}/{samples_n}"
        return record

    @staticmethod
    def func3(df: pd.DataFrame, record: Record) -> Record:
        """
        辅助函数，用来计算真阳数，假阳数和假阴数等结果（SNVs）
        :param df: 输入的数据框
        :param record: 用于更新的结果信息
        :return:
        """
        value_counts = df["JUDGE"].value_counts()
        record.tps_snvs = value_counts.get("TP", 0)
        record.fps_snvs = value_counts.get("FP", 0)
        record.fns_snvs = value_counts.get("FN", 0)
        record.sensitivity_snvs = (record.tps_snvs / (record.tps_snvs + record.fns_snvs)) * 100
        record.precision_snvs = (record.tps_snvs / (record.tps_snvs + record.fps_snvs)) * 100
        record.jaccard_similarity_snvs = (record.tps_snvs / (value_counts.sum())) * 100
        return record

    @staticmethod
    def func4(samples_n: int, df: pd.DataFrame, record: Record) -> Record:
        """
        辅助函数，样本二三代数据是否完全一致 (SNVs)
        :return:
        """
        # 完全一致样本的数量
        n = 0
        for _, g in df.groupby(by="SAMPLE", as_index=False):
            if all(r["JUDGE"] == "TP" for _, r in g.iterrows()):
                n += 1
        record.perfect_concordance_snvs = f"{n}/{samples_n}"
        return record

    def parse(self):
        """
        解析输入文件
        :return:
        """
        self.np_vcf_df = pd.read_table(self.np_vcf_file)
        # 获取三代中的样本，作为标准用于后续比较
        samples = self.np_vcf_df["SAMPLE"].unique().tolist()
        self.samples = samples
        # 筛选二代和三代都有的样本
        self.mgi_vcf_df = pd.read_table(self.mgi_vcf_file).query("SAMPLE in @samples")
        # 获取总的覆盖度表
        mgi_coverage_df = pd.read_table(self.mgi_coverage_file).assign(PLATFORM=PLATFORM2)
        np_coverage_df = pd.read_table(self.np_coverage_file).assign(PLATFORM=PLATFORM3)
        self.all_coverage_df = pd.concat([mgi_coverage_df, np_coverage_df]).query("SAMPLE in @samples")
        #  获取总的变异信息表
        self.all_df = pd.merge(
            left=self.np_vcf_df,
            right=self.mgi_vcf_df,
            on=["SAMPLE", "POS"],
            how="outer",
            suffixes=("_QNome", "_MGI")
        ).assign(
            # 使用二代结果判断三代结果是否准确
            JUDGE=lambda x: x.apply(self.func, axis=1)
        )
        # 总的 minus 表
        self.all_minus_df = pd.merge(
            left=self.np_vcf_df.query("BLACKED == 'NO'"),
            right=self.mgi_vcf_df.query("BLACKED == 'NO'"),
            on=["SAMPLE", "POS"],
            how="outer",
            suffixes=("_QNome", "_MGI")
        ).assign(
            # 使用二代结果判断三代结果是否准确
            JUDGE=lambda x: x.apply(self.func, axis=1)
        )
        # snv 总表
        self.snv_df = pd.merge(
            left=self.np_vcf_df.query("TYPE == 'SNV'"),
            right=self.mgi_vcf_df.query("TYPE == 'SNV'"),
            on=["SAMPLE", "POS"],
            how="outer",
            suffixes=("_QNome", "_MGI")
        ).assign(
            # 使用二代结果判断三代结果是否准确
            JUDGE=lambda x: x.apply(self.func, axis=1)
        )
        # 总的 snvs  minus 表
        self.snv_minus_df = pd.merge(
            left=self.np_vcf_df.query("BLACKED == 'NO' and TYPE == 'SNV'"),
            right=self.mgi_vcf_df.query("BLACKED == 'NO' and TYPE == 'SNV'"),
            on=["SAMPLE", "POS"],
            how="outer",
            suffixes=("_QNome", "_MGI")
        ).assign(
            # 使用二代结果判断三代结果是否准确
            JUDGE=lambda x: x.apply(self.func, axis=1)
        )
        # 解析黑名单文件获取黑名单位点
        with open(BLACKLIST, "r") as fr:
            for line in fr:
                if line := line.strip():
                    items = line.split("\t")
                    self.black_sites += list(range(int(items[1]), int(items[2]) + 1))
        return self

    def get_general(self):
        """
        获取样本数, 覆盖度, 阴性位点数等基本信息
        :return:
        """
        # 获取样本数
        self.medaka_record.cases_analysed = len(self.samples)
        self.medaka_minus_record.cases_analysed = len(self.samples)
        # 保存覆盖度到结果中
        self.medaka_record.genome_coverage = self.all_coverage_df["COV_X"].mean()
        self.medaka_minus_record.genome_coverage = self.all_coverage_df["COV_X_MINUS"].mean()
        # 保存阴性位点数到结果中
        self.medaka_record.negative_positions = (
            REF_LEN * len(self.samples) - self.all_df.shape[0]
        )
        self.medaka_minus_record.negative_positions = (
            (REF_LEN - len(self.black_sites)) * len(self.samples) - self.all_minus_df.shape[0]
        )
        return self

    def get_variants(self):
        """
        统计二三代所有变异位点的信息
        :return:
        """
        # 保存所有突变位点数
        self.medaka_record.mgi_variants = self.mgi_vcf_df.shape[0]
        self.medaka_minus_record.mgi_variants = self.mgi_vcf_df.query("BLACKED == 'NO'").shape[0]
        self.medaka_record.qnome_variants = self.np_vcf_df.shape[0]
        self.medaka_minus_record.qnome_variants = self.np_vcf_df.query("BLACKED == 'NO'").shape[0]
        # 保存真阳、假阳、假阴、灵敏度等结果
        self.medaka_record = self.func1(self.all_df, self.medaka_record)
        self.medaka_minus_record = self.func1(self.all_minus_df, self.medaka_minus_record)
        # 保存是否一致的结果
        self.medaka_record = self.func2(len(self.samples), self.all_df, self.medaka_record)
        self.medaka_minus_record = self.func2(len(self.samples), self.all_minus_df, self.medaka_minus_record)
        return self

    def get_snvs(self):
        """
        统计二三代snv位点的信息
        :return:
        """
        self.medaka_record.mgi_snvs = self.mgi_vcf_df.query("TYPE == 'SNV'").shape[0]
        self.medaka_minus_record.mgi_snvs = self.mgi_vcf_df.query("BLACKED == 'NO' and TYPE == 'SNV'").shape[0]
        self.medaka_record.qnome_snvs = self.np_vcf_df.query("TYPE == 'SNV'").shape[0]
        self.medaka_minus_record.qnome_snvs = self.np_vcf_df.query("BLACKED == 'NO' and TYPE == 'SNV'").shape[0]
        # 保存真阳、假阳、假阴、灵敏度等结果
        self.medaka_record = self.func3(self.snv_df, self.medaka_record)
        self.medaka_minus_record = self.func3(self.snv_minus_df, self.medaka_minus_record)
        # 保存是否一致的结果
        self.medaka_record = self.func4(len(self.samples), self.snv_df, self.medaka_record)
        self.medaka_minus_record = self.func4(len(self.samples), self.snv_minus_df, self.medaka_minus_record)
        return self

    def save(self):
        """
        保存中间文件和结果文件
        :return:
        """
        # 输出中间文件
        temp_dir = self.out_dir.joinpath("temps")
        temp_dir.mkdir(parents=True, exist_ok=True)
        self.all_df.to_csv(temp_dir.joinpath("all.variants.xls"), sep="\t", header=True, index=False)
        self.snv_df.to_csv(temp_dir.joinpath("all.snvs.xls"), sep="\t", header=True, index=False)
        self.all_minus_df.to_csv(temp_dir.joinpath("all.variants.blacked.xls"), sep="\t", header=True, index=False)
        self.snv_minus_df.to_csv(temp_dir.joinpath("all.snvs.blacked.xls"), sep="\t", header=True, index=False)
        # 输出结果文件
        out_df = (
            pd
            .DataFrame([self.medaka_record, self.medaka_minus_record])
            .set_index("analysis")
            .transpose()
        )
        out_df.index = [
            "Cases analysed",
            "Genome coverage(%)",
            "Negative_positions",
            "MGI variants",
            "QNome variants",
            "TPs",
            "FNs",
            "FPs",
            "Sensitivity(%)",
            "Precision(%)",
            "Jaccard similarity(%)",
            "Perfect concordance",
            "MGI snvs",
            "QNome snvs",
            "TPs",
            "FNs",
            "FPs",
            "Sensitivity(%)",
            "Precision(%)",
            "Jaccard similarity(%)",
            "Perfect concordance(%)"
        ]
        out_df.to_csv(self.out_file, sep="\t", header=True, index=True)
        return self


if __name__ == "__main__":
    (
        Pipeline()
        .parse()
        .get_general()
        .get_variants()
        .get_snvs()
        .save()
    )
